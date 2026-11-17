"""Optional local Pfam/HMMER searches. Runtime dependencies: standard library only."""

from functools import lru_cache
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time

VERSION = 1


def executable(name):
    return shutil.which(name) or next((str(path) for path in (
        Path('/opt/homebrew/bin') / name, Path('/usr/local/bin') / name)
        if path.is_file() and os.access(path, os.X_OK)), None)


def discover(root, database=None, library=None):
    """Find an installed, pressed Pfam bundle; never download or index it here."""
    configured = database or os.environ.get('CHIMERAX_PFAM_DB')
    paths = []
    if configured:
        path = Path(configured).expanduser()
        paths = [path / 'Pfam-A.hmm' if path.is_dir() else path]
    else:
        base = Path(root) / '.domain_data'
        libraries = [Path(library).expanduser()] if library else []
        settings = base / 'default-library.json'
        if not library and settings.is_file():
            value = json.loads(settings.read_text()).get('library')
            if isinstance(value, str) and value:
                libraries.append(base / value)
        for folder in libraries:
            manifest_path = folder / 'manifest.json'
            if manifest_path.is_file():
                manifest = json.loads(manifest_path.read_text())
                if manifest.get('profile_database'):
                    paths.append(folder / manifest['profile_database'])
        # Prefer the newest installed release, including an unfinished structural
        # build whose complete sequence-profile bundle is already usable.
        candidates = list(base.glob('*/profiles/Pfam*/Pfam-A.hmm'))
        def release_key(path):
            return tuple(int(p) for p in re.findall(r'\d+', path.parent.name))
        paths.extend(sorted(candidates, key=release_key, reverse=True))
    hmmscan = executable('hmmscan')
    if hmmscan is None:
        if configured:
            raise ValueError('HMMER hmmscan is not installed or available on PATH.')
        return None
    for database_path in paths:
        database_path = database_path.resolve()
        metadata_path = Path(str(database_path) + '.dat.gz')
        files = [database_path, metadata_path, *(Path(str(database_path) + suffix)
                  for suffix in ('.h3f', '.h3i', '.h3m', '.h3p'))]
        if not all(path.is_file() for path in files):
            continue
        if any(path.stat().st_mtime_ns < database_path.stat().st_mtime_ns for path in files[2:]):
            continue
        version_file = database_path.parent / 'Pfam.version.gz'
        release = None
        if version_file.is_file():
            with gzip.open(version_file, 'rt') as stream:
                match = re.search(r'Pfam release\s*:\s*([\d.]+)', stream.read())
                release = match[1] if match else None
            files.append(version_file)
        if release is None:
            match = re.fullmatch(r'Pfam([\d.]+)', database_path.parent.name)
            release = match[1] if match else 'unknown release'
        programs = {'hmmscan': hmmscan, 'hmmfetch': executable('hmmfetch'), 'hmmpress': executable('hmmpress')}
        signature = [VERSION, programs, [(str(path), path.stat().st_size, path.stat().st_mtime_ns)
                     for path in [*files, *(Path(p) for p in programs.values() if p)]]]
        identity = 'local-' + hashlib.sha256(json.dumps(signature).encode()).hexdigest()
        return {'kind': 'local', 'identity': identity, 'database': str(database_path),
                'metadata': str(metadata_path), 'hmmscan': hmmscan, 'release': release,
                'hmmfetch': programs['hmmfetch'], 'hmmpress': programs['hmmpress']}
    if configured:
        raise ValueError('The configured Pfam database needs matching .dat.gz metadata and current hmmpress indexes.')
    return None


@lru_cache(maxsize=4)
def _metadata(path, size, modified):
    families, current = {}, {}
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#=GF '):
                key, value = line[5:].split(None, 1)
                if key == 'NE':
                    current.setdefault('NE', []).append(value.strip())
                else:
                    current[key] = value.strip()
            elif line.strip() == '//':
                accession = current['AC'].split('.')[0]
                ga = [float(v.strip(';')) for v in current['GA'].split()]
                if (not re.fullmatch(r'PF\d{5,}', accession) or accession in families or len(ga) != 2
                        or not all(math.isfinite(v) for v in ga)):
                    raise ValueError('Invalid Pfam family metadata.')
                families[accession] = {'accession': accession, 'version': current['AC'], 'ga': ga,
                    'name': current['ID'], 'description': current.get('DE', ''), 'type': current['TP'],
                    'length': int(current['ML']), 'clan': current.get('CL'), 'nested': current.get('NE', [])}
                current = {}
    if current or not families:
        raise ValueError('Incomplete Pfam family metadata.')
    return families


def metadata(backend):
    path = Path(backend['metadata'])
    stat = path.stat()
    return _metadata(str(path), stat.st_size, stat.st_mtime_ns)


def _run(arguments, output, deadline, cancel, progress):
    if cancel.is_set():
        raise InterruptedError('Pfam search cancelled.')
    if time.monotonic() >= deadline:
        raise TimeoutError('Local Pfam search timed out.')
    with output.open('w') as log:
        process = subprocess.Popen([str(a) for a in arguments], stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=os.name != 'nt')
        last_report = time.monotonic()
        try:
            while process.poll() is None:
                if cancel.is_set():
                    raise InterruptedError('Pfam search cancelled.')
                if time.monotonic() >= deadline:
                    raise TimeoutError('Local Pfam search timed out.')
                if time.monotonic() - last_report >= 1:
                    progress()
                    last_report = time.monotonic()
                cancel.wait(.05)
        finally:
            if process.poll() is None:
                try:
                    if os.name == 'nt':
                        process.terminate()
                    else:
                        os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    try:
                        if os.name == 'nt':
                            process.kill()
                        else:
                            os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
    if process.returncode:
        raise ValueError('Local HMMER failed: ' + output.read_text(errors='replace')[-1500:].strip())


def parse_table(path, families, queries, release):
    hits = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        row = line.split(maxsplit=22)
        if len(row) < 22:
            raise ValueError('Truncated local HMMER domain table.')
        accession = row[1].split('.')[0]
        family, query = families[accession], queries[row[3]]
        if row[1] != family['version'] or int(row[2]) != family['length'] or int(row[5]) != len(query['sequence']):
            raise ValueError('Local Pfam profiles and metadata do not match.')
        first, last = int(row[17]), int(row[18])
        hf, hl = int(row[15]), int(row[16])
        seq_score, score, evalue, bias = (float(row[i]) for i in (7, 13, 12, 14))
        if (not 1 <= first <= last <= len(query['sequence']) or not 1 <= hf <= hl <= family['length']
                or not all(math.isfinite(v) for v in (seq_score, score, evalue, bias)) or evalue < 0):
            raise ValueError('Invalid local Pfam alignment coordinates or scores.')
        offset = query['start'] - 1
        hits.append({'origin': 'pfam', 'source': 'Pfam local (' + release + ')', 'database': 'pfam',
                     'accession': accession, 'name': family['name'], 'type': family['type'],
                     'description': family['description'], 'start': first + offset, 'end': last + offset,
                     'evalue': evalue, 'bitscore': score, 'outcompeted': False,
                     'tentative': seq_score < family['ga'][0] or score < family['ga'][1],
                     'gathering_cutoff': family['ga'][1], 'sequence_gathering_cutoff': family['ga'][0],
                     'sequence_score': seq_score, 'profile_coverage': (hl - hf + 1) / family['length'],
                     'composition_bias': bias, 'pfam_backend': 'local', 'pfam_release': release,
                     'query_region': [query['start'], query['start'] + len(query['sequence']) - 1],
                     '_owner': query['owner']})
    return hits


def resolve(hits, families):
    """Deduplicate retries, then apply same-clan competition with nested exceptions."""
    selected = []
    for original in sorted(hits, key=lambda h: (h['tentative'], h['evalue'], -h['bitscore'])):
        hit = {**original, 'outcompeted': False}
        def overlap(other):
            return max(0, min(hit['end'], other['end']) - max(hit['start'], other['start']) + 1)
        if any(old['accession'] == hit['accession'] and overlap(old) >=
               .5 * min(old['end'] - old['start'] + 1, hit['end'] - hit['start'] + 1) for old in selected):
            continue
        family = families[hit['accession']]
        for old in selected:
            previous = families[old['accession']]
            if (not old['outcompeted'] and old['accession'] != hit['accession'] and family['clan']
                    and family['clan'] == previous['clan'] and overlap(old)
                    and old['name'] not in family['nested'] and hit['name'] not in previous['nested']):
                hit['outcompeted'] = True
                break
        selected.append(hit)
    return selected


def retry_queries(sequence, hits, owner):
    curated = {p for h in hits if not h['tentative'] and not h['outcompeted']
               for p in range(h['start'], h['end'] + 1)}
    regions = set()
    for hit in hits:
        if hit['tentative'] and not hit['outcompeted']:
            regions.add((max(1, hit['start'] - 12), min(len(sequence), hit['end'] + 12)))
    first = None
    for p in range(1, len(sequence) + 2):
        if p <= len(sequence) and p not in curated:
            if first is None:
                first = p
        elif first is not None:
            if p - first >= 20:
                regions.add((max(1, first - 12), min(len(sequence), p + 11)))
            first = None
    return [{'sequence': sequence[a - 1:b], 'start': a, 'owner': owner} for a, b in sorted(regions)
            if (a, b) != (1, len(sequence))]


def search(sequences, backend, timeout, cancel, report):
    """Curated full scan, reported tentative hits, and ordinary targeted retries."""
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Local Pfam timeout must be finite and positive.')
    if cancel.is_set():
        raise InterruptedError('Pfam search cancelled.')
    deadline = time.monotonic() + timeout
    sequences = list(dict.fromkeys(sequences))
    if not sequences:
        return {}
    families = metadata(backend)
    full = {f'chain{i}': {'sequence': seq, 'start': 1, 'owner': i} for i, seq in enumerate(sequences)}
    release = backend['release']
    with tempfile.TemporaryDirectory(prefix='trimmings-pfam-') as temporary:
        folder = Path(temporary)
        def scan(database, queries, name, curated, expected_models):
            report('Local Pfam ' + release + ': ' + name)
            fasta, table, log = (folder / (name + suffix) for suffix in ('.fasta', '.tsv', '.log'))
            fasta.write_text(''.join(f'>{key}\n{q["sequence"]}\n' for key, q in queries.items()))
            options = ['--cut_ga'] if curated else ['-E', '10', '--domE', '10']
            _run([backend['hmmscan'], '--noali', '--cpu', '4', '-Z', len(families), '--domtblout', table,
                  *options, database, fasta], log, deadline, cancel, lambda: report('Local Pfam ' + release + ': ' + name))
            counts = re.findall(r'Target model\(s\):\s+(\d+)', log.read_text())
            if len(counts) != len(queries) or any(int(n) != expected_models for n in counts):
                raise ValueError('Local Pfam database inventory does not match its metadata.')
            hits = parse_table(table, families, queries, release)
            if curated and any(h['tentative'] for h in hits):
                raise ValueError('A local curated match failed its recorded gathering thresholds.')
            return hits
        all_hits = scan(backend['database'], full, 'curated scan', True, len(families))
        reported = scan(backend['database'], full, 'tentative scan', False, len(families))
        all_hits.extend(h for h in reported if h['tentative'] and h['evalue'] <= .01)
        if backend.get('hmmfetch') and backend.get('hmmpress'):
            candidates_by_owner = {owner: {h['accession'] for h in all_hits if h['_owner'] == owner}
                                   for owner in range(len(sequences))}
            candidates = set().union(*candidates_by_owner.values())
            regions = {}
            for owner, seq in enumerate(sequences):
                hits = resolve([h for h in all_hits if h['_owner'] == owner], families)
                for query in retry_queries(seq, hits, owner):
                    regions[f'region{len(regions)}'] = query
            if candidates and regions:
                keys, subset = folder / 'keys', folder / 'candidates.hmm'
                keys.write_text(''.join(families[a]['version'] + '\n' for a in sorted(candidates)))
                report('Local Pfam ' + release + ': preparing region retries')
                _run([backend['hmmfetch'], '-o', subset, '-f', backend['database'], keys], folder / 'fetch.log',
                     deadline, cancel, lambda: report('Preparing local Pfam profiles'))
                _run([backend['hmmpress'], '-f', subset], folder / 'press.log', deadline, cancel,
                     lambda: report('Indexing local Pfam profiles'))
                retries = scan(subset, regions, 'padded region retries', True, len(candidates))
                # Batch profile extraction/scanning without letting one protein's
                # candidate families change annotations of an unrelated protein.
                all_hits.extend(h for h in retries if h['accession'] in candidates_by_owner[h['_owner']])
        results = {}
        for owner, seq in enumerate(sequences):
            results[seq] = [{key: value for key, value in h.items() if key != '_owner'}
                            for h in resolve([h for h in all_hits if h['_owner'] == owner], families)]
        return results
