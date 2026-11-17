"""Build and benchmark a small Foldseek library labeled by CDD family membership.

Building/extracting mmCIF files requires gemmi. Searching an existing library uses
only the standard library and the Foldseek executable.
"""

import argparse
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
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen


# Deliberately small, diverse pilot, not a comprehensive annotation database.
PILOT_MODELS = ('pfam02026', 'cd00009', 'cd00876', 'pfam00042', 'cd01659',
                'cd00051', 'cd00174', 'cd00180', 'cd00200', 'cd00590', 'cd00096',
                'cd00657', 'cd11709', 'cd12877', 'cd23278', 'pfam00515', 'pfam00023')
FIELDS = ('query', 'target', 'evalue', 'bits', 'qstart', 'qend', 'qlen',
          'tstart', 'tend', 'tlen', 'qcov', 'tcov', 'qtmscore', 'ttmscore',
          'alntmscore', 'qaln', 'taln')
SCREEN_FIELDS = tuple(field for field in FIELDS if field not in ('qtmscore', 'ttmscore'))
ALIGNMENT_FIELDS = tuple(field for field in SCREEN_FIELDS if field != 'alntmscore')
DEFAULT_THREADS = 4
INTERPRO = 'https://www.ebi.ac.uk/interpro/api'
SEARCH_VERSION = 12


def fetch(url, cache, maximum=150_000_000):
    """Cache completed public downloads; never replace good data with a failed request."""
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / hashlib.sha256(url.encode()).hexdigest()
    if path.exists():
        return path.read_bytes()
    for attempt in range(3):
        try:
            request = Request(url, headers={'User-Agent': 'ChimeraX-Trimmings/1.0'})
            with urlopen(request, timeout=45) as response:
                if response.status == 204:
                    return b''
                data = response.read(maximum + 1)
            if len(data) > maximum:
                raise ValueError('Public download exceeds the size limit: ' + url)
            with tempfile.NamedTemporaryFile(dir=cache, delete=False) as stream:
                stream.write(data)
                temporary = Path(stream.name)
            temporary.replace(path)
            time.sleep(.15)
            return data
        except HTTPError as error:
            if error.code == 404:
                return b''
            if error.code not in (429, 500, 502, 503, 504) or attempt == 2:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError('Download failed: ' + url)


def fetch_json(url, cache):
    data = fetch(url, cache)
    return json.loads(data) if data else None


def cdd_descriptions(cache):
    data = fetch('https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz', cache)
    rows = {}
    for line in gzip.decompress(data).decode().splitlines():
        fields = line.split('\t')
        if len(fields) >= 5:
            rows[fields[1]] = {'pssm_id': fields[0], 'name': fields[2], 'description': fields[3]}
    return rows


def api_model(accession):
    if re.fullmatch(r'pfam[0-9]{5}', accession):
        return 'pfam', 'PF' + accession[4:]
    if re.fullmatch(r'cd[0-9]{5}', accession):
        return 'cdd', accession
    raise ValueError('Expected a CDD accession such as cd00009 or pfam02026.')


def residue_letter(residue):
    import gemmi

    return gemmi.find_tabulated_residue(residue.name).one_letter_code.upper()


def write_ca(path, residues):
    """One chain, with an explicit sequential index; callers retain the source map."""
    if len(residues) > 9999:
        raise ValueError('A PDB query cannot contain more than 9,999 residues.')
    lines = []
    for index, residue in enumerate(residues, 1):
        atom = residue.find_atom('CA', '*')
        p = atom.pos
        if not all(math.isfinite(v) and -999 <= v <= 9999 for v in (p.x, p.y, p.z)):
            raise ValueError('Coordinates cannot be represented safely as a PDB query.')
        lines.append(f'ATOM  {index:5d}  CA  {residue.name:>3s} A{index:4d}    '
                     f'{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}{1.:6.2f}{atom.b_iso:6.2f}           C\n')
    Path(path).write_text(''.join(lines) + 'TER\nEND\n')


def crop_reference(cif, entry, location):
    import gemmi

    structure = gemmi.read_structure(str(cif))
    fragments = location['fragments']
    selected = {p for fragment in fragments for p in range(fragment['start'], fragment['end'] + 1)}
    chain_id = entry['chain']
    # InterPro structure coordinates use label_seq_id, not author residue numbers.
    groups = {}
    for chain in structure[0]:
        for residue in chain:
            if residue.label_seq in selected and residue.find_atom('CA', '*'):
                groups.setdefault((chain.name, residue.subchain), []).append(residue)
    choices = [rs for (author, label), rs in groups.items() if label == chain_id]
    if not choices:
        choices = [rs for (author, label), rs in groups.items() if author == chain_id]
    valid = []
    for residues in choices:
        residues.sort(key=lambda r: r.label_seq)
        if len(residues) < 25 or len(residues) / len(selected) < .8:
            continue
        sequence = entry['sequence']
        if any(r.label_seq > len(sequence) or (sequence[r.label_seq - 1] not in ('X', residue_letter(r))
                                              and residue_letter(r) != 'X') for r in residues):
            continue
        valid.append(residues)
    if len(valid) != 1:
        raise ValueError('Cannot uniquely verify the PDB chain, domain coordinates, and sequence.')
    return valid[0]


def build_library(output, models=PILOT_MODELS, representatives=3, executable='foldseek'):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    cache, structures = output / 'downloads', output / 'domains'
    structures.mkdir(exist_ok=True)
    descriptions = cdd_descriptions(cache)
    records, failures = [], []
    for accession in models:
        print('Building', accession, flush=True)
        database, member = api_model(accession)
        label = descriptions.get(accession)
        if label is None:
            raise ValueError('Model missing from the current CDD annotation table: ' + accession)
        metadata = fetch_json(f'{INTERPRO}/entry/{database}/{member}/', cache)
        if not metadata:
            failures.append({'accession': accession, 'reason': 'No InterPro member metadata'})
            continue
        data = fetch_json(f'{INTERPRO}/structure/pdb/entry/{database}/{member}/?page_size=10', cache) or {}
        candidates = data.get('results', [])
        representative = (metadata.get('metadata', {}).get('representative_structure') or {}).get('accession')
        if representative and not any(r['metadata']['accession'] == representative for r in candidates):
            single = fetch_json(f'{INTERPRO}/structure/pdb/{representative}/entry/{database}/{member}/', cache)
            if single:
                candidates.insert(0, single)
        candidates.sort(key=lambda r: (r['metadata']['accession'] != representative,
                                        min((e.get('sequence_length', 99999) for e in r.get('entries', [])), default=99999)))
        seen_sequences, accepted = set(), 0
        for candidate in candidates:
            if accepted >= representatives:
                break
            pdb = candidate['metadata']['accession']
            if not re.fullmatch(r'[A-Za-z0-9]{4}', pdb):
                continue
            for entry in candidate.get('entries', []):
                if accepted >= representatives:
                    break
                if entry.get('sequence_length', 99999) > 2000:
                    continue  # Keep pilot downloads modest; larger parents can be added later.
                locations = entry.get('entry_structure_locations', [])
                for occurrence, location in enumerate(locations, 1):
                    if accepted >= representatives:
                        break
                    try:
                        cif_data = fetch(f'https://files.rcsb.org/download/{pdb}.cif', cache)
                        cif = output / 'source.cif'
                        cif.write_bytes(cif_data)
                        residues = crop_reference(cif, entry, location)
                        sequence = ''.join(residue_letter(r) for r in residues)
                        if sequence in seen_sequences:
                            continue
                        seen_sequences.add(sequence)
                        identifier = re.sub(r'[^A-Za-z0-9_]', '_', f'{accession}_{pdb}_{entry["chain"]}_{occurrence}')
                        write_ca(structures / (identifier + '.pdb'), residues)
                        records.append({'id': identifier, 'accession': member, 'cdd_accession': accession,
                                        'database': database, **label, 'pdb': pdb, 'chain': entry['chain'],
                                        'segments': [(f['start'], f['end']) for f in location['fragments']],
                                        'residue_map': [[r.label_seq, r.seqid.num, str(r.seqid.icode).strip()] for r in residues],
                                        'sequence': sequence, 'structure_title': candidate['metadata'].get('name', ''),
                                        'evidence': 'InterPro member-database match on an experimental PDB chain',
                                        'entry_url': f'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={label["pssm_id"]}'})
                        accepted += 1
                    except (ValueError, KeyError, IndexError, RuntimeError) as error:
                        failures.append({'accession': accession, 'pdb': pdb, 'reason': str(error)})
        print('  accepted', accepted, 'reference domains', flush=True)
        if not accepted:
            failures.append({'accession': accession, 'reason': 'No verifiable compact reference domains in pilot candidates'})
    (output / 'source.cif').unlink(missing_ok=True)
    if not records:
        raise ValueError('No verified reference domains were obtained.')
    manifest = {'schema': 1, 'pilot': True, 'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                'description': 'Small experimental CDD-labeled library; includes Pfam models imported by CDD.',
                'requested_models': list(models), 'records': records, 'failures': failures,
                'database': 'foldseek_db', 'foldseek_version': subprocess.check_output([executable, 'version'], text=True).strip()}
    # Build from an explicit list so rerunning with fewer families cannot include stale files.
    with tempfile.TemporaryDirectory(dir=output) as directory:
        inputs = Path(directory) / 'inputs'
        inputs.mkdir()
        for record in records:
            shutil.copy2(structures / (record['id'] + '.pdb'), inputs)
        subprocess.run([executable, 'createdb', str(inputs), str(output / 'foldseek_db'), '--threads', '2'], check=True,
                       stdout=subprocess.DEVNULL)
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def parse_segments(text):
    if not re.fullmatch(r'\d+-\d+(?:_\d+-\d+)*(?:,\d+-\d+(?:_\d+-\d+)*)*', text):
        raise ValueError('Use comma-separated domains, with underscore-separated pieces for discontinuous domains.')
    domains = [[tuple(map(int, pair.split('-'))) for pair in domain.split('_')] for domain in text.split(',')]
    if any(a < 1 or b < a for domain in domains for a, b in domain):
        raise ValueError('Domain ranges must be positive and increasing.')
    if any(b >= c for domain in domains for (_, b), (c, _) in zip(domain, domain[1:])):
        raise ValueError('Pieces within one domain must be ordered and non-overlapping.')
    return domains


def query_domains(path, output, segments=None, chain_id=None):
    import gemmi

    structure = gemmi.read_structure(str(path))
    chains = [c for c in structure[0] if any(r.find_atom('CA', '*') for r in c)]
    if chain_id:
        chains = [c for c in chains if c.name == chain_id]
    if len(chains) != 1:
        raise ValueError('Specify one protein chain with --chain.')
    chain = chains[0]
    residues = [r for r in chain if r.find_atom('CA', '*') and r.label_seq is not None]
    if not residues:
        raise ValueError('Query requires mmCIF polymer sequence positions (label_seq_id).')
    domains = parse_segments(segments) if segments else [[(1, max(r.label_seq for r in residues))]]
    queries = {}
    for index, pieces in enumerate(domains, 1):
        selected = [r for r in residues if any(a <= r.label_seq <= b for a, b in pieces)]
        if len(selected) < 20:
            raise ValueError('Each query domain needs at least 20 modeled C-alpha atoms.')
        identifier = f'domain{index}'
        write_ca(Path(output) / (identifier + '.pdb'), selected)
        queries[identifier] = {'segments': pieces, 'positions': [r.label_seq for r in selected],
                               'author_residues': [[chain.name, r.seqid.num, str(r.seqid.icode).strip()] for r in selected]}
    return queries


def parse_hits(text, queries, manifest, fields=FIELDS):
    targets = {record['id']: record for record in manifest['records']}
    hits = []
    for line in text.splitlines():
        values = line.split('\t')
        if len(values) != len(fields):
            raise ValueError('Invalid Foldseek output columns.')
        row = dict(zip(fields, values))
        # Foldseek releases differ in whether the PDB suffix is retained.
        query, target = row['query'].removesuffix('.pdb'), row['target'].removesuffix('.pdb')
        if query not in queries or target not in targets:
            raise ValueError('Foldseek returned an unknown query or target identifier.')
        for field in ('evalue', 'bits', 'qcov', 'tcov', 'qtmscore', 'ttmscore', 'alntmscore'):
            if field not in fields and field in ('qtmscore', 'ttmscore'):
                row[field] = None  # Display-only scores are calculated for accepted candidates below.
                continue
            row[field] = float(row[field])
            if not math.isfinite(row[field]) or row[field] < 0:
                raise ValueError('Invalid Foldseek score.')
        for field in ('qcov', 'tcov', 'qtmscore', 'ttmscore'):
            if row[field] is not None and row[field] > 1:
                raise ValueError(f'Invalid normalized Foldseek score: {field}={row[field]:g} '
                                 f'for {queries[query].get("label", query)} against {target}.')
        for field in ('qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen'):
            row[field] = int(row[field])
        target_length = targets[target].get('length')
        if target_length is None:
            target_length = len(targets[target]['sequence'])
        if (row['qlen'] != len(queries[query]['positions']) or not 1 <= row['qstart'] <= row['qend'] <= row['qlen']
                or row['tlen'] != target_length or not 1 <= row['tstart'] <= row['tend'] <= row['tlen']
                or len(row['qaln']) != len(row['taln'])
                or len(row['qaln'].replace('-', '')) != row['qend'] - row['qstart'] + 1
                or len(row['taln'].replace('-', '')) != row['tend'] - row['tstart'] + 1):
            raise ValueError('Invalid Foldseek alignment coordinates.')
        q, t, aligned = row['qstart'] - 1, row['tstart'] - 1, []
        for aa, bb in zip(row['qaln'], row['taln']):
            if aa != '-' and bb != '-':
                aligned.append(queries[query]['positions'][q])
            q += aa != '-'
            t += bb != '-'
        if q != row['qend'] or t != row['tend']:
            raise ValueError('Foldseek alignment strings do not match their endpoints.')
        if row['alntmscore'] > 1:
            # structureconvertalis passes min(end-start), excluding one residue,
            # to computeTMscore, even with --exact-tmscore. Only correct the
            # identifiable case: every residue of the shorter span is paired
            # and the score is inside that off-by-one bound (4-digit output).
            span = min(row['qend'] - row['qstart'] + 1, row['tend'] - row['tstart'] + 1)
            if span <= 1 or len(aligned) != span or row['alntmscore'] > span / (span - 1) + .0005:
                raise ValueError(f'Invalid normalized Foldseek score: alntmscore={row["alntmscore"]:g} '
                                 f'for {queries[query].get("label", query)} against {target}.')
            row['alntmscore_raw'] = row['alntmscore']
            row['alntmscore'] = min(1., row['alntmscore'] * (span - 1) / span)
        hits.append({**row, 'query': query, 'target': target, 'reference': targets[target], 'aligned_positions': aligned})
    return sorted(hits, key=lambda hit: (hit['evalue'], -hit['bits']))


def qualifying_hit(hit):
    return (hit['evalue'] <= .01 and hit['tcov'] >= .5 and hit['alntmscore'] >= .5
            and len(hit['aligned_positions']) >= 20)


def select_matches(hits, per_family=True):
    """Keep distinct occurrences; collapse competing representatives of one occurrence."""
    selected, footprints = [], {}
    for hit in sorted(hits, key=lambda h: (h['evalue'], -h['bits'])):
        if not qualifying_hit(hit):
            continue
        key = hit['query'], hit['reference']['accession'] if per_family else None
        positions = set(hit['aligned_positions'])
        if any(len(positions & previous) / min(len(positions), len(previous)) >= .5
               for previous in footprints.get(key, [])):
            continue
        footprints.setdefault(key, []).append(positions)
        selected.append(hit)
    return selected


def library_reference_count(manifest):
    return manifest['reference_count'] if manifest.get('schema') == 2 else len(manifest['records'])


def library_file(library, name):
    path = (Path(library) / name).resolve()
    if not path.is_relative_to(Path(library).resolve()):
        raise ValueError('Reference-library files must be inside the library directory.')
    return path


def hit_metadata(library, manifest, text, fields=FIELDS):
    if manifest['schema'] == 1:
        return manifest
    import sqlite3
    from contextlib import closing

    targets = set()
    for line in text.splitlines():
        values = line.split('\t')
        if len(values) != len(fields):
            raise ValueError('Invalid Foldseek output columns.')
        targets.add(values[1].removesuffix('.pdb'))
    path = library_file(library, manifest['reference_metadata'])
    records = []
    # Load only returned targets, rather than the full collection's metadata.
    with closing(sqlite3.connect(f'file:{path}?mode=ro', uri=True)) as connection:
        identifiers = sorted(targets)
        for first in range(0, len(identifiers), 500):
            batch = identifiers[first:first + 500]
            rows = connection.execute('SELECT record FROM refs WHERE id IN (' + ','.join('?' for _ in batch) + ')', batch)
            records.extend(json.loads(row[0]) for row in rows)
    return {**manifest, 'records': records}


def _run_foldseek(command, log_path, timeout, cancel=None, progress=None):
    """Run one stage with cancellation, progress and a bounded lifetime."""
    if timeout <= 0:
        raise TimeoutError('Foldseek search timed out.')
    if cancel is not None and cancel.is_set():
        raise InterruptedError('Foldseek search cancelled.')
    started = time.monotonic()
    next_progress = started + 1
    with Path(log_path).open('a+') as log:
        process = subprocess.Popen(list(map(str, command)), stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while process.poll() is None:
                if cancel is not None and cancel.is_set():
                    raise InterruptedError('Foldseek search cancelled.')
                if time.monotonic() - started > timeout:
                    raise TimeoutError('Foldseek search timed out.')
                if progress is not None and time.monotonic() >= next_progress:
                    progress()
                    next_progress = time.monotonic() + 1
                time.sleep(.05)
            if process.returncode:
                log.seek(0)
                raise RuntimeError('Foldseek failed: ' + log.read()[-2500:])
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()


def _parallel_scoring(jobs, directory, deadline, cancel, progress):
    """Score independent alignment shards without exceeding the CPU budget."""
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    from threading import Event
    from types import SimpleNamespace

    stopping = Event()
    token = SimpleNamespace(is_set=lambda: stopping.is_set() or (cancel is not None and cancel.is_set()))
    def worker(index, commands):
        for command in commands:
            _run_foldseek(command, directory / f'score-{index}.log', deadline - time.monotonic(), token)
    pool = ThreadPoolExecutor(max_workers=len(jobs))
    pending = {pool.submit(worker, index, commands) for index, commands in enumerate(jobs)}
    next_progress = time.monotonic() + 1
    try:
        while pending:
            if token.is_set():
                raise InterruptedError('Foldseek search cancelled.')
            if time.monotonic() >= deadline:
                raise TimeoutError('Foldseek search timed out.')
            ready, pending = wait(pending, timeout=.05, return_when=FIRST_COMPLETED)
            for future in ready:
                future.result()
            if progress is not None and time.monotonic() >= next_progress:
                progress()
                next_progress = time.monotonic() + 1
    finally:
        stopping.set()
        for future in pending:
            future.cancel()
        pool.shutdown(wait=True, cancel_futures=True)


def _screen_alignments(query_db, database, alignment_db, result, directory, threads,
                       deadline, cancel, progress, executable, convert):
    # Exact TM scoring is parallelized by query inside Foldseek. Split a large
    # result into target shards so one long query cannot leave most cores idle.
    index = Path(str(alignment_db) + '.index')
    size = sum(int(line.split()[2]) for line in index.read_text().splitlines()) if index.is_file() else 0
    workers = min(threads, 8, max(1, size // 16384))
    if workers == 1:
        convert(alignment_db, result, SCREEN_FIELDS)
        return result.read_text()
    ordering = directory / 'ordering.tsv'
    convert(alignment_db, ordering, ALIGNMENT_FIELDS)
    rows = [line.split('\t') for line in ordering.read_text().splitlines()]
    if any(len(row) != len(ALIGNMENT_FIELDS) for row in rows):
        raise ValueError('Invalid Foldseek alignment columns.')
    order = [tuple(row) for row in rows]
    if len(set(order)) != len(order):
        # Unusual custom databases can contain duplicate records. Preserve their
        # established conversion behavior rather than guessing at record identity.
        convert(alignment_db, result, SCREEN_FIELDS)
        return result.read_text()
    jobs, tables = [], []
    for shard in range(workers):
        source, table = directory / f'shard-{shard}', directory / f'shard-{shard}.tsv'
        tables.append(table)
        jobs.append([
            [executable, 'filterdb', alignment_db, source, '--filter-expression',
             f'$1 % {workers} == {shard}', '--threads', '1'],
            [executable, 'convertalis', query_db, database, source, table, '--threads', '1',
             '--format-output', ','.join(SCREEN_FIELDS), '--exact-tmscore', '1'],
        ])
    _parallel_scoring(jobs, directory, deadline, cancel, progress)
    scored = {}
    for table in tables:
        for line in table.read_text().splitlines():
            values = line.split('\t')
            if len(values) != len(SCREEN_FIELDS):
                raise ValueError('Invalid Foldseek score columns.')
            fields = dict(zip(SCREEN_FIELDS, values))
            key = tuple(fields[field] for field in ALIGNMENT_FIELDS)
            if key in scored:
                raise ValueError('Foldseek scored an alignment more than once.')
            scored[key] = line
    if set(scored) != set(order):
        raise ValueError('Foldseek changed or omitted alignments during parallel scoring.')
    # Keep the original order, including exact score ties that guide retries.
    return ''.join(scored[key] + '\n' for key in order)


def _merge_alignment_shards(shards, destination, prefilter, deadline, cancel=None):
    """Restore structurealign's global score/length/target ordering per query."""
    expected = [int(line.split()[0]) for line in Path(str(prefilter) + '.index').read_text().splitlines()]
    rows = {key: [] for key in expected}
    seen, dbtype = set(), None
    for shard in shards:
        kind = Path(str(shard) + '.dbtype').read_bytes()
        if dbtype is not None and kind != dbtype:
            raise ValueError('Candidate alignment shards have different database types.')
        dbtype = kind
        with Path(shard).open('rb') as stream:
            for line in Path(str(shard) + '.index').read_text().splitlines():
                if cancel is not None and cancel.is_set():
                    raise InterruptedError('Foldseek search cancelled.')
                if time.monotonic() >= deadline:
                    raise TimeoutError('Foldseek search timed out.')
                key, offset, length = map(int, line.split())
                if key not in rows or offset < 0 or length < 1:
                    raise ValueError('Invalid candidate alignment shard index.')
                stream.seek(offset)
                data = stream.read(length)
                if len(data) != length or not data.endswith(b'\0'):
                    raise ValueError('Truncated candidate alignment shard.')
                for record in data[:-1].splitlines(keepends=True):
                    fields = record.split()
                    if len(fields) < 10:
                        raise ValueError('Invalid candidate alignment record.')
                    target, score, target_length = int(fields[0]), int(fields[1]), int(fields[9])
                    if (key, target) in seen:
                        raise ValueError('Candidate alignment was calculated more than once.')
                    seen.add((key, target))
                    rows[key].append(((-score, target_length, target), record))
    with Path(destination).open('wb') as data, Path(str(destination) + '.index').open('w') as index:
        for key in expected:
            offset = data.tell()
            for _, record in sorted(rows[key], key=lambda row: row[0]):
                data.write(record)
            data.write(b'\0')
            index.write(f'{key}\t{offset}\t{data.tell() - offset}\n')
    if dbtype is None:
        raise ValueError('No candidate alignment shards were produced.')
    Path(str(destination) + '.dbtype').write_bytes(dbtype)


def _candidate_parallel_alignment(query_db, database, alignment_db, directory, threads,
                                  deadline, cancel, progress, executable, max_seqs, evalue):
    """Keep the full prefilter, then independently align target-candidate shards."""
    directory.mkdir()
    prefilter = directory / 'prefilter'
    _run_foldseek([executable, 'ungappedprefilter', str(query_db) + '_ss', str(database) + '_ss', prefilter,
                  '--threads', threads, '--max-seqs', max_seqs, '--comp-bias-corr', '1',
                  '--comp-bias-corr-scale', '0.15', '--min-ungapped-score', '30',
                  '--sub-mat', 'aa:3di.out,nucl:3di.out', '-c', '0', '--cov-mode', '0',
                  '-e', '1.7976931348623157e+308', '--prefilter-mode', '1'],
                  directory / 'prefilter.log', deadline - time.monotonic(), cancel, progress)
    workers = min(threads, 8)
    jobs, shards = [], []
    for worker in range(workers):
        selected, aligned = directory / f'prefilter-{worker}', directory / f'alignments-{worker}'
        shards.append(aligned)
        jobs.append([
            [executable, 'filterdb', prefilter, selected, '--filter-expression',
             f'$1 % {workers} == {worker}', '--threads', '1'],
            [executable, 'structurealign', query_db, database, selected, aligned,
             '--threads', '1', '-a', '1', '--alignment-mode', '3', '--alignment-type', '2',
             '--sort-by-structure-bits', '1', '--comp-bias-corr', '1', '--comp-bias-corr-scale', '0.5',
             '--sub-mat', 'aa:3di.out,nucl:3di.out', '--min-aln-len', '20',
             '--exact-tmscore', '1', '-e', evalue, '--compressed', '0'],
        ])
    _parallel_scoring(jobs, directory, deadline, cancel, progress)
    _merge_alignment_shards(shards, alignment_db, prefilter, deadline, cancel)


def _run_search_pass(query_dir, queries, library, timeout=120, cancel=None, executable='foldseek', threads=DEFAULT_THREADS,
                     progress=None, manifest=None, max_seqs=1000, evalue=10, candidate_parallel=True):
    library = Path(library).resolve()
    if manifest is None:
        manifest = json.loads((library / 'manifest.json').read_text())
    if manifest.get('schema') not in (1, 2) or (manifest['schema'] == 1 and not isinstance(manifest.get('records'), list)):
        raise ValueError('Unsupported domain library manifest.')
    database = library_file(library, manifest['database'])
    if timeout <= 0 or not math.isfinite(timeout) or threads < 1:
        raise ValueError('Search timeout and thread count must be positive.')
    if cancel is not None and cancel.is_set():
        raise InterruptedError('Foldseek search cancelled.')
    started = time.monotonic()
    # Keep query construction and candidate discovery on their established
    # thread budget. Independent alignment/scoring shards share the CPU limit;
    # results are restored to the original ordering before adaptive passes.
    command_threads = min(4, threads)
    with tempfile.TemporaryDirectory(prefix='domain-foldseek-') as directory:
        directory = Path(directory)
        query_db, alignment_db = directory / 'query', directory / 'alignments'
        result = directory / 'hits.tsv'
        def run(*arguments):
            _run_foldseek([executable, *arguments], directory / 'log',
                          timeout - (time.monotonic() - started), cancel, progress)
        def convert(source, destination, fields):
            run('convertalis', query_db, database, source, destination, '--threads', command_threads,
                '--format-output', ','.join(fields), '--exact-tmscore', '1')
        # These are the same public stages and search options used by easy-search.
        # Keep exact structural scoring and the complete candidate budget intact.
        run('createdb', query_dir, query_db, '--threads', command_threads, '--chain-name-mode', '0')
        # Foldseek's alignment loop parallelizes queries. For a few long queries,
        # divide their *candidate lists*, retaining the full database/statistics.
        parallel = (candidate_parallel and command_threads > 1 and len(queries) <= command_threads // 2
                    and max(len(query['positions']) for query in queries.values()) >= 256
                    and library_reference_count(manifest) >= 512
                    and Path(str(query_db) + '_ss').is_file() and Path(str(database) + '_ss').is_file()
                    and not any(Path(str(database) + suffix).exists() for suffix in ('.idx', '_ss.idx')))
        if parallel:
            _candidate_parallel_alignment(query_db, database, alignment_db, directory / 'candidate-shards',
                                          command_threads, started + timeout, cancel, progress, executable, max_seqs, evalue)
        else:
            run('search', query_db, database, alignment_db, directory / 'work', '--threads', command_threads,
                '-a', '1', '--alignment-mode', '3', '--max-seqs', max_seqs, '-e', evalue,
                '-s', '9.5', '-k', '6', '--prefilter-mode', '1', '--min-aln-len', '20', '--exact-tmscore', '1')
        text = _screen_alignments(query_db, database, alignment_db, result, directory, threads,
                                  started + timeout, cancel, progress, executable, convert)
        metadata = hit_metadata(library, manifest, text, SCREEN_FIELDS)
        hits = parse_hits(text, queries, metadata, SCREEN_FIELDS)
        accepted = {hit['target'] for hit in hits if qualifying_hit(hit)}
        if accepted:
            # qTM and tTM are display-only. Compute them exactly for a superset
            # of accepted pairs; do not compute them for thousands of weak seeds.
            source, keys, found = alignment_db, [], set()
            lookup = Path(str(database) + '.lookup')
            if lookup.is_file():
                with lookup.open() as stream:
                    for line in stream:
                        values = line.rstrip('\n').split('\t')
                        if len(values) >= 2 and values[1].removesuffix('.pdb') in accepted:
                            keys.append(values[0])
                            found.add(values[1].removesuffix('.pdb'))
            if found == accepted:
                key_file = directory / 'accepted-targets'
                key_file.write_text(''.join(key + '\n' for key in keys))
                source = directory / 'accepted-alignments'
                run('filterdb', alignment_db, source, '--filter-file', key_file,
                    '--filter-column', '1', '--threads', command_threads)
            # A custom library without a usable lookup retains the old full
            # conversion behavior rather than losing annotations or score fields.
            detailed = directory / 'detailed.tsv'
            convert(source, detailed, FIELDS)
            full_hits = parse_hits(detailed.read_text(), queries, metadata)
            def key(hit):
                return tuple(hit[field] for field in ('query', 'target', 'qstart', 'qend', 'tstart', 'tend'))
            complete = {key(hit): hit for hit in full_hits}
            for index, hit in enumerate(hits):
                full = complete.get(key(hit))
                if full is not None:
                    if any(hit[field] != full[field] for field in SCREEN_FIELDS if field not in ('query', 'target')):
                        raise ValueError('Foldseek changed an alignment while calculating detailed scores.')
                    hits[index] = full
                elif qualifying_hit(hit):
                    raise ValueError('Foldseek omitted detailed scores for an accepted candidate.')
    return {'elapsed_seconds': time.monotonic() - started, 'pilot': manifest.get('pilot', False),
            'references': library_reference_count(manifest), 'queries': queries, 'hits': hits,
            'catalog': manifest.get('catalog', 'CDD')}


def review_windows(hits, queries, index_maps, tested):
    """Revisit provisional assignments without letting them block their own region."""
    stronger = {root: set() for root in queries}
    proposals = []
    for hit in select_matches(hits, per_family=False):
        root = hit['query']
        indices = {index_maps[root][p] for p in hit['aligned_positions']}
        first = max(0, min(indices) - 12)
        last = min(len(queries[root]['positions']), max(indices) + 13)
        # Keep context but stop at stronger neighboring placements. An overlap
        # near a boundary must not force the weaker placement to keep that cut.
        runs, start = [], None
        for index in range(first, last + 1):
            available = index < last and index not in stronger[root]
            if available and start is None:
                start = index
            elif not available and start is not None:
                runs.append((start, index))
                start = None
        if runs:
            a, b = max(runs, key=lambda span: sum(span[0] <= p < span[1] for p in indices))
            window = root, a, b
            if sum(a <= p < b for p in indices) >= 20 and (*window, 5000) not in tested:
                proposals.append((hit['evalue'], hit['bits'], window))
        stronger[root].update(indices)
    # Spend the review budget on the least certain existing assignments first.
    counts, selected = {}, []
    for _, _, window in sorted(proposals, key=lambda row: (-row[0], row[1])):
        root = window[0]
        if counts.get(root, 0) < 4 and window not in selected:
            selected.append(window)
            counts[root] = counts.get(root, 0) + 1
    return selected


def search_windows(queries, index_maps, covered, pass_hits, tested, depth):
    """Propose unassigned stretches and bounded searches around plausible hits."""
    windows = {}
    for root, query in queries.items():
        start = None
        for index in range(len(query['positions']) + 1):
            unassigned = index < len(query['positions']) and index not in covered[root]
            if unassigned and start is None:
                start = index
            elif not unassigned and start is not None:
                window = root, start, index
                if index - start >= 20 and (*window, depth) not in tested:
                    windows[f'residual{len(windows)}'] = window
                start = None
    focused = {root: [] for root in queries}
    for hit in sorted(pass_hits, key=lambda h: (h['evalue'], -h['bits'])):
        if (qualifying_hit(hit) or hit['evalue'] > 10 or hit['tcov'] < .5
                or hit['alntmscore'] < .5 or len(hit['aligned_positions']) < 20):
            continue
        root = hit['query']
        indices = {index_maps[root][p] for p in hit['aligned_positions']}
        if len(indices & covered[root]) >= .5 * len(indices) or len(focused[root]) >= 4:
            continue
        first, last = max(0, min(indices) - 12), min(len(queries[root]['positions']), max(indices) + 13)
        window = root, first, last
        if (*window, depth) in tested or window in windows.values():
            continue
        if any(max(0, min(last, b) - max(first, a)) / max(last - first, b - a) >= .7
               for a, b in focused[root]):
            continue
        focused[root].append((first, last))
        windows[f'focus{len(windows)}'] = window
    return windows


def guided_windows(queries, guidance):
    """Pad confident sequence boundaries and cover every remaining modeled residue."""
    from bisect import bisect_left, bisect_right

    windows = {}
    useful = False
    for root, query in queries.items():
        positions = query['positions']
        spans = set()
        for a, b in (guidance or {}).get(root, []):
            first, last = bisect_left(positions, a), bisect_right(positions, b)
            if last - first >= 20:
                spans.add((max(0, first - 12), min(len(positions), last + 12)))
        if not spans:
            spans.add((0, len(positions)))
        useful |= spans != {(0, len(positions))}
        covered = set(i for a, b in spans for i in range(a, b))
        start = None
        for i in range(len(positions) + 1):
            if i < len(positions) and i not in covered:
                if start is None:
                    start = i
            elif start is not None:
                a, b = max(0, start - 12), min(len(positions), i + 12)
                a = max(0, min(a, b - 20))
                b = min(len(positions), max(b, a + 20))
                spans.add((a, b))
                start = None
        for first, last in sorted(spans):
            windows[f'guided{len(windows)}'] = root, first, last
    return windows if useful else {}


def run_search(query_dir, queries, library, timeout=120, cancel=None, executable='foldseek', threads=DEFAULT_THREADS, progress=None,
               guidance=None, candidate_parallel=True):
    """Revisit unassigned stretches and focus plausible below-cutoff placements."""
    started = time.monotonic()
    # Keep all repeat passes on one immutable generation while a build may publish.
    manifest = json.loads((Path(library) / 'manifest.json').read_text())
    source_lines = {}
    for identifier, query in queries.items():
        lines = (Path(query_dir) / (identifier + '.pdb')).read_text().splitlines(keepends=True)
        source_lines[identifier] = [line for line in lines if line.startswith('ATOM  ') and line[12:16].strip() == 'CA']
        if len(source_lines[identifier]) != len(query['positions']):
            raise ValueError('Query coordinate file does not match its residue index map.')
    index_maps = {identifier: {p: i for i, p in enumerate(query['positions'])} for identifier, query in queries.items()}
    covered = {identifier: set() for identifier in queries}
    windows = {identifier: (identifier, 0, len(query['positions'])) for identifier, query in queries.items()}
    tested, hits, result, final_review = set(), [], None, False
    with tempfile.TemporaryDirectory(prefix='domain-repeat-') as directory:
        def prepare(next_windows, name):
            folder = Path(directory) / str(name)
            folder.mkdir()
            active = {}
            for identifier, (root, first, last) in next_windows.items():
                lines = [f'{line[:6]}{number:5d}{line[11:22]}{number:4d}{line[26:]}'
                         for number, line in enumerate(source_lines[root][first:last], 1)]
                (folder / (identifier + '.pdb')).write_text(''.join(lines) + 'TER\nEND\n')
                active[identifier] = {'positions': queries[root]['positions'][first:last],
                                      'label': queries[root].get('label', root)}
            return folder, active

        folder = Path(query_dir)
        active_queries = queries
        initial = guided_windows(queries, guidance)
        if initial:
            windows = initial
            folder, active_queries = prepare(windows, 'guided')
        review_after = 4 if initial else 3
        search_round = 0
        while True:
            search_round += 1
            guided = bool(initial) and search_round == 1
            max_seqs = 5000 if final_review else 1000
            remaining = timeout - (time.monotonic() - started)
            if remaining <= 0:
                raise TimeoutError('Foldseek search timed out.')
            phase = 'Pfam-guided search' if guided else 'final review' if final_review else 'search'
            message = f'{phase} pass {search_round}, {len(active_queries)} query region(s)'
            report = (lambda message=message: progress(message)) if progress is not None else None
            if report is not None:
                report()
            result = _run_search_pass(folder, active_queries, library, remaining, cancel, executable, threads,
                                      report, manifest, max_seqs=max_seqs, evalue=.01 if final_review or guided else 10,
                                      **({} if candidate_parallel else {'candidate_parallel': False}))
            for identifier, window in windows.items():
                if not guided:  # Strict guided screening must not block loose fallback searches.
                    tested.add((*window, max_seqs))
            pass_hits = []
            for hit in result['hits']:
                root, first, last = windows[hit['query']]
                hit['search_query'] = hit['query']
                hit['query'] = root
                hit['search_round'] = search_round
                hit['prefilter_limit'] = max_seqs
                hit['search_region'] = [queries[root]['positions'][first], queries[root]['positions'][last - 1]]
                hit['original_query_length'] = len(queries[root]['positions'])
                if guided:
                    hit['search_strategy'] = 'pfam_guided'
                hits.append(hit)
                pass_hits.append(hit)
            # Recompute from the best placements across all passes. Superseded
            # weaker assignments must not keep reserving extra query residues.
            covered = {identifier: set() for identifier in queries}
            for hit in select_matches(hits, per_family=False):
                root = hit['query']
                covered[root].update(index_maps[root][position] for position in hit['aligned_positions'])
            if final_review:
                # A target contributes one alignment per query region. Repeated
                # domains can therefore need more than four passes. Continue on
                # genuinely untested residuals at the review depth until exhausted;
                # the shared deadline/cancellation still bound all work.
                next_windows = {key: window for key, window in
                                search_windows(queries, index_maps, covered, pass_hits, tested, 5000).items()
                                if key.startswith('residual')}
                if not next_windows:
                    break
                folder, active_queries = prepare(next_windows, search_round)
                windows = next_windows
                continue
            next_windows = search_windows(queries, index_maps, covered, pass_hits, tested, 1000)
            if not next_windows or search_round >= review_after:
                final_review = True
                next_windows = search_windows(queries, index_maps, covered, pass_hits, tested, 5000)
                for window in review_windows(hits, queries, index_maps, tested):
                    if window not in next_windows.values():
                        next_windows[f'review{len(next_windows)}'] = window
            if not next_windows:
                break
            folder, active_queries = prepare(next_windows, search_round)
            windows = next_windows
    return {**result, 'queries': queries, 'hits': sorted(hits, key=lambda h: (h['evalue'], -h['bits'])),
            'elapsed_seconds': time.monotonic() - started, 'search_rounds': search_round}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('build')
    build.add_argument('--out', required=True)
    build.add_argument('--models', nargs='+', default=PILOT_MODELS)
    build.add_argument('--representatives', type=int, default=3)
    catalog = sub.add_parser('catalog', help='Index the official CDD seed archive across all model sources')
    catalog.add_argument('--out', required=True)
    catalog.add_argument('--archive', required=True, help='Downloaded NCBI fasta.tar.gz')
    catalog.add_argument('--cache', default='.domain_data/cdd-pilot/downloads')
    catalog.add_argument('--candidates', type=int, default=8, help='Maximum seed candidates per model per structure source')
    expand = sub.add_parser('expand', help='Resume building CDD-wide, sequence-verified domain references')
    expand.add_argument('--out', required=True)
    expand.add_argument('--catalog', required=True)
    expand.add_argument('--representatives', type=int, default=3)
    expand.add_argument('--limit', type=int, help='Maximum families with candidates to process this run; omit for all')
    expand.add_argument('--seed-library')
    expand.add_argument('--models', nargs='+', help='Process only these CDD accessions; other checkpoints remain intact')
    expand.add_argument('--no-alphafold', action='store_true')
    expand.add_argument('--min-plddt', type=float, default=70)
    expand.add_argument('--retry-failed', action='store_true')
    status = sub.add_parser('status', help='Show coverage from the last published library generation')
    status.add_argument('--library', required=True)
    benchmark = sub.add_parser('benchmark', help='Run a diverse public-protein structural control panel')
    benchmark.add_argument('--library', required=True)
    benchmark.add_argument('--out', required=True)
    ecod = sub.add_parser('install-ecod', help='Install the complete published ECOD F40 Foldseek database')
    ecod.add_argument('--out', required=True)
    ecod.add_argument('--release', default='v295.2')
    ecod.add_argument('--cdd-catalog', help='Optional catalog.sqlite for an explicit CDD model coverage audit')
    merge = sub.add_parser('merge', help='Combine two complete libraries, retaining per-reference provenance')
    merge.add_argument('--libraries', nargs=2, required=True)
    merge.add_argument('--out', required=True)
    merge.add_argument('--replace', action='store_true', help='Atomically extend an existing input library; retain its previous manifest')
    search = sub.add_parser('search')
    search.add_argument('structure')
    search.add_argument('--library', required=True)
    search.add_argument('--chain')
    search.add_argument('--segments')
    search.add_argument('--out', required=True)
    args = parser.parse_args()
    if args.command == 'build':
        manifest = build_library(args.out, args.models, args.representatives)
        print('Built', len(manifest['records']), 'references covering', len({r['cdd_accession'] for r in manifest['records']}), 'models')
    elif args.command == 'catalog':
        from domain_foldseek_catalog import catalog
        print(json.dumps(catalog(args.out, args.archive, args.cache, args.candidates), indent=2))
    elif args.command == 'expand':
        from domain_foldseek_catalog import expand
        manifest = expand(args.catalog, args.out, args.representatives, args.limit, args.seed_library,
                          not args.no_alphafold, args.min_plddt, args.retry_failed, models=args.models)
        print(json.dumps(manifest['coverage'], indent=2))
    elif args.command == 'status':
        print((Path(args.library) / 'coverage.json').read_text())
    elif args.command == 'benchmark':
        from domain_foldseek_catalog import benchmark
        benchmark(args.library, args.out)
    elif args.command == 'install-ecod':
        from domain_foldseek_ecod import install
        manifest = install(args.out, args.release, args.cdd_catalog)
        coverage = manifest['coverage']
        print(json.dumps({k: v for k, v in coverage.items() if k not in ('unmapped_target_ids', 'pfam_accessions')}, indent=2))
    elif args.command == 'merge':
        from domain_foldseek_ecod import merge
        manifest = merge(args.libraries, args.out, replace=args.replace)
        print('Merged', manifest['reference_count'], 'references')
    else:
        started = time.monotonic()
        with tempfile.TemporaryDirectory(prefix='domain-query-') as directory:
            queries = query_domains(args.structure, directory, args.segments, args.chain)
            result = run_search(directory, queries, args.library)
        result['total_seconds'] = time.monotonic() - started
        Path(args.out).write_text(json.dumps(result, indent=2) + '\n')
        print('Search:', result['elapsed_seconds'], 'seconds; total:', result['total_seconds'])
        for query in queries:
            print(query, [(h['reference']['name'], h['evalue'], h['qcov'], h['qtmscore']) for h in result['hits'] if h['query'] == query][:5])


if __name__ == '__main__':
    main()
