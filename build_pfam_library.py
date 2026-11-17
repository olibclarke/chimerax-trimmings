#!/usr/bin/env python3
"""Build a release-pinned Pfam structural library, with a strict completeness gate.

One invocation downloads the full profile inventory, annotates existing Foldseek
references in resumable HMMER shards, and fills gaps from PDB, AlphaFold, and an
optional local predictor. A missing structure never counts as structural coverage.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextlib import closing
from email.utils import parsedate_to_datetime
import csv
import gzip
import hashlib
from http.client import HTTPException
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import sqlite3
import subprocess
import tarfile
import tempfile
from threading import Event, Lock, BoundedSemaphore
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import domain_foldseek as fs
from domain_foldseek_catalog import atomic_json, build_lock
from domain_foldseek_ecod import checksum, checksums, staged_generation

VERSION = 2
PFAM = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{release}/'
SIFTS = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_pfam.csv.gz'
ACCESSION = re.compile(r'PF\d{5,}')
FOLDED_PROFILE_TYPES = frozenset(('Domain', 'Family', 'Repeat'))
UNIPROT = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?')


class IncompleteLibrary(RuntimeError):
    pass


class RangeUnsupported(RuntimeError):
    pass


def say(message):
    print(time.strftime('%H:%M:%S'), message, flush=True)


def run(command, log, cancel=None, timeout=None):
    """Terminate our process group on cancellation; preserve the diagnostic log."""
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with Path(log).open('w') as stream:
        process = subprocess.Popen(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            while process.poll() is None:
                if cancel is not None and cancel.is_set():
                    raise InterruptedError('Build cancelled.')
                if timeout and time.monotonic() - started > timeout:
                    raise TimeoutError('Command timed out; see ' + str(log))
                time.sleep(.1)
            if process.returncode:
                raise RuntimeError(f'{command[0]} exited {process.returncode}; see {log}')
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()


def download(url, path, expected=None):
    """Checksum-verified, resumable bulk download; corrupt partials start afresh."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and (expected is None or checksum(path) == expected):
        return path
    partial = path.with_name(path.name + '.part')
    if path.exists():
        path.unlink()
    if partial.exists() and expected and checksum(partial) == expected:
        partial.replace(path)
        return path
    if expected:
        try:
            if parallel_download(url, path, expected):
                return path
        except RangeUnsupported:
            say('Server did not honor byte ranges; using the resumable single transfer')
    say('Downloading ' + path.name)
    result = subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--connect-timeout', '20',
                             '--continue-at', '-', '--output', str(partial), url])
    if result.returncode:
        # Range requests can be refused by a proxy, or a completed stale partial.
        if result.returncode in (22, 33) and partial.exists():
            partial.unlink()
            subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--connect-timeout', '20',
                            '--output', str(partial), url], check=True)
        else:
            result.check_returncode()
    if expected and checksum(partial) != expected:
        partial.unlink()
        raise ValueError('Checksum mismatch: ' + str(path))
    partial.replace(path)
    return path


def parallel_download(url, path, expected, chunk_size=4 * 1024 * 1024, workers=8, minimum_size=64 * 1024 * 1024):
    """Resume checked byte ranges; the final upstream checksum is authoritative."""
    try:
        request = Request(url, method='HEAD')
        with urlopen(request, timeout=15) as response:
            size = int(response.headers.get('Content-Length', '0'))
            ranges = response.headers.get('Accept-Ranges', '').lower()
            etag = response.headers.get('ETag')
    except (OSError, ValueError, HTTPException):
        return False
    if size < minimum_size or ranges != 'bytes':
        return False
    path = Path(path)
    folder = path.with_name(path.name + '.ranges')
    folder.mkdir(parents=True, exist_ok=True)
    metadata = {'url': url, 'size': size, 'md5': expected, 'chunk_size': chunk_size}
    marker = folder / 'source.json'
    if marker.exists() and json.loads(marker.read_text()) != metadata:
        shutil.rmtree(folder)
        folder.mkdir()
    atomic_json(marker, metadata)
    partial = path.with_name(path.name + '.part')
    pieces = [(i, first, min(size, first + chunk_size) - 1)
              for i, first in enumerate(range(0, size, chunk_size))]
    def destination(piece):
        return folder / f'{piece[0]:08d}.bin'
    if partial.exists():
        # Reuse the complete ranges in an earlier single-stream partial.
        with partial.open('rb') as source:
            for piece in pieces:
                length = piece[2] - piece[1] + 1
                if piece[1] + length > partial.stat().st_size:
                    break
                target = destination(piece)
                if not target.exists() or target.stat().st_size != length:
                    source.seek(piece[1])
                    target.write_bytes(source.read(length))
    pending = [piece for piece in pieces if not destination(piece).exists()
               or destination(piece).stat().st_size != piece[2] - piece[1] + 1]
    cancel = Event()
    started, complete = time.monotonic(), len(pieces) - len(pending)
    say(f'Parallel download {path.name}: {complete}/{len(pieces)} ranges cached; {workers} connections')
    def transfer(piece):
        _, first, last = piece
        target = destination(piece)
        for attempt in range(4):
            if cancel.is_set():
                raise InterruptedError('Download cancelled.')
            temporary = target.with_suffix('.part')
            try:
                headers = {'Range': f'bytes={first}-{last}'}
                if etag:
                    headers['If-Range'] = etag
                with urlopen(Request(url, headers=headers), timeout=30) as response:
                    if (response.status != 206 or response.headers.get('Content-Range') != f'bytes {first}-{last}/{size}'):
                        raise RangeUnsupported('Server ignored the requested byte range.')
                    with temporary.open('wb') as stream:
                        while True:
                            if cancel.is_set():
                                raise InterruptedError('Download cancelled.')
                            block = response.read(256 * 1024)
                            if not block:
                                break
                            stream.write(block)
                if temporary.stat().st_size != last - first + 1:
                    raise OSError('Incomplete download range.')
                temporary.replace(target)
                return
            except (OSError, HTTPException):
                if attempt == 3:
                    raise
                if cancel.wait(2 ** attempt):
                    raise InterruptedError('Download cancelled.')
            finally:
                temporary.unlink(missing_ok=True)
    for _, _ in parallel_jobs(pending, transfer, workers, cancel):
        complete += 1
        say(f'{path.name}: {complete}/{len(pieces)} ranges; {time.monotonic() - started:.0f}s elapsed')
    with partial.open('wb') as output:
        for piece in pieces:
            with destination(piece).open('rb') as source:
                shutil.copyfileobj(source, output)
    if partial.stat().st_size != size or checksum(partial) != expected:
        partial.unlink()
        shutil.rmtree(folder)
        raise ValueError('Parallel download checksum mismatch: ' + str(path))
    partial.replace(path)
    shutil.rmtree(folder)
    return True


def release_assets(output, release, supplied=None, audit=False):
    folder = Path(output) / 'downloads' / ('Pfam' + release)
    folder.mkdir(parents=True, exist_ok=True)
    base = PFAM.format(release=release)
    names = ['md5_checksums', 'Pfam.version.gz', 'Pfam-A.hmm.dat.gz']
    if not audit:
        names += ['Pfam-A.hmm.gz', 'Pfam-A.seed.gz']
    expected = {}
    for name in names:
        path = folder / name
        source = Path(supplied) / name if supplied else None
        if source is not None and source.is_file() and not path.exists():
            shutil.copyfile(source, path)
        download(base + name, path, expected.get(name))
        if name == 'md5_checksums':
            expected = checksums(path.read_text())
            if set(names[1:]) - expected.keys():
                raise ValueError('Release checksum manifest omits required Pfam files.')
    with gzip.open(folder / 'Pfam.version.gz', 'rt') as stream:
        version = stream.read()
    if not re.search(r'Pfam release\s*:\s*' + re.escape(release) + r'\s*$', version, re.M):
        raise ValueError('Pfam.version does not match the requested release.')
    count = re.search(r'Pfam-A families\s*:\s*(\d+)', version)
    if count is None:
        raise ValueError('Pfam.version omits the family count.')
    return folder, expected, int(count[1])


def catalog(path, expected_count):
    families, record = {}, {}
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#=GF '):
                key, value = line[5:].split(None, 1)
                record[key] = value.strip()
            elif line.strip() == '//':
                version = record['AC']
                accession = version.split('.')[0]
                if not ACCESSION.fullmatch(accession) or accession in families:
                    raise ValueError('Invalid/duplicate Pfam accession: ' + accession)
                ga = [float(value.strip(';')) for value in record['GA'].split()]
                if len(ga) != 2 or not all(math.isfinite(value) for value in ga):
                    raise ValueError('Invalid gathering thresholds: ' + accession)
                families[accession] = {'accession': accession, 'version': version, 'name': record['ID'],
                    'description': record.get('DE', ''), 'type': record['TP'], 'length': int(record['ML']), 'ga': ga}
                record = {}
    if record or len(families) != expected_count:
        raise ValueError(f'Pfam catalog incomplete: found {len(families)}, expected {expected_count}.')
    return families


class Library:
    def __init__(self, folder):
        self.folder = Path(folder).resolve()
        text = (self.folder / 'manifest.json').read_bytes()
        self.manifest = json.loads(text)
        self.digest = hashlib.sha256(text).hexdigest()
        self.database = fs.library_file(self.folder, self.manifest['database'])
        self.metadata = (fs.library_file(self.folder, self.manifest['reference_metadata'])
                         if self.manifest['schema'] == 2 else None)
        self.small_records = {record['id']: record for record in self.manifest.get('records', [])}

    def records(self):
        if self.metadata is None:
            yield from self.small_records.values()
        else:
            with closing(sqlite3.connect(f'file:{self.metadata}?mode=ro', uri=True)) as connection:
                for (text,) in connection.execute('SELECT record FROM refs'):
                    yield json.loads(text)

    def record(self, identifier):
        if self.metadata is None:
            return self.small_records[identifier]
        with closing(sqlite3.connect(f'file:{self.metadata}?mode=ro', uri=True)) as connection:
            row = connection.execute('SELECT record FROM refs WHERE id=?', (identifier,)).fetchone()
        if row is None:
            raise ValueError('Reference lacks metadata: ' + identifier)
        return json.loads(row[0])


def mapping_audit(output, families, libraries, release):
    mapped = set()
    for library in libraries:
        for record in library.records():
            mapped.update(record.get('pfam_accessions', []))
            if record.get('database') == 'pfam':
                mapped.add(record['accession'])
    covered = families.keys() & mapped
    summary = {'pfam_release': release, 'families': len(families), 'explicit_family_mappings': len(covered),
               'missing_mappings': len(families) - len(covered),
               'meaning': 'Cross-references, not newly validated profile matches or proof of absent folds.'}
    atomic_json(Path(output) / 'mapping-audit.json', summary)
    with (Path(output) / 'mapping-audit.tsv').open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['accession', 'name', 'type', 'status'])
        for accession, family in sorted(families.items()):
            writer.writerow([accession, family['name'], family['type'], 'mapped' if accession in covered else 'unmapped'])
    return summary


def state_database(path):
    connection = sqlite3.connect(path)
    connection.executescript('''
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT);
        CREATE TABLE IF NOT EXISTS targets(id TEXT PRIMARY KEY,library INTEGER,dbkey INTEGER,original TEXT,sequence TEXT);
        CREATE INDEX IF NOT EXISTS target_source ON targets(library,original);
        CREATE TABLE IF NOT EXISTS target_groups(target TEXT PRIMARY KEY,group_id TEXT);
        CREATE TABLE IF NOT EXISTS hits(family TEXT,target TEXT,first INTEGER,last INTEGER,score REAL,evalue REAL,
            coverage REAL,tried INTEGER DEFAULT 0,PRIMARY KEY(family,target,first,last));
        CREATE INDEX IF NOT EXISTS hits_family ON hits(family,tried,evalue,score);
        CREATE TABLE IF NOT EXISTS seeds(family TEXT,uniprot TEXT,header TEXT,sequence TEXT,
            UNIQUE(family,uniprot,sequence));
        CREATE INDEX IF NOT EXISTS seeds_family ON seeds(family);
        CREATE TABLE IF NOT EXISTS pdbs(family TEXT,pdb TEXT,chain TEXT,UNIQUE(family,pdb,chain));
        CREATE INDEX IF NOT EXISTS pdbs_family ON pdbs(family);
        CREATE TABLE IF NOT EXISTS refs(id TEXT PRIMARY KEY,family TEXT,sequence TEXT,record TEXT,path TEXT,sha256 TEXT,
            UNIQUE(family,sequence));
        CREATE TABLE IF NOT EXISTS runs(family TEXT PRIMARY KEY,status TEXT,errors TEXT);
    ''')
    return connection


def done(state, key):
    return state.execute('SELECT 1 FROM meta WHERE key=?', (key,)).fetchone() is not None


def mark(state, key):
    state.execute('INSERT OR REPLACE INTO meta VALUES(?,?)', (key, 'done'))
    state.commit()


def profiles(assets, folder, families, selected, size=250):
    """Keep every HMM in the bundle; shard only the requested scan scope."""
    folder.mkdir(parents=True, exist_ok=True)
    library = folder / 'Pfam-A.hmm'
    if not library.exists():
        partial = library.with_suffix('.part')
        with gzip.open(assets / 'Pfam-A.hmm.gz', 'rb') as source, partial.open('wb') as target:
            shutil.copyfileobj(source, target)
        partial.replace(library)
    marker = folder / 'shards.json'
    signature = hashlib.sha256(json.dumps([sorted(selected), size]).encode()).hexdigest()
    if marker.exists():
        saved = json.loads(marker.read_text())
        if checksum(library, 'sha256') != saved.get('library_sha256'):
            raise ValueError('The decompressed Pfam profile library changed; remove it and rerun to restore it from the verified archive.')
        if saved.get('signature') == signature:
            paths = [folder / path for path in saved['paths']]
            if all(path.is_file() for path in paths):
                return library, paths
    shard_folder = folder / ('shards-' + signature[:16])
    shard_folder.mkdir(exist_ok=True)
    block, seen, paths, count = [], set(), [], 0
    output = None
    try:
        with library.open() as stream:
            for line in stream:
                block.append(line)
                if line.strip() != '//':
                    continue
                accessions = [line.split()[1] for line in block if line.startswith('ACC ')]
                if len(accessions) != 1:
                    raise ValueError('HMM profile missing a unique accession.')
                version = accessions[0]
                accession = version.split('.')[0]
                if accession not in families or version != families[accession]['version'] or accession in seen:
                    raise ValueError('Profile library/catalog version mismatch: ' + version)
                seen.add(accession)
                if accession in selected:
                    if count % size == 0:
                        if output:
                            output.close()
                        path = shard_folder / f'{len(paths):05d}.hmm'
                        paths.append(path)
                        output = path.open('w')
                    output.writelines(block)
                    count += 1
                block = []
        if block or seen != families.keys():
            raise ValueError('The HMM library does not contain exactly the complete Pfam catalog.')
    finally:
        if output:
            output.close()
    atomic_json(marker, {'signature': signature, 'library_sha256': checksum(library, 'sha256'),
                         'paths': [str(p.relative_to(folder)) for p in paths]})
    return library, paths


def fasta_records(path):
    name, sequence = None, []
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(sequence)
                name, sequence = line[1:].split()[0], []
            else:
                sequence.append(line.strip())
        if name is not None:
            yield name, ''.join(sequence)


def export_references(state, libraries, work, executable, cancel):
    fasta = work / 'references.fasta'
    if done(state, 'export') and fasta.exists():
        return fasta
    state.execute('DELETE FROM targets')
    with fasta.open('w') as output:
        for index, library in enumerate(libraries):
            source = work / f'source-{index}.fasta'
            run([executable, 'convert2fasta', library.database, source], work / f'export-{index}.log', cancel)
            with Path(str(library.database) + '.lookup').open() as lookup:
                keys = {row[1]: int(row[0]) for row in csv.reader(lookup, delimiter='\t')}
            for identifier, sequence in fasta_records(source):
                identifier = identifier.removesuffix('.pdb')
                if identifier not in keys or not re.fullmatch(r'[A-Za-z0-9_.:-]+', identifier):
                    raise ValueError('Invalid/unindexed source reference ID: ' + identifier)
                sid = f'r{index}_{keys[identifier]}'
                state.execute('INSERT INTO targets VALUES(?,?,?,?,?)', (sid, index, keys[identifier], identifier, sequence))
                output.write(f'>{sid}\n{sequence}\n')
    mark(state, 'export')
    return fasta


def parse_domtbl(path, families, min_coverage):
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('#') or not line.strip():
                continue
            row = line.split(maxsplit=22)
            if len(row) < 22:
                raise ValueError('Truncated HMMER domain table.')
            accession = row[4].split('.')[0]
            if accession not in families or row[4] != families[accession]['version']:
                raise ValueError('Unexpected profile accession/version in HMMER output: ' + row[4])
            length, hlen = int(row[2]), int(row[5])
            hfirst, hlast, first, last = map(int, (row[15], row[16], row[19], row[20]))
            evalue, score = float(row[12]), float(row[13])
            coverage = (hlast - hfirst + 1) / hlen
            if not (1 <= hfirst <= hlast <= hlen and 1 <= first <= last <= length and hlen == families[accession]['length']
                    and math.isfinite(evalue) and evalue >= 0 and math.isfinite(score)):
                raise ValueError('Invalid HMMER scores or domain coordinates.')
            if coverage >= min_coverage:
                yield accession, row[0], first, last, score, evalue, coverage


def parallel_jobs(items, function, workers, cancel):
    """Bound queued work as well as running threads; yield results to one DB writer."""
    iterator = iter(items)
    executor = ThreadPoolExecutor(max_workers=workers)
    pending = {}
    started = last_report = time.monotonic()
    try:
        for _ in range(workers * 2):
            item = next(iterator, None)
            if item is not None:
                pending[executor.submit(function, item)] = item
        while pending:
            ready, _ = wait(pending, return_when=FIRST_COMPLETED, timeout=.5)
            if time.monotonic() - last_report >= 10:
                say(f'Working: {len(pending)} queued/running jobs; {time.monotonic() - started:.0f}s elapsed')
                last_report = time.monotonic()
            if cancel.is_set():
                raise InterruptedError('Build cancelled.')
            for future in ready:
                item = pending.pop(future)
                yield item, future.result()
                following = next(iterator, None)
                if following is not None:
                    pending[executor.submit(function, following)] = following
    except BaseException:
        cancel.set()
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)


def scan_references(state, shards, fasta, work, families, args, cancel, prefix=''):
    if not state.execute('SELECT 1 FROM targets LIMIT 1').fetchone():
        return
    jobs = [(index, shard) for index, shard in enumerate(shards) if not done(state, f'{prefix}scan-{index}')]
    def search(item):
        index, shard = item
        table = work / f'scan-{index}.tsv'
        if not table.exists():
            partial = table.with_suffix('.part')
            run([args.hmmsearch, '--cut_ga', '--noali', '--cpu', max(1, args.threads // args.scan_jobs),
                 '--domtblout', partial, shard, fasta], work / f'scan-{index}.log', cancel)
            partial.replace(table)
        return table
    for (index, _), table in parallel_jobs(jobs, search, args.scan_jobs, cancel):
        with state:
            for row in parse_domtbl(table, families, args.min_coverage):
                if not state.execute('SELECT 1 FROM targets WHERE id=?', (row[1],)).fetchone():
                    raise ValueError('HMMER returned an unknown reference.')
                state.execute('INSERT OR IGNORE INTO hits(family,target,first,last,score,evalue,coverage) VALUES(?,?,?,?,?,?,?)', row)
            state.execute('INSERT INTO meta VALUES(?,?)', (f'{prefix}scan-{index}', 'done'))
        completed = sum(done(state, f'{prefix}scan-{i}') for i in range(len(shards)))
        say(f'Profile scan {prefix}: {completed}/{len(shards)} shards complete')
        atomic_json(Path(args.out) / 'progress.json', {'phase': 'reference scan', 'batch': prefix,
                    'completed': completed, 'total': len(shards), 'covered_families': len(represented(state))})


def reference_batches(state, libraries, work, initial_size):
    """Visit distinct source families first, then exhaust every remaining target."""
    marker = work / 'reference-batches.json'
    if marker.exists():
        data = json.loads(marker.read_text())
        paths = [work / name for name in data['paths']]
        if all(path.is_file() for path in paths):
            return paths
    with state:
        state.execute('DELETE FROM target_groups')
        for index, library in enumerate(libraries):
            for record in library.records():
                group = str(record.get('database', '')) + ':' + str(record.get('accession', record['id']))
                group = hashlib.sha256(group.encode()).hexdigest()
                state.execute('INSERT INTO target_groups SELECT id,? FROM targets WHERE library=? AND original=?',
                              (group, index, record['id']))
        if state.execute('SELECT count(*) FROM target_groups').fetchone()[0] != state.execute('SELECT count(*) FROM targets').fetchone()[0]:
            raise ValueError('Source metadata does not account for every exported target.')
    rows = state.execute('''SELECT id,sequence FROM (
        SELECT targets.id,targets.sequence,target_groups.group_id,
               row_number() OVER (PARTITION BY group_id ORDER BY length(sequence) DESC,targets.id) AS rank
        FROM targets JOIN target_groups ON targets.id=target_groups.target)
        ORDER BY rank,group_id,id''')
    folder = work / 'reference-batches'
    folder.mkdir(exist_ok=True)
    paths, stream, count, size = [], None, 0, initial_size
    try:
        for identifier, sequence in rows:
            if stream is None or count == size:
                if stream is not None:
                    stream.close()
                    size *= 2
                path = folder / f'{len(paths):04d}.fasta'
                paths.append(path)
                stream, count = path.open('w'), 0
            stream.write(f'>{identifier}\n{sequence}\n')
            count += 1
    finally:
        if stream is not None:
            stream.close()
    atomic_json(marker, {'paths': [str(path.relative_to(work)) for path in paths],
                         'strategy': 'Breadth first by source family; doubling disjoint batches; no targets omitted.'})
    return paths


def subset_shards(library, selected, families, folder, args, cancel):
    """Use the HMMER index to retrieve only still-unrepresented profiles."""
    marker = folder / 'shards.json'
    if marker.exists():
        paths = [folder / name for name in json.loads(marker.read_text())['paths']]
        if all(path.is_file() for path in paths):
            return paths
    keys, subset = folder / 'keys.txt', folder / 'selected.hmm'
    keys.write_text(''.join(families[accession]['version'] + '\n' for accession in sorted(selected)))
    run([args.hmmfetch, '-f', '-o', subset, library, keys], folder / 'fetch.log', cancel)
    paths, block, seen, stream = [], [], set(), None
    try:
        with subset.open() as source:
            for line in source:
                block.append(line)
                if line.strip() != '//':
                    continue
                version = next(line.split()[1] for line in block if line.startswith('ACC '))
                accession = version.split('.')[0]
                if accession not in selected or families[accession]['version'] != version or accession in seen:
                    raise ValueError('Unexpected/duplicate profile in adaptive scan.')
                if len(seen) % args.shard_size == 0:
                    if stream is not None:
                        stream.close()
                    path = folder / f'{len(paths):05d}.hmm'
                    paths.append(path)
                    stream = path.open('w')
                stream.writelines(block)
                seen.add(accession)
                block = []
        if block or seen != selected:
            raise ValueError('Adaptive profile extraction omitted requested models.')
    finally:
        if stream is not None:
            stream.close()
    subset.unlink()
    atomic_json(marker, {'paths': [str(path.relative_to(folder)) for path in paths]})
    return paths


def adaptive_scan(state, libraries, profile_library, original_shards, work, families, selected, args, cancel):
    paths = reference_batches(state, libraries, work, args.reference_batch_size)
    for index, fasta in enumerate(paths):
        key = f'reference-batch-{index}'
        if done(state, key):
            continue
        folder = work / 'adaptive-scans' / str(index)
        folder.mkdir(parents=True, exist_ok=True)
        snapshot = folder / 'active.json'
        if snapshot.exists():
            active = set(json.loads(snapshot.read_text()))
        else:
            counts = represented(state)
            active = {accession for accession in selected if counts.get(accession, 0) < args.representatives}
            atomic_json(snapshot, sorted(active))
        if not active:
            break
        say(f'Reference batch {index + 1}/{len(paths)}: searching {len(active)} still-needed Pfam profiles')
        shards = original_shards if active == selected else subset_shards(profile_library, active, families, folder, args, cancel)
        scan_references(state, shards, fasta, folder, families, args, cancel, prefix=f'batch-{index}-')
        reuse_structures(state, libraries, work, {accession: families[accession] for accession in active}, args, cancel)
        mark(state, key)
        coverage_report(Path(args.out), state, families, selected, args)
        say(f'Reference batch {index + 1} complete: {len(represented(state))}/{len(families)} families covered')


def pfam_reference_record(record):
    """Label the validated family as Pfam while retaining coordinate provenance."""
    record = dict(record)
    if record.get('database') != 'pfam':
        raise ValueError('Expected a sequence-validated Pfam reference.')
    if record.get('catalog') and record['catalog'] != 'Pfam':
        record.setdefault('parent_catalog', record['catalog'])
    record['catalog'] = 'Pfam'
    return record


def save_reference(family, residues, output, source, args, **provenance):
    import gemmi

    sequence = ''.join(fs.residue_letter(residue) for residue in residues)
    if not args.min_residues <= len(sequence) <= 9999:
        raise ValueError(f'Coordinate crop has {len(sequence)} residues; required {args.min_residues}–9999.')
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    geometry = hashlib.sha256(repr([(tuple(r.find_atom('CA', '*').pos), r.find_atom('CA', '*').b_iso)
                                   for r in residues]).encode()).hexdigest()
    identifier = family['accession'] + '_' + digest[:16] + '_' + geometry[:12]
    path = output / 'domains' / (identifier + '.pdb')
    path.parent.mkdir(parents=True, exist_ok=True)
    coords = [tuple(residue.find_atom('CA', '*').pos) for residue in residues]
    if not all(math.isfinite(value) for point in coords for value in point):
        raise ValueError('Non-finite reference coordinates.')
    lower = tuple(min(point[axis] for point in coords) for axis in range(3))
    # Reference placement is arbitrary. Translate each crop into the positive
    # PDB coordinate range so large source assemblies do not cause false gaps.
    lines = []
    for number, (residue, point, letter) in enumerate(zip(residues, coords, sequence), 1):
        x, y, z = (value - offset for value, offset in zip(point, lower))
        if max(x, y, z) > 9999:
            raise ValueError('Reference domain extent exceeds the coordinate export limit.')
        name = gemmi.expand_one_letter(letter, gemmi.ResidueKind.AA)
        atom = residue.find_atom('CA', '*')
        lines.append(f'ATOM  {number:5d}  CA  {name:>3s} A{number:4d}    '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{atom.b_iso:6.2f}           C\n')
    path.write_text(''.join(lines) + 'TER\nEND\n')
    accession = family['accession']
    record = {**provenance, 'id': identifier, 'accession': accession, 'catalog_accession': accession,
              'database': 'pfam', 'name': family['name'], 'description': family['description'],
              'pfam_accessions': [accession], 'pfam_version': family['version'], 'pfam_release': args.release,
              'sequence': sequence, 'length': len(sequence), 'structure_source': source,
              'coordinate_translation': [-value for value in lower],
              'entry_url': f'https://www.ebi.ac.uk/interpro/entry/pfam/{accession}/',
              'structures_url': f'https://www.ebi.ac.uk/interpro/structure/pdb/entry/pfam/{accession}/'}
    return pfam_reference_record(record), path


def add_reference(state, item, *, commit=True):
    record, path = item
    cursor = state.execute('INSERT OR IGNORE INTO refs VALUES(?,?,?,?,?,?)',
                  (record['id'], record['accession'], record['sequence'], json.dumps(record), str(path), checksum(path, 'sha256')))
    if commit:
        state.commit()
    return cursor.rowcount == 1


def represented(state):
    return dict(state.execute('SELECT family,count(*) FROM refs GROUP BY family'))


def reuse_structures(state, libraries, work, families, args, cancel):
    """Extract the next best coordinate target for each uncovered family in batches."""
    import gemmi
    while True:
        counts = represented(state)
        selected = []
        for accession in families:
            if counts.get(accession, 0) >= args.representatives:
                continue
            row = state.execute('SELECT family,target,first,last,score,evalue,coverage FROM hits '
                                'WHERE family=? AND tried=0 ORDER BY evalue,score DESC LIMIT 1', (accession,)).fetchone()
            if row:
                selected.append(row)
        if not selected:
            return
        by_library = defaultdict(dict)
        for row in selected:
            target = state.execute('SELECT library,dbkey,original,sequence FROM targets WHERE id=?', (row[1],)).fetchone()
            by_library[target[0]][target[1]] = target
        with tempfile.TemporaryDirectory(dir=work, prefix='reuse-') as temporary:
            folder = Path(temporary)
            for index, targets in by_library.items():
                keys = folder / f'{index}.keys'
                keys.write_text(''.join(f'{key}\n' for key in targets))
                db, pdbs = folder / f'db-{index}', folder / f'pdbs-{index}'
                run([args.foldseek, 'createsubdb', keys, libraries[index].database, db, '-v', '1'], folder / f'sub-{index}.log', cancel)
                run([args.foldseek, 'convert2pdb', db, pdbs, '--pdb-output-mode', '1', '--threads', args.threads, '-v', '1'],
                    folder / f'pdb-{index}.log', cancel)
            for processed, (accession, sid, first, last, score, evalue, coverage) in enumerate(selected, 1):
                if cancel.is_set():
                    raise InterruptedError('Build cancelled.')
                index, key, original, sequence = state.execute('SELECT library,dbkey,original,sequence FROM targets WHERE id=?', (sid,)).fetchone()
                try:
                    pdb = folder / f'pdbs-{index}' / (original + '.pdb')
                    structure = gemmi.read_structure(str(pdb))
                    residues = [r for chain in structure[0] for r in chain if r.find_atom('CA', '*')]
                    if ''.join(fs.residue_letter(r) for r in residues) != sequence:
                        raise ValueError('Exported coordinates do not match the source database sequence.')
                    metadata = libraries[index].record(original)
                    source = metadata.get('structure_source', 'reused')
                    if source not in ('pdb', 'alphafold'):
                        source = 'reused'
                    item = save_reference(families[accession], residues[first - 1:last], work, source, args,
                        catalog=metadata.get('catalog', 'Pfam'), pdb=metadata.get('pdb'), chain=metadata.get('chain', 'A'),
                        ecod_domain=metadata.get('ecod_domain'), uniprot=metadata.get('uniprot'),
                        structure_title=metadata.get('structure_title', original),
                        structure_url=metadata.get('structure_url') or metadata.get('entry_url'),
                        evidence='Pfam curated gathering-threshold sequence-profile match to an existing structural reference',
                        parent_library=str(libraries[index].folder), parent_manifest_sha256=libraries[index].digest,
                        parent_target=original, segments=[(first, last)],
                        profile_match={'bitscore': score, 'evalue': evalue, 'hmm_coverage': coverage})
                    add_reference(state, item, commit=False)
                except (OSError, ValueError, RuntimeError) as error:
                    say(f'{accession}: reuse candidate rejected: {error}')
                state.execute('UPDATE hits SET tried=1 WHERE family=? AND target=? AND first=? AND last=?',
                              (accession, sid, first, last))
                # Checkpoint bounded batches rather than doing two durable disk
                # flushes per crop. A crash can repeat at most 127 cheap crops.
                if processed % 128 == 0:
                    state.commit()
                    covered = len(represented(state))
                    say(f'Coordinate verification: {processed}/{len(selected)} candidates; {covered} families covered')
                    atomic_json(Path(args.out) / 'progress.json', {'phase': 'coordinate verification',
                                'completed': processed, 'total': len(selected), 'covered_families': covered})
            state.commit()
        say(f'Reuse stage: {len(represented(state))}/{len(families)} families have verified coordinates')


def index_seeds(state, path, families):
    if done(state, 'seeds'):
        return
    accession, sequences, accessions, seen = None, {}, {}, set()
    with state:
        state.execute('DELETE FROM seeds')
        with gzip.open(path, 'rt') as stream:
            for line in stream:
                if line.startswith('#=GF AC '):
                    version = line.split()[2]
                    accession = version.split('.')[0]
                    if accession not in families or version != families[accession]['version'] or accession in seen:
                        raise ValueError('Seed/catalog accession or version mismatch: ' + version)
                    seen.add(accession)
                elif line.startswith('#=GS '):
                    row = line.split()
                    if len(row) >= 4 and row[2] == 'AC':
                        accessions[row[1]] = row[3].rstrip(';').split('.')[0]
                elif line.strip() == '//':
                    if accession in families:
                        for header, pieces in sequences.items():
                            sequence = ''.join(pieces).replace('.', '').replace('-', '').upper()
                            uniprot = accessions.get(header, header.split('/')[0])
                            if UNIPROT.fullmatch(uniprot) and sequence and re.fullmatch('[A-Z]+', sequence):
                                state.execute('INSERT OR IGNORE INTO seeds VALUES(?,?,?,?)', (accession, uniprot, header, sequence))
                    accession, sequences, accessions = None, {}, {}
                elif line.strip() and not line.startswith('#'):
                    row = line.split()
                    if len(row) != 2:
                        raise ValueError('Invalid Stockholm seed alignment.')
                    sequences.setdefault(row[0], []).append(row[1])
        if accession is not None or sequences or seen != families.keys():
            raise ValueError('Truncated Pfam seed archive.')
    mark(state, 'seeds')


def index_pdbs(state, path, families):
    if done(state, 'sifts'):
        return
    with state, gzip.open(path, 'rt') as stream:
        state.execute('DELETE FROM pdbs')
        rows = csv.DictReader(line for line in stream if not line.startswith('#'))
        for row in rows:
            accession = row['PFAM_ID'].split('.')[0]
            pdb = row['PDB'].lower()
            if accession in families and re.fullmatch(r'[a-z0-9]{4}', pdb):
                state.execute('INSERT OR IGNORE INTO pdbs VALUES(?,?,?)', (accession, pdb, row['CHAIN']))
    mark(state, 'sifts')


class PublicCache:
    """Deduplicate concurrent downloads, retry transient failures, cache 404s briefly."""
    def __init__(self, folder, cancel):
        self.folder, self.cancel = Path(folder), cancel
        self.folder.mkdir(parents=True, exist_ok=True)
        # Fixed-size striped locks keep memory bounded even across millions of
        # candidate URLs. Identical URLs always share the same lock.
        self.locks = [Lock() for _ in range(1024)]
        self.host_lock, self.host_next, self.host_cooldown = Lock(), {}, {}
        self.request_starts = 0

    def wait_for_host(self, host):
        while True:
            if self.cancel.is_set():
                raise InterruptedError('Build cancelled.')
            with self.host_lock:
                now = time.monotonic()
                delay = max(self.host_next.get(host, 0), self.host_cooldown.get(host, 0)) - now
                if delay <= 0:
                    self.host_next[host] = now + .5
                    self.request_starts += 1
                    return
            if self.cancel.wait(delay):
                raise InterruptedError('Build cancelled.')

    def get(self, url):
        if not url.startswith('https://'):
            raise ValueError('Public resource URLs must use HTTPS.')
        key = hashlib.sha256(url.encode()).hexdigest()
        lock = self.locks[int(key[:8], 16) % len(self.locks)]
        with lock:
            path, negative = self.folder / key, self.folder / (key + '.404')
            if path.exists() and '/api/' in urlparse(url).path and time.time() - path.stat().st_mtime >= 86400:
                path.unlink()  # Dynamic candidate/model availability must not remain stale forever.
            if path.exists():
                return path.read_bytes()
            if negative.exists() and time.time() - negative.stat().st_mtime < 86400:
                return b''
            for attempt in range(4):
                delay = 2 ** attempt
                if self.cancel.is_set():
                    raise InterruptedError('Build cancelled.')
                try:
                    request = Request(url, headers={'User-Agent': 'Trimmings-Pfam-builder/1.0', 'Accept': 'application/json,*/*'})
                    self.wait_for_host(urlparse(url).hostname)
                    with urlopen(request, timeout=30) as response:
                        data = response.read(150_000_001)
                    if len(data) > 150_000_000:
                        raise ValueError('Public resource exceeds 150 MB: ' + url)
                    if data.lstrip().lower().startswith((b'<!doctype html', b'<html')):
                        raise ValueError('Public service returned HTML instead of data: ' + url)
                    # Unique temporaries are also safe when separate builds share a cache.
                    temporary = None
                    try:
                        with tempfile.NamedTemporaryFile(dir=self.folder, delete=False) as stream:
                            temporary = Path(stream.name)
                            stream.write(data)
                        temporary.replace(path)
                    finally:
                        if temporary is not None:
                            temporary.unlink(missing_ok=True)
                    return data
                except HTTPError as error:
                    status = error.code
                    retry_after = error.headers.get('Retry-After', '0') if error.headers else '0'
                    error.close()
                    if status in (404, 204):
                        negative.touch()
                        return b''
                    try:
                        delay = max(delay, float(retry_after))
                    except ValueError:
                        try:
                            delay = max(delay, parsedate_to_datetime(retry_after).timestamp() - time.time())
                        except (ValueError, TypeError, OverflowError):
                            pass
                    if status in (429, 503):
                        with self.host_lock:
                            host = urlparse(url).hostname
                            self.host_cooldown[host] = max(self.host_cooldown.get(host, 0), time.monotonic() + delay)
                    if status not in (429, 500, 502, 503, 504) or attempt == 3:
                        raise
                except (URLError, TimeoutError, HTTPException, ConnectionError):
                    if attempt == 3:
                        raise
                if self.cancel.wait(delay):
                    raise InterruptedError('Build cancelled.')
        raise RuntimeError('Public resource download failed.')

    def json(self, url):
        data = self.get(url)
        try:
            return json.loads(data) if data else None
        except ValueError:
            # A transient bad JSON response must not poison a resumable build.
            (self.folder / hashlib.sha256(url.encode()).hexdigest()).unlink(missing_ok=True)
            raise


def polymer_chains(data, chain_id=''):
    import gemmi
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(data.decode()).sole_block())
    if not len(structure):
        raise ValueError('Coordinate file contains no structural models.')
    for chain in structure[0]:
        if chain_id and chain.name != chain_id:
            continue
        for subchain in dict.fromkeys(r.subchain for r in chain if r.label_seq is not None):
            entity = structure.get_entity_of(chain.get_subchain(subchain))
            if entity is None or not entity.full_sequence:
                continue
            sequence = gemmi.one_letter_code(entity.full_sequence).upper()
            residues = [r for r in chain if r.subchain == subchain and r.label_seq is not None and r.find_atom('CA', '*')]
            yield sequence, residues, chain.name


def crop_positions(sequence, residues, first, last, args, plddt=False):
    if not 1 <= first <= last <= len(sequence):
        raise ValueError('Reference domain positions exceed the source polymer sequence.')
    selected = sorted((r for r in residues if first <= r.label_seq <= last), key=lambda r: r.label_seq)
    if len({r.label_seq for r in selected}) != len(selected):
        raise ValueError('Duplicate polymer residue positions.')
    if len(selected) / (last - first + 1) < args.min_coverage:
        raise ValueError('Insufficient modeled domain coverage.')
    if any(fs.residue_letter(r) not in (sequence[r.label_seq - 1], 'X') for r in selected):
        raise ValueError('Coordinate/sequence mismatch.')
    mean = sum(r.find_atom('CA', '*').b_iso for r in selected) / max(1, len(selected))
    if plddt and mean < args.min_plddt:
        raise ValueError(f'Mean domain pLDDT {mean:.1f} is below {args.min_plddt:g}.')
    return selected, mean


def family_profile(library, family, folder, args, cancel):
    path = folder / (family['accession'] + '.hmm')
    if not path.exists():
        partial = path.with_suffix('.part')
        # hmmfetch supports -o; never interpolate family names into a shell.
        run([args.hmmfetch, '-o', partial, library, family['version']], folder / 'hmmfetch.log', cancel)
        partial.replace(path)
    return path


def profile_hits(sequence, family, profile, folder, args, cancel):
    fasta, table = folder / 'chain.fasta', folder / 'chain.tsv'
    fasta.write_text('>chain\n' + sequence + '\n')
    run([args.hmmsearch, '--cut_ga', '--noali', '--cpu', '1', '--domtblout', table, profile, fasta],
        folder / 'chain.log', cancel)
    return sorted(parse_domtbl(table, {family['accession']: family}, args.min_coverage), key=lambda row: (row[5], -row[4]))


def profile_crop(data, chain, family, profile, folder, args, cancel, plddt=False):
    for sequence, residues, author_chain in polymer_chains(data, chain):
        hits = profile_hits(sequence, family, profile, folder, args, cancel)
        for _, _, first, last, score, evalue, coverage in hits:
            try:
                crop, _ = crop_positions(sequence, residues, first, last, args, plddt=plddt)
                return crop, author_chain, (first, last), {'bitscore': score, 'evalue': evalue, 'hmm_coverage': coverage}
            except ValueError:
                continue
    raise ValueError('No current Pfam gathering-threshold match with sufficient modeled coverage.')


def alphafold_seed(uniprot, seed, family, cache, work, profile, folder, args, cancel):
    predictions = cache.json(f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot}') or []
    if not isinstance(predictions, list) or any(not isinstance(item, dict) for item in predictions):
        raise ValueError('Unrecognized AlphaFold prediction response.')
    for prediction in predictions:
        url = prediction.get('cifUrl', '')
        if prediction.get('uniprotAccession') != uniprot or urlparse(url).hostname != 'alphafold.ebi.ac.uk':
            continue
        data = cache.get(url)
        for sequence, residues, chain in polymer_chains(data):
            first = sequence.find(seed)
            if first < 0:
                continue
            try:
                crop, chain, segment, match = profile_crop(data, chain, family, profile, folder, args, cancel, plddt=True)
                mean = sum(r.find_atom('CA', '*').b_iso for r in crop) / len(crop)
                return save_reference(family, crop, work, 'alphafold', args, catalog='Pfam', pdb=None,
                    uniprot=uniprot, chain=chain, segments=[segment], mean_plddt=mean, profile_match=match,
                    structure_url=f'https://alphafold.ebi.ac.uk/entry/{uniprot}', structure_title=uniprot,
                    source_coordinate_url=url, source_coordinate_sha256=hashlib.sha256(data).hexdigest(),
                    evidence='Pfam gathering-threshold match on an AlphaFold polymer containing an exact current Pfam seed/member fragment',
                    source_sequence_sha256=hashlib.sha256(seed.encode()).hexdigest())
            except ValueError:
                continue
    raise ValueError('No sequence-verified AlphaFold profile crop meeting coverage/confidence requirements.')


def interpro_candidates(accession, cache):
    """Exhaust pages only for still-uncovered families; never silently sample."""
    url = f'{fs.INTERPRO}/protein/uniprot/entry/pfam/{accession}/?page_size=200&extra_fields=sequence'
    seen = set()
    while url:
        parsed = urlparse(url)
        if parsed.scheme != 'https' or parsed.hostname != 'www.ebi.ac.uk' or url in seen:
            raise ValueError('Invalid/cyclic InterPro pagination URL.')
        seen.add(url)
        data = cache.json(url)
        if data is None:
            return
        if not isinstance(data, dict) or not isinstance(data.get('results'), list):
            raise ValueError('Unrecognized InterPro protein response.')
        for item in data['results']:
            if not isinstance(item, dict) or not isinstance(item.get('metadata'), dict):
                raise ValueError('Invalid InterPro protein metadata.')
            metadata = item.get('metadata', {})
            if metadata.get('in_alphafold') is False:
                continue
            identifier = metadata.get('accession', '')
            sequence = metadata.get('sequence') or item.get('extra_fields', {}).get('sequence')
            if not UNIPROT.fullmatch(identifier):
                continue
            if isinstance(sequence, dict):
                sequence = sequence.get('sequence')
            if not isinstance(sequence, str) or not sequence:
                continue
            for entry in item.get('entries', []):
                if entry.get('accession') != accession:
                    continue
                for location in entry.get('entry_protein_locations', []):
                    fragments = location.get('fragments', [])
                    if fragments:
                        # Keep the natural intervening sequence, never concatenate
                        # distant fragments. The pinned HMM validates the resulting span.
                        first = min(fragment['start'] for fragment in fragments)
                        last = max(fragment['end'] for fragment in fragments)
                        if 1 <= first <= last <= len(sequence):
                            yield identifier, sequence[first - 1:last]
        url = data.get('next')


def predict_seed(seed, family, work, profile, args, cancel):
    import gemmi
    command_digest = hashlib.sha256(args.predict_command.encode()).hexdigest()[:12]
    folder = work / 'predictions' / family['accession'] / (hashlib.sha256(seed.encode()).hexdigest()[:16] + '-' + command_digest)
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / 'input.fasta'
    fasta.write_text('>' + family['accession'] + '\n' + seed + '\n')
    hits = profile_hits(seed, family, profile, folder, args, cancel)
    if not hits:
        raise ValueError('Seed lacks sufficient current-profile coverage at gathering thresholds.')
    paths = [path for path in folder.rglob('*') if path.suffix in ('.cif', '.pdb')]
    if not paths:
        values = {'fasta': str(fasta), 'outdir': str(folder), 'accession': family['accession']}
        command = [token.format(**values) for token in shlex.split(args.predict_command)]
        run(command, folder / 'predict.log', cancel, args.prediction_timeout)
        paths = [path for path in folder.rglob('*') if path.suffix in ('.cif', '.pdb')]
    for path in sorted(paths):
        structure = gemmi.read_structure(str(path))
        for chain in structure[0]:
            residues = [r for r in chain if r.find_atom('CA', '*')]
            sequence = ''.join(fs.residue_letter(r) for r in residues)
            if sequence != seed:
                continue
            for _, _, first, last, score, evalue, coverage in hits:
                crop = residues[first - 1:last]
                # Predictor adapters must preserve pLDDT in the C-alpha B-factor field.
                mean = sum(r.find_atom('CA', '*').b_iso for r in crop) / max(1, len(crop))
                if mean < args.min_plddt:
                    continue
                return save_reference(family, crop, work, 'prediction', args, catalog='Pfam', pdb=None,
                    chain=chain.name, structure_title=family['accession'] + ' predicted seed',
                    structure_url=f'https://www.ebi.ac.uk/interpro/entry/pfam/{family["accession"]}/',
                    prediction_command=args.predict_command, coordinate_sha256=checksum(path, 'sha256'),
                    source_sequence_sha256=hashlib.sha256(seed.encode()).hexdigest(),
                    segments=[(first, last)], mean_plddt=mean,
                    profile_match={'bitscore': score, 'evalue': evalue, 'hmm_coverage': coverage},
                    evidence='Current Pfam profile match on a locally predicted, exact Pfam seed sequence')
    raise ValueError('Predictor produced no sequence-exact model with sufficient pLDDT.')


def fill_family(accession, state_path, families, library, cache, work, args, cancel):
    family, records, errors = families[accession], [], []
    folder = work / 'families' / accession
    folder.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(f'file:{state_path}?mode=ro', uri=True)) as state:
        existing = {row[0] for row in state.execute('SELECT sequence FROM refs WHERE family=?', (accession,))}
        pdbs = state.execute('SELECT pdb,chain FROM pdbs WHERE family=? ORDER BY pdb,chain', (accession,)).fetchall()
        seeds = state.execute('SELECT uniprot,sequence FROM seeds WHERE family=? ORDER BY rowid', (accession,)).fetchall()
    def enough(item):
        if item[0]['sequence'] not in existing:
            records.append(item)
            existing.add(item[0]['sequence'])
        return len(existing) >= args.representatives
    def error(source, identifier, exception):
        if isinstance(exception, InterruptedError):
            raise exception
        retryable = (isinstance(exception, (URLError, TimeoutError, HTTPException, ConnectionError))
                     or isinstance(exception, HTTPError) and exception.code in (429, 500, 502, 503, 504))
        # HTTPError subclasses URLError; permanent HTTP failures must not become
        # retryable merely because of their base class.
        if isinstance(exception, HTTPError):
            retryable = exception.code in (429, 500, 502, 503, 504)
        detail = {'source': source, 'id': identifier, 'error': str(exception)}
        if retryable:
            detail['retryable'] = True
        with (folder / 'attempts.jsonl').open('a') as stream:
            stream.write(json.dumps(detail) + '\n')
        if len(errors) < 25:
            errors.append(detail)
        elif retryable:
            errors[-1] = detail
        return retryable
    profile = None
    def af(uniprot, sequence):
        nonlocal profile
        if profile is None:
            profile = family_profile(library, family, folder, args, cancel)
        return alphafold_seed(uniprot, sequence, family, cache, work, profile, folder, args, cancel)
    for pdb, chain in pdbs:
        try:
            if profile is None:
                profile = family_profile(library, family, folder, args, cancel)
            data = cache.get(f'https://files.rcsb.org/download/{pdb}.cif.gz')
            if not data:
                continue
            residues, author_chain, segment, match = profile_crop(data, chain, family, profile, folder, args, cancel)
            item = save_reference(family, residues, work, 'pdb', args, catalog='Pfam', pdb=pdb,
                chain=author_chain, segments=[segment], profile_match=match,
                structure_title=pdb, structure_url=f'https://www.rcsb.org/structure/{pdb}',
                source_coordinate_url=f'https://files.rcsb.org/download/{pdb}.cif.gz',
                source_coordinate_sha256=hashlib.sha256(data).hexdigest(),
                evidence='Current Pfam gathering-threshold match on a SIFTS-linked PDB polymer chain')
            if enough(item):
                return records, errors
        except (OSError, ValueError, RuntimeError, KeyError) as exception:
            if error('pdb', pdb, exception):
                return records, errors
    attempted = set()
    for uniprot, sequence in seeds:
        if cancel.is_set():
            raise InterruptedError('Build cancelled.')
        attempted.add((uniprot, sequence))
        try:
            if enough(af(uniprot, sequence)):
                return records, errors
        except (OSError, ValueError, RuntimeError, KeyError) as exception:
            if error('alphafold-seed', uniprot, exception):
                return records, errors
    if not args.no_interpro:
        try:
            for uniprot, sequence in interpro_candidates(accession, cache):
                if (uniprot, sequence) in attempted:
                    continue
                attempted.add((uniprot, sequence))
                try:
                    if enough(af(uniprot, sequence)):
                        return records, errors
                except (OSError, ValueError, RuntimeError, KeyError) as exception:
                    if error('alphafold-member', uniprot, exception):
                        return records, errors
        except (OSError, ValueError, RuntimeError, KeyError) as exception:
            error('interpro', accession, exception)
    if args.predict_command and seeds:
        # Additional seeds are tried if the first representative cannot be modeled.
        for _, sequence in seeds:
            if sequence in existing:
                continue
            try:
                if profile is None:
                    profile = family_profile(library, family, folder, args, cancel)
                with args._prediction_gate:
                    if enough(predict_seed(sequence, family, work, profile, args, cancel)):
                        return records, errors
            except (OSError, ValueError, RuntimeError, KeyError) as exception:
                error('prediction', accession, exception)
    return records, errors


def coverage_report(output, state, families, selected, args):
    counts = represented(state)
    covered = families.keys() & counts.keys()
    complete = len(covered) == len(families) and selected == families.keys()
    runs = {family: (status, json.loads(errors)) for family, status, errors in state.execute('SELECT * FROM runs')}
    summary = {'pfam_release': args.release, 'families': len(families), 'selected_families': len(selected),
               'families_with_verified_structures': len(covered), 'missing_families': len(families) - len(covered),
               'structural_coverage_complete': complete, 'profile_coverage_complete': True,
               'references': state.execute('SELECT count(*) FROM refs').fetchone()[0],
               'criteria': {'minimum_modeled_and_profile_coverage': args.min_coverage,
                            'minimum_predicted_mean_plddt': args.min_plddt, 'minimum_modeled_residues': args.min_residues},
               'by_type': {kind: {'total': sum(f['type'] == kind for f in families.values()),
                                  'covered': sum(families[a]['type'] == kind for a in covered)}
                           for kind in sorted({f['type'] for f in families.values()})}}
    if getattr(args, 'scope', 'all') == 'folded':
        required = {a for a, family in families.items() if family['type'] in FOLDED_PROFILE_TYPES}
        excluded = set(families) - required
        summary.update(target_scope='folded', required_families=len(required),
                       covered_required_families=len(required & covered),
                       missing_required_families=len(required - covered),
                       target_scope_complete=required <= covered and required <= selected,
                       excluded_family_accessions=sorted(excluded),
                       excluded_by_type=dict(Counter(families[a]['type'] for a in excluded)))
    atomic_json(Path(output) / 'build-coverage.json', summary)
    with (Path(output) / 'build-families.tsv').open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['accession', 'version', 'name', 'type', 'references', 'status', 'failure_sample'])
        for accession, family in sorted(families.items()):
            status, errors = runs.get(accession, ('unattempted' if accession in selected else 'outside requested subset', []))
            if getattr(args, 'scope', 'all') == 'folded' and family['type'] not in FOLDED_PROFILE_TYPES:
                status = 'outside folded-domain structural scope; retained for sequence annotation'
            elif accession in counts:
                status = 'verified structural representative'
            writer.writerow([accession, family['version'], family['name'], family['type'], counts.get(accession, 0),
                             status, json.dumps(errors)])
    return summary


def publish(output, state, coverage, profile_library, assets, args, cancel):
    complete = coverage.get('target_scope_complete', coverage['structural_coverage_complete'])
    if not complete and not args.publish_partial:
        missing = coverage.get('missing_required_families', coverage['missing_families'])
        raise IncompleteLibrary(f'{missing} Pfam families in the requested scope still lack verified coordinates; '
                                'no incomplete release was published. See build-families.tsv. '
                                'Rerun to resume, optionally with --predict-command.')
    if not coverage['references']:
        raise IncompleteLibrary('No verified structural references; no database published.')
    with staged_generation(output) as generation:
        inputs = generation / 'inputs'
        inputs.mkdir()
        metadata = generation / 'references.sqlite'
        with closing(sqlite3.connect(metadata)) as destination:
            destination.execute('CREATE TABLE refs(id TEXT PRIMARY KEY,record TEXT)')
            excluded = set(coverage.get('excluded_family_accessions', []))
            for identifier, family, record, path, digest in state.execute('SELECT id,family,record,path,sha256 FROM refs ORDER BY id'):
                if family in excluded:
                    continue
                if checksum(path, 'sha256') != digest:
                    raise ValueError('Reference coordinates changed before publication: ' + identifier)
                target = inputs / (identifier + '.pdb')
                try:
                    os.link(path, target)
                except OSError:
                    shutil.copyfile(path, target)
                # Older checkpoints retain the coordinate source's catalog in
                # this field, although the exported assignment is now Pfam.
                record = pfam_reference_record(json.loads(record))
                destination.execute('INSERT INTO refs VALUES(?,?)', (identifier, json.dumps(record)))
            destination.commit()
        database = generation / 'foldseek_db'
        run([args.foldseek, 'createdb', inputs, database, '--threads', args.threads, '--chain-name-mode', '0'],
            generation / 'createdb.log', cancel)
        with Path(str(database) + '.lookup').open() as stream:
            targets = {row[1].removesuffix('.pdb') for row in csv.reader(stream, delimiter='\t')}
        expected = {identifier for identifier, family in state.execute('SELECT id,family FROM refs') if family not in excluded}
        if not expected:
            raise IncompleteLibrary('No verified structural references in the requested scope.')
        if targets != expected:
            raise ValueError('Foldseek dropped/renamed structural references; publication refused.')
        shutil.rmtree(inputs)
        shutil.copyfile(Path(output) / 'build-families.tsv', generation / 'families.tsv')
        shutil.copyfile(assets / 'md5_checksums', generation / 'upstream-md5.txt')
        provenance = Path(state.execute('PRAGMA database_list').fetchone()[2]).parent / 'provenance.json'
        if provenance.is_file():
            shutil.copyfile(provenance, generation / 'build-provenance.json')
        atomic_json(generation / 'coverage.json', coverage)
        manifest = {'schema': 2, 'catalog': 'Pfam', 'release': args.release,
                    'pilot': not complete,
                    'description': 'Release-pinned, sequence-validated Pfam structural reference library.',
                    'database': str(database.relative_to(output)), 'reference_metadata': str(metadata.relative_to(output)),
                    'profile_database': str(profile_library.relative_to(output)), 'reference_count': len(targets),
                    'coverage': coverage}
        if coverage.get('target_scope') == 'folded':
            manifest.update(scope='folded', description='Sequence-validated Pfam Domain, Family and Repeat structural templates; motif, disordered and coiled-coil types excluded.')
        atomic_json(Path(output) / 'coverage.json', coverage)
        atomic_json(Path(output) / 'manifest.json', manifest)
    return manifest


def distribution_archive(output, manifest, destination, *, prefix=None):
    """Package published data and notices, omitting build state and downloads."""
    output, destination = Path(output).resolve(), Path(destination).expanduser().resolve()
    if not str(destination).endswith(('.tar.gz', '.tgz')):
        raise ValueError('Distribution archives must use .tar.gz or .tgz.')
    generation = fs.library_file(output, manifest['database']).parent
    profile = fs.library_file(output, manifest['profile_database']) if manifest.get('profile_database') else None
    files = {output / 'manifest.json', output / 'coverage.json'}
    if manifest.get('reference_metadata'):
        reference_metadata = fs.library_file(output, manifest['reference_metadata'])
        files.add(reference_metadata)
        if (reference_metadata.parent / 'source-manifest.json').is_file():
            files.add(reference_metadata.parent / 'source-manifest.json')
    for name in ('audit', 'families'):
        if manifest.get('pfam_mapping', {}).get(name):
            files.add(fs.library_file(output, manifest['pfam_mapping'][name]))
    files.update(path for path in generation.rglob('*') if path.is_file())
    if profile is not None:
        files.update(Path(str(profile) + suffix) for suffix in ('', '.ssi', '.h3f', '.h3i', '.h3m', '.h3p'))
        files.update(profile.parent / name for name in ('Pfam.version.gz', 'Pfam-A.hmm.dat.gz', 'profile-manifest.json'))
    notices = [Path(__file__).parent / name for name in ('LICENSE', 'DATA_NOTICE.md', 'THIRD_PARTY.md')]
    if any(not path.is_file() or path.is_symlink() for path in notices):
        raise ValueError('Distribution requires LICENSE, DATA_NOTICE.md and THIRD_PARTY.md beside the builder.')
    if any(path.is_symlink() or not path.resolve().is_relative_to(output) for path in files):
        raise ValueError('Distribution input escapes the output library.')
    if destination in files or destination.is_relative_to(generation) or (profile and destination.is_relative_to(profile.parent)):
        raise ValueError('The archive must not replace a published library file.')
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.part')
    release = manifest.get('pfam_release', manifest.get('release'))
    if profile is not None and not release:
        raise ValueError('Distribution manifest is missing the Pfam release.')
    if prefix is None:
        prefix = ('pfam-' + release + ('-folded' if manifest.get('scope') == 'folded' else '')
                  + ('-pilot' if manifest['pilot'] else '')) if profile else 'domain-library'
    if not prefix or Path(prefix).name != prefix or prefix in ('.', '..'):
        raise ValueError('Distribution prefix must be one directory name.')
    say('Packaging verified release: ' + str(destination))
    with tarfile.open(partial, 'w:gz', compresslevel=3) as archive:
        for path in notices:
            archive.add(path, arcname=str(Path(prefix) / path.name), recursive=False)
        for path in sorted(files):
            archive.add(path, arcname=str(Path(prefix) / path.relative_to(output)), recursive=False)
    partial.replace(destination)
    digest = checksum(destination, 'sha256')
    destination.with_name(destination.name + '.sha256').write_text(digest + '  ' + destination.name + '\n')
    return destination


def integration_request(args):
    """Identify a completed request independently of execution resource limits."""
    names = ('release', 'families', 'representatives', 'min_coverage', 'min_plddt',
             'min_residues', 'predict_command', 'no_interpro', 'publish_partial')
    return {**{name: getattr(args, name) for name in names},
            **({'scope': 'folded'} if getattr(args, 'scope', 'all') == 'folded' else {}),
            'reuse': [str(Path(path).expanduser().resolve()) for path in args.reuse],
            'target': str(Path(args.integrate_into).expanduser().resolve())}


def default_library(target, *, write=False):
    root = Path(__file__).resolve().parent / '.domain_data'
    if not target.is_relative_to(root):
        raise ValueError('The default library must be inside .domain_data beside the Trimmings script.')
    if write:
        atomic_json(root / 'default-library.json', {'library': str(target.relative_to(root))})


def completed_integration(output, args):
    """Do not scan and merge an already-integrated request a second time."""
    marker = output / 'integration.json'
    if not args.integrate_into or not marker.is_file():
        return None
    info = json.loads(marker.read_text())
    if info.get('request') != integration_request(args):
        return None
    target = Path(args.integrate_into).expanduser().resolve()
    if (checksum(output / 'manifest.json', 'sha256') != info.get('publication_sha256') or
            checksum(target / 'manifest.json', 'sha256') != info.get('manifest_sha256')):
        return None
    for source in info.get('input_manifests', []):
        if Path(source['path']) != target and checksum(Path(source['path']) / 'manifest.json', 'sha256') != source['sha256']:
            return None
    result = Library(output).manifest
    integrated = Library(target).manifest
    profile = fs.library_file(target, integrated['profile_database'])
    if not all(Path(str(profile) + suffix).is_file() for suffix in ('', '.ssi', '.h3f', '.h3i', '.h3m', '.h3p')):
        raise ValueError('The completed integration is missing its profile database files.')
    if args.make_default:
        default_library(target, write=True)
    if args.archive:
        archive = Path(args.archive).expanduser().resolve()
        digest = archive.with_name(archive.name + '.sha256')
        if (not archive.is_file() or not digest.is_file() or
                checksum(archive, 'sha256') != digest.read_text().split()[0]):
            distribution_archive(target, integrated, archive)
    say('This request is already integrated at ' + str(target))
    return result


def build(args):
    output = Path(args.out).expanduser().resolve()
    cancel = Event()
    if args.make_default:
        if not args.integrate_into:
            raise ValueError('--make-default requires --integrate-into.')
        default_library(Path(args.integrate_into).expanduser().resolve())
    if args.integrate_into:
        with build_lock(output):
            previous = completed_integration(output, args)
        if previous is not None:
            return previous
    libraries = [Library(path) for path in args.reuse]
    if any(library.folder == output for library in libraries):
        raise ValueError('Choose a separate output directory; --reuse libraries are read-only inputs.')
    if not args.audit_only:
        import gemmi
        for name in ('foldseek', 'hmmsearch', 'hmmfetch', 'hmmpress'):
            if not shutil.which(getattr(args, name)):
                raise ValueError(f'Missing executable {getattr(args, name)}; install Foldseek and HMMER 3.')
    with build_lock(output):
        assets, md5s, count = release_assets(output, args.release, args.assets, args.audit_only)
        families = catalog(assets / 'Pfam-A.hmm.dat.gz', count)
        selected = set(args.families) if args.families else set(families)
        if not selected or selected - families.keys():
            raise ValueError('Unknown/empty Pfam family selection.')
        if args.scope == 'folded':
            selected = {a for a in selected if families[a]['type'] in FOLDED_PROFILE_TYPES}
            if not selected:
                raise ValueError('The selection contains no folded-domain candidate types.')
        if args.audit_only:
            result = mapping_audit(output, families, libraries, args.release)
            say(json.dumps(result))
            return result
        tools = {'foldseek': subprocess.check_output([args.foldseek, 'version'], text=True).strip(),
                 'hmmer': subprocess.check_output([args.hmmsearch, '-h'], text=True).splitlines()[:2],
                 'gemmi': gemmi.__version__}
        fingerprint = hashlib.sha256(json.dumps([VERSION, md5s, [(str(lib.folder), lib.digest) for lib in libraries],
            sorted(selected), args.min_coverage, args.min_plddt, args.min_residues, args.shard_size,
            args.reference_batch_size, args.representatives, tools],
            sort_keys=True).encode()).hexdigest()
        work = output / 'work' / fingerprint[:24]
        work.mkdir(parents=True, exist_ok=True)
        atomic_json(work / 'provenance.json', {'tools': tools, 'upstream_md5': md5s,
                    'reuse': [{'path': str(lib.folder), 'manifest_sha256': lib.digest} for lib in libraries]})
        state_path = work / 'state.sqlite'
        with closing(state_database(state_path)) as state:
            say('Validating saved reference checkpoints')
            # Lost/changed files are not valid checkpoints.
            for identifier, path, digest in state.execute('SELECT id,path,sha256 FROM refs').fetchall():
                if not Path(path).is_file() or checksum(path, 'sha256') != digest:
                    state.execute('DELETE FROM refs WHERE id=?', (identifier,))
                    state.execute('UPDATE hits SET tried=0')
            state.commit()
            say(f'Preparing the complete {len(families)}-model Pfam profile database')
            profile_library, shards = profiles(assets, output / 'profiles' / ('Pfam' + args.release), families, selected, args.shard_size)
            index_marker = profile_library.parent / 'indexed.json'
            index_state = json.loads(index_marker.read_text()) if index_marker.exists() else {}
            ssi = Path(str(profile_library) + '.ssi')
            if not ssi.exists() or checksum(ssi, 'sha256') != index_state.get('ssi_sha256'):
                say('Indexing profiles for individual-family retrieval')
                Path(str(profile_library) + '.ssi').unlink(missing_ok=True)
                run([args.hmmfetch, '--index', profile_library], work / 'hmm-index.log', cancel)
                index_state['ssi_sha256'] = checksum(ssi, 'sha256')
                atomic_json(index_marker, index_state)
            pressed = [Path(str(profile_library) + suffix) for suffix in ('.h3f', '.h3i', '.h3m', '.h3p')]
            if not all(path.is_file() and checksum(path, 'sha256') == index_state.get('pressed_sha256', {}).get(path.name)
                       for path in pressed):
                say('Preparing the complete hmmscan search index')
                run([args.hmmpress, '-f', profile_library], work / 'hmm-press.log', cancel)
                index_state['pressed_sha256'] = {path.name: checksum(path, 'sha256') for path in pressed}
                atomic_json(index_marker, index_state)
            shutil.copyfile(assets / 'Pfam.version.gz', profile_library.parent / 'Pfam.version.gz')
            shutil.copyfile(assets / 'Pfam-A.hmm.dat.gz', profile_library.parent / 'Pfam-A.hmm.dat.gz')
            atomic_json(profile_library.parent / 'profile-manifest.json',
                        {'pfam_release': args.release, 'models': len(families), 'profile_coverage_complete': True,
                         'hmm_compressed_md5': md5s['Pfam-A.hmm.gz']})
            args._prediction_gate = BoundedSemaphore(args.prediction_workers)
            try:
                if libraries:
                    say('Scanning existing reference sequences with current Pfam gathering thresholds')
                    fasta = export_references(state, libraries, work, args.foldseek, cancel)
                    adaptive_scan(state, libraries, profile_library, shards, work, families, selected, args, cancel)
                say('Indexing the full Pfam seed inventory and bulk PDB mappings')
                index_seeds(state, assets / 'Pfam-A.seed.gz', families)
                sifts = output / 'downloads' / 'pdb_chain_pfam.csv.gz'
                download(SIFTS, sifts)
                index_pdbs(state, sifts, families)
                cache = PublicCache(args.cache or output / 'downloads' / 'structures', cancel)
                counts = represented(state)
                pending = sorted(a for a in selected if counts.get(a, 0) < args.representatives)
                say(f'Filling {len(pending)} families from PDB/AlphaFold; {len(counts)} already represented')
                started, completed = time.monotonic(), 0
                def fill(accession):
                    return fill_family(accession, state_path, families, profile_library, cache, work, args, cancel)
                for accession, (records, errors) in parallel_jobs(pending, fill, args.workers, cancel):
                    for item in records:
                        if add_reference(state, item):
                            counts[accession] = counts.get(accession, 0) + 1
                    status = ('retryable source failure' if any(e.get('retryable') for e in errors)
                              else 'no suitable public coordinates') if not records else 'verified structural representative'
                    state.execute('INSERT OR REPLACE INTO runs VALUES(?,?,?)', (accession, status, json.dumps(errors)))
                    state.commit()
                    completed += 1
                    say(f'[{completed}/{len(pending)}] {accession}: {counts.get(accession, 0)} references; '
                        f'{len(counts)}/{len(families)} families covered; {time.monotonic() - started:.0f}s')
                    atomic_json(output / 'progress.json', {'phase': 'fill', 'completed': completed, 'pending': len(pending),
                                'covered_families': len(counts), 'total_families': len(families)})
            except BaseException:
                cancel.set()
                coverage_report(output, state, families, selected, args)
                raise
            report = coverage_report(output, state, families, selected, args)
            result = publish(output, state, report, profile_library, assets, args, cancel)
            archive_root, archive_manifest = output, result
            if args.integrate_into:
                from domain_foldseek_ecod import merge
                target = Path(args.integrate_into).expanduser().resolve()
                if target == output:
                    raise ValueError('The Pfam build and integrated library need separate directories.')
                say('Integrating verified Pfam references into ' + str(target))
                integrated = merge([target, output], target, executable=args.foldseek, replace=True)
                archive_root, archive_manifest = target, integrated
                integration = {'library': str(target), 'reference_count': integrated['reference_count'],
                               'pfam_coverage': report,
                               'manifest_sha256': checksum(target / 'manifest.json', 'sha256'),
                               'publication_sha256': checksum(output / 'manifest.json', 'sha256'),
                               'request': integration_request(args),
                               'input_manifests': [{'path': str(lib.folder), 'sha256': lib.digest} for lib in libraries]}
                atomic_json(output / 'integration.json', integration)
                if args.make_default:
                    default_library(target, write=True)
            if args.archive:
                distribution_archive(archive_root, archive_manifest, args.archive)
            atomic_json(output / 'progress.json', {'phase': 'published', **report})
            say(f'Published {report["families_with_verified_structures"]}/{len(families)} families; '
                f'{report["references"]} references at {output}')
            return result


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', required=True, help='Separate output library directory; the same command resumes')
    p.add_argument('--release', default='39.0', help='Pinned official Pfam release (default: 39.0)')
    p.add_argument('--scope', choices=('all', 'folded'), default='all',
                   help='Structural completeness target; folded includes Domain, Family and Repeat types')
    p.add_argument('--reuse', action='append', default=[], metavar='LIBRARY', help='Existing Trimmings Foldseek library; repeatable')
    p.add_argument('--assets', help='Optional local directory of release files and md5_checksums')
    p.add_argument('--cache', help='Shared coordinate/API download cache')
    p.add_argument('--threads', type=int, default=max(1, (os.cpu_count() or 2) - 1), help='Total HMMER/Foldseek CPU budget')
    p.add_argument('--workers', type=int, default=4, help='Concurrent family/download workers (default: 4)')
    p.add_argument('--scan-jobs', type=int, default=2, help='Concurrent profile-scan shards, sharing --threads')
    p.add_argument('--shard-size', type=int, default=250, help='Profiles per resumable HMMER checkpoint')
    p.add_argument('--reference-batch-size', type=int, default=10000, help='Initial diverse reference batch; later batches double, without sampling')
    p.add_argument('--representatives', type=int, default=1, help='Desired distinct crop sequences per family; completeness requires at least one')
    p.add_argument('--min-coverage', type=float, default=.8, help='Minimum modeled and profile coverage (default: 0.8)')
    p.add_argument('--min-plddt', type=float, default=70, help='Minimum mean predicted-domain pLDDT (default: 70)')
    p.add_argument('--min-residues', type=int, default=20, help='Minimum modeled C-alpha residues (default: 20)')
    p.add_argument('--predict-command', help='Optional predictor adapter argv template using {fasta}, {outdir}, {accession}; no shell')
    p.add_argument('--prediction-timeout', type=float, default=3600, help='Seconds allowed per optional predictor invocation')
    p.add_argument('--prediction-workers', type=int, default=1, help='Concurrent local predictors; default 1 avoids GPU oversubscription')
    p.add_argument('--archive', help='Also write a distributable .tar.gz and SHA256, omitting downloads/work files')
    p.add_argument('--integrate-into', help='After publication, atomically merge verified references into this existing library')
    p.add_argument('--make-default', action='store_true', help='Use the integrated .domain_data library by default in ChimeraX')
    p.add_argument('--no-interpro', action='store_true', help='Skip exhaustive non-seed member fallback; gaps still fail completeness')
    p.add_argument('--families', nargs='+', help='Explicit subset for pilot/testing; cannot claim full-release completeness')
    p.add_argument('--publish-partial', action='store_true', help='Explicitly publish a pilot despite gaps; default refuses publication')
    p.add_argument('--audit-only', action='store_true', help='Only download the small catalog and audit existing family mappings')
    for name in ('foldseek', 'hmmsearch', 'hmmfetch', 'hmmpress'):
        p.add_argument('--' + name, default=name, help='Executable path')
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if not re.fullmatch(r'\d+\.\d+', args.release):
        p.error('Use a numeric pinned release such as 39.0.')
    if min(args.threads, args.workers, args.scan_jobs, args.shard_size, args.representatives, args.min_residues,
           args.prediction_workers, args.reference_batch_size) < 1:
        p.error('CPU, worker, shard, representative and residue counts must be positive.')
    if not 20 <= args.min_residues <= 9999:
        p.error('Trimmings structural searches require at least 20 paired residues; the PDB crop limit is 9999.')
    if (not 0 < args.min_coverage <= 1 or not 0 <= args.min_plddt <= 100
            or not math.isfinite(args.prediction_timeout) or args.prediction_timeout <= 0):
        p.error('Invalid coverage, pLDDT, or predictor timeout.')
    args.scan_jobs = min(args.scan_jobs, args.threads)
    args.workers = min(args.workers, args.threads)
    try:
        build(args)
    except IncompleteLibrary as error:
        say(str(error))
        return 2
    except KeyboardInterrupt:
        say('Interrupted. Completed shards, downloads, and family references are saved; rerun to resume.')
        return 130
    except (OSError, ValueError, RuntimeError, sqlite3.Error, KeyError, ImportError, subprocess.SubprocessError) as error:
        say(f'Build failed: {error}')
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
