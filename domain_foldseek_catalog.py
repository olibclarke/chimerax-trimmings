"""Build a resumable, CDD-wide structural reference library from public seed MSAs.

Only the build commands import this module. ChimeraX searches need neither
Gemmi nor the seed archive. Full downloaded structures remain in the cache;
the published reference database uses cropped C-alpha coordinates.
"""

from collections import Counter, defaultdict
from contextlib import closing, contextmanager
import hashlib
import json
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tarfile
import tempfile
import time
from urllib.error import URLError
from uuid import uuid4

import domain_foldseek as fs

SEED_URL = 'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz'
BUILD_VERSION = 1
GENERAL_CONTROLS = (
    ('P68871', 'hemoglobin beta', ['pfam00042']),
    ('P01112', 'H-Ras', ['cd00876', 'smart00010']),
    ('P0AA25', 'thioredoxin', ['cd01659']),
    ('P0DP23', 'calmodulin', ['cd00051', 'pfam00036']),
    ('P02794', 'ferritin heavy chain', ['cd00657']),
    ('P00760', 'trypsin', ['smart00020']),
    ('P12931', 'Src; SH2 and SH3 controls', ['pfam00017', 'pfam00018', 'cd00174']),
    ('P00698', 'lysozyme; unrepresented fold control', []),
)


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


@contextmanager
def build_lock(directory):
    import fcntl

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / '.build.lock').open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another builder is already using ' + str(directory)) from None
        yield


def fasta_records(text):
    header, sequence = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(sequence).replace('-', '').replace('.', '').upper()
            header, sequence = line[1:].strip(), []
        elif header is not None:
            sequence.append(''.join(line.split()))
    if header is not None:
        yield header, ''.join(sequence).replace('-', '').replace('.', '').upper()


def seed_identifier(header):
    match = re.search(r'(?:^|\|)pdb\|([A-Za-z0-9]{4})\|([^\s|]+)', header)
    if match:
        return 'pdb', match[1].lower(), match[2]
    match = re.search(r'(?:^|\|)(?:sp|tr)\|([A-Z0-9]+)(?:\.\d+)?\|', header)
    if match and re.fullmatch(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})', match[1]):
        return 'alphafold', match[1], ''
    return None


def catalog_summary(connection):
    total = connection.execute('SELECT count(*) FROM models WHERE source != "cl"').fetchone()[0]
    covered = connection.execute('SELECT count(DISTINCT model) FROM candidates').fetchone()[0]
    return {'models': total,
            'derived_superfamilies': connection.execute('SELECT count(*) FROM models WHERE source="cl"').fetchone()[0],
            'models_with_candidates': covered,
            'models_without_usable_direct_seed_candidates': total - covered,
            'models_by_source': dict(connection.execute('SELECT source,count(*) FROM models GROUP BY source')),
            'candidates_by_structure_source': dict(connection.execute('SELECT source,count(*) FROM candidates GROUP BY source'))}


def catalog(directory, archive, cache, per_source=8):
    """Index all sources without extracting archive members onto the filesystem."""
    directory, archive = Path(directory), Path(archive)
    if per_source < 1:
        raise ValueError('The candidate count must be positive.')
    with build_lock(directory):
        descriptions = fs.cdd_descriptions(cache)
        temporary = directory / 'catalog.sqlite.tmp'
        temporary.unlink(missing_ok=True)
        connection = sqlite3.connect(temporary)
        try:
            connection.executescript('''
                CREATE TABLE models(accession TEXT PRIMARY KEY, pssm_id TEXT, name TEXT, description TEXT, source TEXT);
                CREATE TABLE candidates(model TEXT, source TEXT, identifier TEXT, chain TEXT,
                                        sequence TEXT, header TEXT, digest TEXT,
                                        UNIQUE(model,source,identifier,chain,digest));
                CREATE INDEX candidate_models ON candidates(model);
                CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
            ''')
            connection.executemany('INSERT INTO models VALUES(?,?,?,?,?)',
                [(accession, label['pssm_id'], label['name'], label['description'], re.match(r'[A-Za-z]+', accession)[0])
                 for accession, label in descriptions.items()])
            seen_models, counts = set(), Counter()
            with tarfile.open(archive, 'r|gz') as stream:
                for member in stream:
                    accession = Path(member.name).name.removesuffix('.FASTA')
                    if not member.isfile() or accession not in descriptions or accession.startswith('cl'):
                        continue
                    if member.size > 50_000_000:
                        continue
                    text = stream.extractfile(member).read().decode('ascii')
                    seen_models.add(accession)
                    seen = set()
                    for header, sequence in fasta_records(text):
                        identifier = seed_identifier(header)
                        if (identifier is None or not 25 <= len(sequence) <= 9999
                                or not sequence.isalpha() or not sequence.isascii()):
                            continue
                        source, key, chain = identifier
                        digest = hashlib.sha256(sequence.encode()).hexdigest()
                        if (source, digest) in seen or counts[accession, source] >= per_source:
                            continue
                        seen.add((source, digest))
                        connection.execute('INSERT OR IGNORE INTO candidates VALUES(?,?,?,?,?,?,?)',
                                           (accession, source, key, chain, sequence, header, digest))
                        counts[accession, source] += 1
                    if len(seen_models) % 5000 == 0:
                        connection.commit()
                        print('Indexed', len(seen_models), 'CDD seed alignments', flush=True)
            with archive.open('rb') as stream:
                digest = hashlib.file_digest(stream, 'sha256').hexdigest()
            metadata = {'schema': 1, 'seed_url': SEED_URL, 'archive_sha256': digest,
                        'indexed_alignments': len(seen_models), 'candidates_per_source': per_source,
                        'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
            connection.executemany('INSERT INTO metadata VALUES(?,?)', [(k, json.dumps(v)) for k, v in metadata.items()])
            connection.commit()
            summary = {**metadata, **catalog_summary(connection)}
        finally:
            connection.close()
        temporary.replace(directory / 'catalog.sqlite')
        atomic_json(directory / 'coverage.json', summary)
    return summary


def benchmark(library, output, executable='foldseek'):
    """Small general-purpose control panel, not a calibrated accuracy estimate."""
    from contextlib import closing

    library, output = Path(library).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((library / 'manifest.json').read_text())
    if manifest['schema'] == 1:
        records = manifest['records']
    else:
        metadata = fs.library_file(library, manifest['reference_metadata'])
        with closing(sqlite3.connect(f'file:{metadata}?mode=ro', uri=True)) as connection:
            records = [json.loads(row[0]) for row in connection.execute('SELECT record FROM refs ORDER BY id')]
    report = {'library': str(library), 'library_manifest': manifest, 'controls': [],
              'note': 'Exact crop sequences contained in each query are excluded. This small panel is not a remote-homology or calibrated specificity benchmark.'}
    for accession, name, expected in GENERAL_CONTROLS:
        print('Control', accession, name, flush=True)
        try:
            predictions = fs.fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{accession}', output / 'downloads') or []
            prediction = next(p for p in predictions if p.get('uniprotAccession') == accession)
            query_sequence = prediction['uniprotSequence']
            url = prediction['cifUrl']
            if not url.startswith('https://alphafold.ebi.ac.uk/'):
                raise ValueError('Unexpected prediction source.')
            cif = output / (accession + '.cif')
            cif.write_bytes(fs.fetch(url, output / 'downloads'))
            references = [r for r in records if r.get('uniprot') != accession and r['sequence'] not in query_sequence]
            if not references:
                raise ValueError('No independent reference sequences remain.')
            with tempfile.TemporaryDirectory(prefix='cdd-control-') as directory:
                directory = Path(directory)
                inputs, query_dir = directory / 'references', directory / 'queries'
                inputs.mkdir(); query_dir.mkdir()
                for record in references:
                    shutil.copy2(library / 'domains' / (record['id'] + '.pdb'), inputs / (record['id'] + '.pdb'))
                subprocess.run([executable, 'createdb', str(inputs), str(directory / 'db'), '--threads', '2'],
                               check=True, stdout=subprocess.DEVNULL)
                atomic_json(directory / 'manifest.json', {'schema': 1, 'pilot': True, 'database': 'db', 'records': references})
                queries = fs.query_domains(cif, query_dir)
                result = fs.run_search(query_dir, queries, directory, executable=executable)
            matches = [{'cdd_accession': h['reference']['cdd_accession'], 'name': h['reference']['name'],
                        'start': min(h['aligned_positions']), 'end': max(h['aligned_positions']),
                        'evalue': h['evalue'], 'reference_coverage': h['tcov'], 'alignment_tm': h['alntmscore']}
                       for h in fs.select_matches(result['hits'])]
            found = {h['cdd_accession'] for h in matches}
            row = {'accession': accession, 'name': name, 'expected_models': expected,
                   'excluded_identical_references': len(records) - len(references),
                   'matches': matches, 'search_seconds': result['elapsed_seconds'],
                   'expected_family_found': bool(found.intersection(expected)) if expected else None}
        except (OSError, ValueError, KeyError, IndexError, StopIteration, RuntimeError) as error:
            row = {'accession': accession, 'name': name, 'error': str(error)}
        report['controls'].append(row)
        atomic_json(output / 'results.json', report)
        print(json.dumps(row), flush=True)
    return report


def crop_seed(cif, chain_id, sequence, minimum_plddt=None):
    """Require an exact, unique seed-fragment placement in a verified polymer chain."""
    import gemmi

    structure = gemmi.read_structure(str(cif))
    sequences = {e.name: ''.join(gemmi.one_letter_code([name]) for name in e.full_sequence)
                 for e in structure.entities if e.polymer_type != gemmi.PolymerType.Unknown}
    choices = []
    for chain in structure[0]:
        groups = defaultdict(list)
        for residue in chain:
            if residue.label_seq is not None:
                groups[residue.subchain].append(residue)
        # A PDB seed header specifies an author chain; AlphaFold may supply no chain.
        if chain_id and chain.name != chain_id:
            continue
        for subchain, residues in groups.items():
            entity = structure.get_entity_of(chain.get_subchain(subchain))
            parent = sequences.get(entity.name, '') if entity is not None else ''
            first = parent.find(sequence)
            if first < 0 or parent.find(sequence, first + 1) >= 0:
                continue
            selected = [r for r in residues if first < r.label_seq <= first + len(sequence) and r.find_atom('CA', '*')]
            selected.sort(key=lambda r: r.label_seq)
            if len(selected) < 25 or len(selected) / len(sequence) < .8:
                continue
            if any(fs.residue_letter(r) not in (parent[r.label_seq - 1], 'X') for r in selected):
                continue
            plddt = sum(r.find_atom('CA', '*').b_iso for r in selected) / len(selected)
            if minimum_plddt is not None and plddt < minimum_plddt:
                continue
            choices.append((selected, chain.name, (first + 1, first + len(sequence)), plddt))
    if not choices:
        raise ValueError('No exact seed placement with sufficient modeled coverage and confidence.')
    # Identical chains in a symmetric assembly are interchangeable references.
    return choices[0]


def candidate_record(candidate, label, cache, structures, minimum_plddt=70):
    model, source, identifier, chain, sequence, header, digest = candidate
    if source == 'pdb':
        choices = [{'url': f'https://files.rcsb.org/download/{identifier}.cif.gz', 'chain': chain,
                    'structure_url': f'https://www.rcsb.org/structure/{identifier}', 'title': identifier}]
    else:
        predictions = fs.fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{identifier}', cache) or []
        choices = [{'url': p['cifUrl'], 'chain': '', 'structure_url': f'https://alphafold.ebi.ac.uk/entry/{identifier}',
                    'title': p.get('gene', '') or p.get('uniprotDescription', 'AlphaFold prediction')}
                   for p in predictions if p.get('uniprotAccession') == identifier and p.get('cifUrl', '').startswith('https://alphafold.ebi.ac.uk/')]
    errors = []
    for choice in choices:
        data = fs.fetch(choice['url'], cache)
        if not data:
            continue
        suffix = '.cif.gz' if choice['url'].endswith('.gz') else '.cif'
        with tempfile.NamedTemporaryFile(suffix=suffix, dir=structures) as stream:
            stream.write(data)
            stream.flush()
            try:
                residues, author_chain, segment, plddt = crop_seed(stream.name, choice['chain'], sequence,
                                                                 minimum_plddt if source == 'alphafold' else None)
            except ValueError as error:
                errors.append(str(error))
                continue
            target = re.sub(r'[^A-Za-z0-9_]', '_', f'{model}_{source}_{identifier}_{author_chain}_{digest[:12]}')
            fs.write_ca(structures / (target + '.pdb'), residues)
            try:
                database, member = fs.api_model(model)
            except ValueError:
                database, member = 'cdd', model
            structures_url = (f'https://www.ebi.ac.uk/interpro/structure/pdb/entry/{database}/{member}/'
                              if model.startswith(('cd', 'pfam')) else None)
            return {'id': target, 'accession': member, 'cdd_accession': model, 'database': database,
                    'pssm_id': label[0], 'name': label[1], 'description': label[2],
                    'sequence': ''.join(fs.residue_letter(r) for r in residues),
                    'pdb': identifier if source == 'pdb' else None, 'chain': author_chain,
                    'uniprot': identifier if source == 'alphafold' else None,
                    'segments': [segment], 'structure_source': source, 'structure_url': choice['structure_url'],
                    'structure_title': choice['title'], 'structures_url': structures_url,
                    'mean_plddt': plddt if source == 'alphafold' else None,
                    'seed_header': header, 'seed_sequence_sha256': digest,
                    'evidence': 'Exact CDD seed-fragment match to ' + ('experimental PDB polymer sequence' if source == 'pdb' else 'AlphaFold polymer sequence'),
                    'entry_url': f'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={label[0]}'}
    raise ValueError('; '.join(errors) or 'No available prediction or coordinate file.')


def open_build_state(output):
    connection = sqlite3.connect(Path(output) / 'build.sqlite')
    connection.execute('CREATE TABLE IF NOT EXISTS families(model TEXT PRIMARY KEY, signature TEXT, status TEXT, records TEXT, errors TEXT)')
    return connection


def import_pilot(state, seed_library, output):
    seed_library = Path(seed_library)
    manifest = json.loads((seed_library / 'manifest.json').read_text())
    groups = defaultdict(list)
    for record in manifest['records']:
        identifier = record['id']
        if not re.fullmatch(r'[A-Za-z0-9_]+', identifier):
            raise ValueError('Invalid reference identifier in seed library.')
        shutil.copy2(seed_library / 'domains' / (identifier + '.pdb'), output / 'domains' / (identifier + '.pdb'))
        groups[record['cdd_accession']].append(record)
    for model, records in groups.items():
        state.execute('INSERT OR IGNORE INTO families VALUES(?,?,?,?,?)', (model, '', 'seeded', json.dumps(records), '[]'))
    state.commit()


def publish(output, state, inventory, catalog_metadata, executable):
    """Publish one immutable generation; the old manifest remains usable until ready."""
    rows = state.execute('SELECT model,status,records FROM families ORDER BY model').fetchall()
    records = [record for _, _, text in rows for record in json.loads(text)]
    if not records:
        raise ValueError('No verified references are available to publish.')
    generation = output / 'databases' / uuid4().hex
    generation.mkdir(parents=True)
    with tempfile.TemporaryDirectory(dir=output) as directory:
        inputs = Path(directory)
        for record in records:
            shutil.copy2(output / 'domains' / (record['id'] + '.pdb'), inputs / (record['id'] + '.pdb'))
        subprocess.run([executable, 'createdb', str(inputs), str(generation / 'foldseek_db'), '--threads', '2'],
                       check=True, stdout=subprocess.DEVNULL)
    with closing(sqlite3.connect(generation / 'references.sqlite')) as connection:
        connection.execute('CREATE TABLE refs(id TEXT PRIMARY KEY,record TEXT)')
        connection.executemany('INSERT INTO refs VALUES(?,?)',
                               [(r['id'], json.dumps(r)) for r in records])
        connection.commit()
    covered = len({r['cdd_accession'] for r in records})
    coverage = {**inventory, 'covered_models': covered, 'references': len(records),
                'attempted_models': len(rows), 'build_status': dict(Counter(status for _, status, _ in rows)),
                'reference_sources': dict(Counter(r.get('structure_source', 'pdb') for r in records))}
    manifest = {'schema': 2, 'pilot': True, 'description': 'CDD-wide reference build in progress; exact seed-fragment crops with documented provenance.',
                'database': str((generation / 'foldseek_db').relative_to(output)),
                'reference_metadata': str((generation / 'references.sqlite').relative_to(output)),
                'reference_count': len(records), 'coverage': coverage, 'catalog': catalog_metadata,
                'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                'foldseek_version': subprocess.check_output([executable, 'version'], text=True).strip()}
    atomic_json(output / 'manifest.json', manifest)
    atomic_json(output / 'coverage.json', coverage)
    return manifest


def expand(catalog_path, output, representatives=3, limit=None, seed_library=None,
           allow_alphafold=True, minimum_plddt=70, retry_failed=False, executable='foldseek', models=None):
    output = Path(output).resolve()
    if representatives < 1 or (limit is not None and limit < 1) or not 0 <= minimum_plddt <= 100:
        raise ValueError('Representative and batch limits must be positive; pLDDT must be between 0 and 100.')
    with build_lock(output), closing(sqlite3.connect(f'file:{Path(catalog_path).resolve()}?mode=ro', uri=True)) as inventory:
        structures, cache = output / 'domains', output / 'downloads'
        structures.mkdir(exist_ok=True)
        state = open_build_state(output)
        try:
            if seed_library:
                import_pilot(state, seed_library, output)
            metadata = {k: json.loads(v) for k, v in inventory.execute('SELECT key,value FROM metadata')}
            signature = hashlib.sha256(json.dumps([BUILD_VERSION, metadata, representatives, allow_alphafold, minimum_plddt], sort_keys=True).encode()).hexdigest()
            model_rows = inventory.execute('SELECT accession,pssm_id,name,description,source FROM models WHERE source != "cl" ORDER BY accession').fetchall()
            if models is not None:
                requested = set(models)
                known = {row[0] for row in model_rows}
                if requested - known:
                    raise ValueError('Unknown searchable CDD models: ' + ', '.join(sorted(requested - known)))
                model_rows = [row for row in model_rows if row[0] in requested]
            # Interleave sources so bounded batches grow coverage across the collection.
            groups = defaultdict(list)
            for row in model_rows:
                groups[row[4]].append(row)
            ordered = [rows[index] for index in range(max(map(len, groups.values()), default=0))
                       for _, rows in sorted(groups.items()) if index < len(rows)]
            seeded = {row[0] for row in state.execute('SELECT model FROM families WHERE status="seeded"')}
            ordered.sort(key=lambda row: row[0] not in seeded)
            attempted = 0
            try:
                for model, pssm, name, description, _ in ordered:
                    existing = state.execute('SELECT signature,status,records FROM families WHERE model=?', (model,)).fetchone()
                    if existing and existing[0] == signature and not (retry_failed and existing[1] == 'error'):
                        continue
                    candidates = inventory.execute('SELECT * FROM candidates WHERE model=? ORDER BY CASE source WHEN "pdb" THEN 0 ELSE 1 END,rowid', (model,)).fetchall()
                    if not allow_alphafold:
                        candidates = [row for row in candidates if row[1] == 'pdb']
                    records = json.loads(existing[2]) if existing else []
                    if any(r.get('structure_source', 'pdb') == 'pdb' for r in records):
                        candidates.sort(key=lambda row: row[1] != 'alphafold')
                    errors, transient_error = [], False
                    for candidate in candidates:
                        if len(records) >= representatives:
                            break
                        if any(r['sequence'] == candidate[4] for r in records):
                            continue
                        try:
                            record = candidate_record(candidate, (pssm, name, description), cache, structures, minimum_plddt)
                            if not any(r['sequence'] == record['sequence'] for r in records):
                                records.append(record)
                        except (OSError, URLError, ValueError, KeyError, IndexError, RuntimeError) as error:
                            transient_error |= isinstance(error, (OSError, URLError))
                            errors.append({'source': candidate[1], 'identifier': candidate[2], 'reason': str(error)})
                    status = 'complete' if len(records) >= representatives else 'error' if transient_error else 'partial' if records else 'unavailable' if candidates else 'no_candidates'
                    state.execute('INSERT OR REPLACE INTO families VALUES(?,?,?,?,?)', (model, signature, status, json.dumps(records), json.dumps(errors)))
                    state.commit()
                    if candidates or records:
                        attempted += 1
                        print(model, name, ':', len(records), 'references;', status, flush=True)
                    if limit is not None and attempted >= limit:
                        break
            except KeyboardInterrupt:
                print('Interrupted; publishing completed families. Rerun the same command to resume.', flush=True)
            return publish(output, state, catalog_summary(inventory), metadata, executable)
        finally:
            state.close()
