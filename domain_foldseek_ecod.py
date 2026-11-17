"""Install the complete published ECOD F40 Foldseek database with a coverage audit.

Uses bulk release files, without per-domain API requests or Gemmi. Publication
is atomic; the original archive and upstream checksums are retained.
"""

from collections import Counter
from contextlib import closing, contextmanager
import csv
import hashlib
import json
import mmap
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tarfile
from urllib.request import urlopen
from urllib.error import HTTPError
from uuid import uuid4

from domain_foldseek_catalog import atomic_json, build_lock

BASE = 'http://prodata.swmed.edu/ecod/distributions/'


@contextmanager
def staged_generation(output):
    path = Path(output) / 'databases' / uuid4().hex
    path.mkdir(parents=True)
    try:
        yield path
    except BaseException:
        shutil.rmtree(path)
        raise


def checksum(path, algorithm='md5'):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def download(url, path, expected=None):
    """Resume large transfers; checksum failures never become completed downloads."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and (expected is None or checksum(path) == expected):
        return path
    partial = path.with_name(path.name + '.part')
    if path.is_file():
        path.replace(partial)
    print('Downloading', path.name, flush=True)
    subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--continue-at', '-',
                    '--output', str(partial), url], check=True)
    if expected is not None and checksum(partial) != expected:
        raise ValueError('Checksum mismatch: ' + str(partial))
    partial.replace(path)
    return path


def checksums(text):
    values = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 2 and re.fullmatch(r'[a-f0-9]{32}', fields[0]):
            name = fields[1].lstrip('*')
            if Path(name).name != name:
                raise ValueError('Unexpected checksum filename.')
            values[name] = fields[0]
    if not values:
        raise ValueError('No checksums in the ECOD release manifest.')
    return values


def extract_database(archive, output, expected):
    """Extract only ordinary files explicitly named in the pinned release manifest."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    extracted = []
    with tarfile.open(archive, 'r|gz') as stream:
        for member in stream:
            if member.isdir():
                continue
            name = Path(member.name).name
            auxiliary = name in ('README.txt', 'missing_structures.tsv', 'excluded_uids.tsv') and member.size < 1_000_000
            if not member.isfile() or (name not in expected and not auxiliary) or name in extracted:
                raise ValueError('Unexpected ECOD archive member: ' + member.name)
            path = output / name
            with stream.extractfile(member) as source, path.open('wb') as target:
                shutil.copyfileobj(source, target)
            if name in expected and checksum(path) != expected[name]:
                raise ValueError('Extracted file checksum mismatch: ' + name)
            extracted.append(name)
    return extracted


def table_rows(path):
    with Path(path).open() as stream:
        yield from csv.DictReader((line for line in stream if not line.startswith('#')), delimiter='\t')


def pfam_mapping(path):
    with Path(path).open() as stream:
        return {row['f_id']: [p for p in row['pfam_acc'].split(',') if re.fullmatch(r'PF\d+', p)]
                for row in csv.DictReader(line for line in stream if not line.startswith('#'))}


def index_metadata(database, domains, families, destination):
    """Account for every actual Foldseek target, including chain-split targets."""
    mapping = pfam_mapping(families)
    counts, pfams, missing, indexed = Counter(), set(), [], 0
    with Path(str(database) + '.lookup').open() as stream:
        target_uids = {int(line.split('\t')[1].split('_')[0]) for line in stream}
    with closing(sqlite3.connect(destination)) as connection:
        connection.executescript('''
            CREATE TABLE domains(uid INTEGER PRIMARY KEY,record TEXT);
            CREATE TABLE refs(id TEXT PRIMARY KEY,record TEXT);
        ''')
        for row in table_rows(domains):
            if int(row['uid']) in target_uids:
                connection.execute('INSERT INTO domains VALUES(?,?)', (int(row['uid']), json.dumps(row)))
        connection.commit()
        offsets = {}
        with Path(str(database) + '.index').open() as stream:
            for line in stream:
                key, offset, length = line.split()
                offsets[int(key)] = (int(offset), int(length))
        with database.open('rb') as sequence_file, mmap.mmap(sequence_file.fileno(), 0, access=mmap.ACCESS_READ) as sequences, \
                Path(str(database) + '.lookup').open() as stream:
            for line in stream:
                key, identifier, _ = line.rstrip().split('\t')
                match = re.fullmatch(r'(\d+)(?:_(.+))?', identifier)
                if match is None:
                    raise ValueError('Unrecognized ECOD target identifier: ' + identifier)
                uid, chain_suffix = int(match[1]), match[2]
                row = connection.execute('SELECT record FROM domains WHERE uid=?', (uid,)).fetchone()
                domain = json.loads(row[0]) if row else {}
                offset, length = offsets[int(key)]
                sequence = sequences[offset:offset + length].rstrip(b'\n\0').decode('ascii')
                if not sequence or not sequence.isalpha():
                    raise ValueError('Invalid reference sequence: ' + identifier)
                family = domain.get('f_id', '')
                family_pfams = mapping.get(family, [])
                pfams.update(family_pfams)
                pdb = domain.get('pdb', '')
                source = 'pdb' if re.fullmatch(r'[A-Za-z0-9]{4}', pdb) else 'predicted_or_other'
                counts[source] += 1
                counts['split_chain_targets'] += bool(chain_suffix)
                entry_url = f'http://prodata.swmed.edu/ecod/domain/{uid}'
                name = ((domain.get('t_name') or domain.get('h_name')) if family.endswith('.0') else domain.get('f_name'))
                name = name or domain.get('t_name') or domain.get('h_name')
                name = name or (f'ECOD topology {family.rsplit(".", 1)[0]}' if domain else 'ECOD domain (metadata unavailable)')
                accession = 'ECOD:' + family if family else 'ECOD-UID:' + str(uid)
                description = '; '.join(f'{label}: {domain[field]}' for field, label in
                    [('architecture_name', 'Architecture'), ('x_name', 'Possible homology'),
                     ('h_name', 'Homology group'), ('t_name', 'Topology'), ('f_name', 'Family')]
                    if domain.get(field))
                if not domain:
                    missing.append(identifier)
                record = {'id': identifier, 'catalog': 'ECOD', 'catalog_accession': accession,
                          'accession': accession, 'database': 'ecod', 'name': name,
                          'description': description, 'length': len(sequence),
                          'ecod_uid': uid, 'ecod_domain': domain.get('ecod_domain_id'),
                          'pfam_accessions': family_pfams, 'entry_url': entry_url, 'structures_url': None,
                          'pdb': pdb if source == 'pdb' else None, 'chain': chain_suffix or domain.get('chain', ''),
                          'structure_source': source, 'structure_url': entry_url,
                          'structure_title': domain.get('ecod_domain_id') or identifier,
                          'ecod_seqid_range': domain.get('seqid_range'), 'ecod_pdb_range': domain.get('pdb_range'),
                          'evidence': 'Published ECOD F40 Foldseek representative, with upstream domain boundaries',
                          'metadata_status': 'mapped' if domain else 'unmapped'}
                connection.execute('INSERT INTO refs VALUES(?,?)', (identifier, json.dumps(record)))
                indexed += 1
                if indexed % 50000 == 0:
                    connection.commit()
                    print('Indexed', indexed, 'ECOD targets', flush=True)
        if indexed != len(offsets):
            raise ValueError('The sequence index and target lookup do not contain the same entries.')
        connection.execute('DROP TABLE domains')
        connection.commit()
        connection.execute('VACUUM')
    split_chains = counts.pop('split_chain_targets', 0)
    return {'targets': indexed, 'metadata_mapped': indexed - len(missing),
            'metadata_unmapped': len(missing), 'unmapped_target_ids': missing,
            'split_chain_targets': split_chains, 'structure_sources': dict(counts), 'pfam_accessions': sorted(pfams)}


def install(output, release='v295.2', cdd_catalog=None):
    """Install the complete upstream database; never silently sample its targets."""
    if not re.fullmatch(r'v\d+(?:\.\d+)*', release):
        raise ValueError('Expected an ECOD release such as v295.2.')
    output = Path(output).resolve()
    with build_lock(output), staged_generation(output) as generation:
        cache = output / 'downloads'
        cache.mkdir(exist_ok=True)
        manifest_file = cache / f'ecod.{release}.md5'
        if not manifest_file.is_file():
            data = None
            for name in (f'ecod.{release}.md5', 'ecod.latest.md5'):
                try:
                    with urlopen(BASE + name, timeout=45) as response:
                        candidate = response.read(1_000_000).decode('ascii')
                except HTTPError as error:
                    if error.code == 404:
                        continue
                    raise
                if all(item.startswith(f'ecod.{release}.') for item in checksums(candidate)):
                    data = candidate
                    break
            if data is None:
                raise ValueError('No checksums for the requested release; supply its pinned checksum file in downloads/.')
            manifest_file.write_text(data)
        expected = checksums(manifest_file.read_text())
        names = [f'ecod.{release}.F40.foldseek_db.tar.gz', f'ecod.{release}.F40.domains.txt',
                 f'ecod.{release}.f_id_pfam_acc.txt', f'ecod.{release}.F40.pdb.empty.tsv',
                 f'ecod.{release}.F40.pdb.missing.tsv']
        for name in names:
            download(BASE + name, cache / name, expected[name])
        extract_database(cache / names[0], generation, expected)
        database = generation / f'ecod.{release}.foldseek'
        coverage = index_metadata(database, cache / names[1], cache / names[2], generation / 'references.sqlite')
        if coverage['metadata_unmapped']:
            full = f'ecod.{release}.domains.txt'
            download(BASE + full, cache / full, expected[full])
            (generation / 'references.sqlite').unlink()
            coverage = index_metadata(database, cache / full, cache / names[2], generation / 'references.sqlite')
        coverage['release'] = release
        coverage['scope'] = 'All targets in the published ECOD F40 Foldseek database; not all CDD models.'
        coverage['upstream_empty_coordinate_records'] = sum(1 for _ in (cache / names[3]).open())
        coverage['upstream_missing_coordinate_records'] = sum(1 for _ in (cache / names[4]).open())
        if cdd_catalog:
            with closing(sqlite3.connect(f'file:{Path(cdd_catalog).resolve()}?mode=ro', uri=True)) as catalog:
                models = catalog.execute('SELECT accession,name FROM models WHERE source != "cl" ORDER BY accession').fetchall()
        else:
            from domain_foldseek import cdd_descriptions
            models = [(accession, label['name']) for accession, label in cdd_descriptions(cache).items()
                      if not accession.startswith('cl')]
        pfams, counts = set(coverage['pfam_accessions']), Counter()
        with (generation / 'cdd_coverage.tsv').open('w') as report:
            report.write('cdd_accession\tname\tstatus\n')
            for accession, name in models:
                mapped = accession.startswith('pfam') and 'PF' + accession[4:] in pfams
                status = 'Pfam family cross-reference present' if mapped else 'no direct CDD mapping'
                counts[status] += 1
                report.write(f'{accession}\t{name}\t{status}\n')
        coverage['cdd_model_audit'] = dict(counts)
        manifest = {'schema': 2, 'pilot': False, 'catalog': 'ECOD', 'release': release,
                    'description': 'Complete published ECOD F40 Foldseek database; CDD cross-references audited separately.',
                    'database': str(database.relative_to(output)),
                    'reference_metadata': str((generation / 'references.sqlite').relative_to(output)),
                    'reference_count': coverage['targets'], 'coverage': coverage,
                    'upstream_checksums_sha256': checksum(manifest_file, 'sha256')}
        atomic_json(output / 'coverage.json', coverage)
        atomic_json(output / 'manifest.json', manifest)
        return manifest


def merge(libraries, output, executable='foldseek', replace=False):
    """Combine complete libraries and their metadata without rebuilding coordinates."""
    import domain_foldseek as fs

    libraries, output = [Path(p).resolve() for p in libraries], Path(output).resolve()
    if len(libraries) != 2:
        raise ValueError('Supply two libraries; additional merges can be chained.')
    if output in libraries and not replace:
        raise ValueError('Use a separate output directory for the merged library.')
    manifest_texts = [(p / 'manifest.json').read_bytes() for p in libraries]
    manifests = [json.loads(text) for text in manifest_texts]
    with build_lock(output), staged_generation(output) as generation:
        for library, text in zip(libraries, manifest_texts):
            if (library / 'manifest.json').read_bytes() != text:
                raise ValueError('An input library changed while the integration was starting.')
        metadata = generation / 'references.sqlite'
        with closing(sqlite3.connect(metadata)) as connection:
            connection.execute('CREATE TABLE refs(id TEXT PRIMARY KEY,record TEXT)')
            for library, manifest in zip(libraries, manifests):
                if manifest['schema'] == 2:
                    path = fs.library_file(library, manifest['reference_metadata'])
                    with closing(sqlite3.connect(f'file:{path}?mode=ro', uri=True)) as source:
                        connection.executemany('INSERT INTO refs VALUES(?,?)', source.execute('SELECT id,record FROM refs'))
                else:
                    connection.executemany('INSERT INTO refs VALUES(?,?)',
                        [(r['id'], json.dumps(r)) for r in manifest['records']])
            connection.commit()
            count = connection.execute('SELECT count(*) FROM refs').fetchone()[0]
        databases = [fs.library_file(library, manifest['database']) for library, manifest in zip(libraries, manifests)]
        for suffix in ('', '_h', '_ss', '_ca'):
            subprocess.run([executable, 'concatdbs', str(databases[0]) + suffix,
                            str(databases[1]) + suffix, str(generation / 'foldseek_db') + suffix,
                            '--threads', '1', '-v', '1'], check=True)
        # Copy indexed profiles into this immutable generation, rather than
        # leaving an external path or a hard link to mutable build checkpoints.
        profile_manifest = None
        for library, source in zip(libraries, manifests):
            if not source.get('profile_database'):
                continue
            if profile_manifest is not None:
                raise ValueError('Merge one profile-bearing release at a time.')
            source_profile = fs.library_file(library, source['profile_database'])
            profile_folder = generation / 'profiles' / source_profile.parent.name
            profile_folder.mkdir(parents=True)
            for suffix in ('', '.ssi', '.h3f', '.h3i', '.h3m', '.h3p'):
                path = Path(str(source_profile) + suffix)
                shutil.copyfile(path, profile_folder / path.name)
            for name in ('Pfam.version.gz', 'Pfam-A.hmm.dat.gz', 'profile-manifest.json'):
                shutil.copyfile(source_profile.parent / name, profile_folder / name)
            profile_manifest = {
                'profile_database': str((profile_folder / source_profile.name).relative_to(output)),
                'pfam_coverage': source.get('coverage', {}).get('pfam', source.get('coverage', {})),
                'pfam_release': source.get('pfam_release', source.get('release')),
            }
        manifest = {'schema': 2, 'catalog': 'merged', 'pilot': any(m.get('pilot') for m in manifests),
                    'database': str((generation / 'foldseek_db').relative_to(output)),
                    'reference_metadata': str(metadata.relative_to(output)), 'reference_count': count,
                    'sources': [{'path': str(p), 'manifest': m} for p, m in zip(libraries, manifests)]}
        coverage = {'scope': 'Union of the two supplied libraries, preserving reference provenance.',
                    'targets': count, 'sources': [{'library': str(p), 'targets': fs.library_reference_count(m)}
                                                for p, m in zip(libraries, manifests)]}
        if profile_manifest:
            manifest['profile_database'] = profile_manifest['profile_database']
            manifest['pfam_release'] = profile_manifest['pfam_release']
            coverage['pfam'] = profile_manifest['pfam_coverage']
        manifest['coverage'] = coverage
        atomic_json(generation / 'coverage.json', coverage)
        if (output / 'manifest.json').is_file():
            previous = (output / 'manifest.json').read_bytes()
            history = output / 'history'
            history.mkdir(exist_ok=True)
            atomic_json(history / (hashlib.sha256(previous).hexdigest() + '.json'), json.loads(previous))
        atomic_json(output / 'coverage.json', coverage)
        atomic_json(output / 'manifest.json', manifest)
        return manifest
