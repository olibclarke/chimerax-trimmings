"""CCD density screening and its native ChimeraX browser.

The GUI owns models. Background jobs exchange arrays, component records and
small messages only; no map or AtomicStructure is passed to a fitting worker.
"""

from pathlib import Path
from contextlib import closing
import json
import math
import os
import queue
import sqlite3
import threading
import time

CCD_URL = 'https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz'
PDB_IDS_URL = 'https://data.rcsb.org/rest/v1/holdings/current/entry_ids'
PDB_DATA_URL = 'https://data.rcsb.org/graphql'
CHEBI_URL = 'https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_core.obo.gz'
CHEMONT_SERVICE = 'https://classyfire2.gnps2.org'
GROUP_SCHEMA = '4'
# Stable command keys, GUI labels, ChEBI roots, and chemical-class/role provenance.
CHEMICAL_GROUPS = (
    ('nucleotides', 'Nucleotides / derivatives', ('CHEBI:36976', 'CHEBI:231540'), 'class'),
    ('nucleosides', 'Nucleosides', ('CHEBI:33838',), 'class'),
    ('nucleobases', 'Nucleobases', ('CHEBI:18282',), 'class'),
    ('amino_acids', 'Amino acids / derivatives', ('CHEBI:33709', 'CHEBI:83821'), 'class'),
    ('peptides', 'Peptides', ('CHEBI:16670',), 'class'),
    ('lipids', 'Lipids', ('CHEBI:18059',), 'class'),
    ('fatty_acids', 'Fatty acids / anions', ('CHEBI:35366', 'CHEBI:28868'), 'class'),
    ('phospholipids', 'Phospholipids', ('CHEBI:16247',), 'class'),
    ('glycerolipids', 'Glycerolipids', ('CHEBI:35741',), 'class'),
    ('sphingolipids', 'Sphingolipids', ('CHEBI:26739',), 'class'),
    ('steroids', 'Steroids', ('CHEBI:35341',), 'class'),
    ('sterols', 'Sterols', ('CHEBI:15889',), 'class'),
    ('isoprenoids', 'Isoprenoids / terpenoids', ('CHEBI:24913', 'CHEBI:26873'), 'class'),
    ('carbohydrates', 'Carbohydrates / derivatives', ('CHEBI:16646', 'CHEBI:63299', 'CHEBI:24400'), 'class'),
    ('glycosides', 'Glycosides', ('CHEBI:24400',), 'class'),
    ('porphyrins', 'Porphyrins', ('CHEBI:8337',), 'class'),
    ('metabolites', 'Metabolites (role)', ('CHEBI:25212',), 'role'),
    ('drugs', 'Drugs (role)', ('CHEBI:23888',), 'role'),
    ('inhibitors', 'Inhibitors (role)', ('CHEBI:35222',), 'role'),
    ('agonists', 'Agonists (role)', ('CHEBI:48705',), 'role'),
    ('antagonists', 'Antagonists (role)', ('CHEBI:48706',), 'role'),
    ('toxins', 'Toxins (role)', ('CHEBI:27026',), 'role'),
    ('cofactors', 'Cofactors (role)', ('CHEBI:23357',), 'role'),
    ('vitamins', 'Vitamins (role)', ('CHEBI:33229',), 'role'),
    ('detergents', 'Detergents / surfactants (role)', ('CHEBI:27780', 'CHEBI:35195'), 'role'),
    ('buffers', 'Buffers (role)', ('CHEBI:35225',), 'role'),
    ('benzenoids', 'Benzenoids', (), 'class'),
    ('heterocycles', 'Organoheterocyclic compounds', (), 'class'),
    ('organic_acids', 'Organic acids / derivatives', (), 'class'),
    ('organic_oxygen', 'Organic oxygen compounds', (), 'class'),
    ('organic_nitrogen', 'Organic nitrogen compounds', (), 'class'),
    ('organohalogens', 'Organohalogen compounds', (), 'class'),
    ('organosulfur', 'Organosulfur compounds', (), 'class'),
    ('organophosphorus', 'Organophosphorus compounds', (), 'class'),
    ('hydrocarbons', 'Hydrocarbons', (), 'class'),
    ('alkaloids', 'Alkaloids / derivatives', (), 'class'),
    ('phenylpropanoids', 'Phenylpropanoids / polyketides', (), 'class'),
    ('organometallics', 'Organometallic compounds', (), 'class'),
    ('inorganic', 'Inorganic compounds', (), 'class'),
)
# Exact ChemOnt categories, not compound-name keyword guesses. Detailed classes
# remain independently selectable even when no convenient preset covers them.
CHEMONT_PRESETS = {
    'nucleotides': ('0001506', '0001509', '0001297', '0001329', '0003394', '0003468', '0004351'),
    'nucleosides': ('0000479', '0000480', '0004174', '0001134', '0003738'),
    'amino_acids': ('0000347',), 'peptides': ('0000348',),
    'lipids': ('0000012',), 'fatty_acids': ('0000262',),
    'phospholipids': ('0000256', '0003392'), 'glycerolipids': ('0000175',),
    'sphingolipids': ('0000257',), 'steroids': ('0000258',),
    'sterols': ('0001469', '0001403', '0003061'), 'carbohydrates': ('0000011',),
    'isoprenoids': ('0000259',), 'glycosides': ('0002105',),
    'porphyrins': ('0000212', '0000452'), 'benzenoids': ('0002448',),
    'heterocycles': ('0000002',), 'organic_acids': ('0000264',),
    'organic_oxygen': ('0004603',), 'organic_nitrogen': ('0000278',),
    'organohalogens': ('0000267',), 'organosulfur': ('0000004',),
    'organophosphorus': ('0000400',), 'hydrocarbons': ('0002837',),
    'alkaloids': ('0000279',), 'phenylpropanoids': ('0000261',),
    'organometallics': ('0000462',), 'inorganic': ('0000001',),
}
SCHEMA = '1'
VOLUME_LEVELS = (.15, .30, .50)
VOLUME_RESOLUTION = 2.6
STANDARD_EXTRA_STARTS, STANDARD_SAMPLES = 24, 3000
FAST_EXTRA_STARTS, FAST_SAMPLES = 8, 1500
DEFAULT_CONFORMERS, MAX_CONFORMERS = 8, 32
FINAL_CONFORMERS, MAX_FINAL_CONFORMERS = 64, 256
DEFAULT_LIBRARY_SOURCE = 'ccd_chebi'
CONFORMER_RMSD = .25  # Ignore torsion variants indistinguishable after rigid alignment.
BLOB_HIGHLIGHT_RGB = (80, 230, 120)
INTERNAL_OVERLAP = 1.6  # Reject newly introduced gross heavy-atom overlaps, in Angstroms.
SHAPE = None
CHEBI = None
IMPLEMENTATION = object()
# Pin the download to an immutable data snapshot, not the latest software release.
LIBRARY_RELEASE = {
    'url': ('https://github.com/olibclarke/chimerax-trimmings/releases/download/'
            'ligands-chebi257/trimmings-ligands-chebi257.tar.gz'),
    'bytes': 61_249_832,
    'sha256': '4763005e51354ac15efaddd01355779d96f65ba59ada144a7e5b77fd47e4dcc9',
}
LIBRARY_FILES = ('ccd.sqlite', 'chebi.sqlite', 'groups.sqlite', 'frequencies.json')


def shape_api():
    global SHAPE
    if SHAPE is None:
        from types import SimpleNamespace
        import runpy

        SHAPE = SimpleNamespace(**runpy.run_path(str(Path(__file__).with_name('chimerax_trimmings.py'))))
    return SHAPE


def default_cache():
    from chimerax.core.fetch import cache_directories

    return Path(cache_directories()[0]) / 'LigandSearch'


def chebi_api():
    """Load a sibling helper even when runscript was called from another folder."""
    import importlib.util
    import sys

    global CHEBI
    path = Path(__file__).with_name('ligand_chebi.py')
    if not path.is_file():
        raise ValueError('Download ligand_chebi.py beside ligand_search.py to use the ligand browser.')
    stamp = path.stat().st_mtime_ns
    if CHEBI is None or CHEBI._source_stamp != stamp:
        spec = importlib.util.spec_from_file_location('_trimmings_ligand_chebi', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ls = sys.modules[__name__]
        module._source_stamp = stamp
        CHEBI = module
    return CHEBI


def component_record(path):
    """Read complete heavy-atom conformers; never mix partial coordinate sets."""
    import numpy as np
    from chimerax.atomic import Element
    from chimerax.geometry import align_points
    from chimerax.mmcif import get_cif_tables

    comp, atoms, bonds = get_cif_tables(str(path), ['chem_comp', 'chem_comp_atom', 'chem_comp_bond'])
    if comp is None or atoms is None:
        return None
    ident, name, formula, kind, status = comp.fields(
        ['id', 'name', 'formula', 'type', 'pdbx_release_status'], allow_missing_fields=True)[0]
    if status in ('OBS', 'HOLD', 'DEL'):
        return None
    rows = atoms.fields(['atom_id', 'type_symbol', 'pdbx_model_Cartn_x_ideal',
                         'pdbx_model_Cartn_y_ideal', 'pdbx_model_Cartn_z_ideal',
                         'model_Cartn_x', 'model_Cartn_y', 'model_Cartn_z'], allow_missing_fields=True)
    rows = [r for r in rows if r[1] not in ('H', 'D')]
    if len(rows) < 4:
        return None  # Rigid shape orientation is undefined for ions/tiny fragments.
    names = [r[0] for r in rows]
    elements = [Element.get_element(r[1]).number for r in rows]
    if min(elements) <= 0 or len(set(names)) != len(names):
        return None
    conformers = []
    for source, start in (('ideal', 2), ('observed', 5)):
        try:
            xyz = np.array([r[start:start+3] for r in rows], dtype=np.float64)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(xyz).all() or np.linalg.matrix_rank(xyz - xyz.mean(axis=0)) < 2:
            continue
        if any(align_points(xyz, np.array(c['xyz']))[1] < .15 for c in conformers):
            continue
        conformers.append(dict(source=source, xyz=xyz.tolist()))
    if not conformers:
        return None
    index = {name: i for i, name in enumerate(names)}
    connections = []
    if bonds is not None:
        for a, b, order in bonds.fields(['atom_id_1', 'atom_id_2', 'value_order'], allow_missing_fields=True):
            if a in index and b in index and a != b:
                connections.append((index[a], index[b], order))
    record = dict(id=ident, name=name, formula=formula, type=kind, names=names,
                  elements=elements, bonds=connections, conformers=conformers)
    return describe_component(record)


def describe_component(record):
    """The same density-volume prior for CCD and externally sourced conformers."""
    import numpy as np

    # A small, reusable descriptor permits immediate queueing of a full library.
    # It is a broad scheduling prior, not a molecular identity or fit score.
    xyz = np.array(record['conformers'][0]['xyz'], dtype=np.float32)
    xyz -= xyz.mean(axis=0)
    try:
        grid = shape_api()._shape_grid(xyz, np.array(record['elements'], np.float32), VOLUME_RESOLUTION)
        values = grid.matrix()
        peak, voxel = float(values.max()), float(np.prod(grid.step))
        volumes = [float(np.count_nonzero(values >= peak * level) * voxel) for level in VOLUME_LEVELS]
    except Exception:
        # Large CCD components remain eligible, but use a conservative ellipsoid
        # descriptor if their full numerical density would exceed a grid limit.
        axes = np.sqrt(np.maximum(np.linalg.eigvalsh(np.cov(xyz.T)), .1))
        volume = float(4 * np.pi / 3 * np.prod(np.sqrt(5) * axes + 1.5))
        volumes = [volume * factor for factor in (1.3, 1., .7)]
    record['volumes'] = volumes
    return record


def cif_blocks(stream):
    """Stream CCD blocks without treating data_ inside semicolon text as a header."""
    lines, text = [], False
    for line in stream:
        if not text and line.lower().startswith('data_'):
            if lines:
                yield ''.join(lines)
            lines = [line]
        elif lines:
            lines.append(line)
        if line.startswith(';'):
            text = not text
    if lines:
        yield ''.join(lines)


def _cancelled(cancel):
    if cancel.is_set():
        raise InterruptedError('Stopped by user.')


def download_ccd(cache, cancel, report, refresh=False):
    """Resume an interrupted public archive download, checking its version."""
    import urllib.request

    archive = cache / 'components.cif.gz'
    if archive.is_file() and not refresh:
        return archive
    part, info = cache / 'components.cif.gz.part', cache / 'download.json'
    metadata = json.loads(info.read_text()) if info.is_file() else {}
    start = part.stat().st_size if part.exists() and metadata.get('etag') else 0
    headers = {'User-Agent': 'ChimeraX-ligandsearch/1'}
    if start:
        headers.update(Range=f'bytes={start}-', **{'If-Range': metadata['etag']})
    if archive.is_file() and not start and metadata.get('etag'):
        headers['If-None-Match'] = metadata['etag']
    try:
        response = urllib.request.urlopen(urllib.request.Request(CCD_URL, headers=headers), timeout=30)
    except urllib.error.HTTPError as error:
        if error.code == 304 and archive.is_file():
            return archive
        raise
    with response:
        append = start > 0 and response.status == 206
        if append and not response.headers.get('Content-Range', '').startswith(f'bytes {start}-'):
            raise RuntimeError('Unexpected CCD download range; remove the partial archive and retry.')
        if not append:
            start = 0
        info.write_text(json.dumps(dict(etag=response.headers.get('ETag'))))
        total = start + int(response.headers.get('Content-Length', 0))
        done, next_report = start, 0.
        with part.open('ab' if append else 'wb') as output:
            while True:
                _cancelled(cancel)
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
                done += len(chunk)
                if time.monotonic() >= next_report:
                    report(f'Downloading CCD: {done / 2**20:.0f}' +
                           (f'/{total / 2**20:.0f}' if total else '') + ' MB')
                    next_report = time.monotonic() + .25
        if total and done != total:
            raise RuntimeError('CCD download was incomplete; retry to resume it.')
    os.replace(part, archive)
    return archive


def build_library(cache, source=None, cancel=None, report=lambda message: None, refresh=False):
    """Atomically replace an index only after the entire source was parsed."""
    import gzip
    import tempfile
    import uuid
    import zlib

    cache = Path(cache).expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)
    cancel = cancel or threading.Event()
    lock = cache / 'build.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        try:
            pid = int(lock.read_text())
            os.kill(pid, 0)
        except (ValueError, ProcessLookupError):
            lock.unlink(missing_ok=True)
            return build_library(cache, source, cancel, report, refresh)
        raise RuntimeError('Another CCD library build is running; wait for it or stop that build.') from None
    with os.fdopen(fd, 'w') as f:
        f.write(str(os.getpid()))
    temporary = cache / ('ccd-' + uuid.uuid4().hex + '.building')
    connection = None
    managed_source = source is None
    try:
        source = Path(source).expanduser() if source is not None else download_ccd(cache, cancel, report, refresh)
        connection = sqlite3.connect(temporary)
        connection.executescript('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);'
                                 'CREATE TABLE components(id TEXT PRIMARY KEY,name TEXT,atoms INTEGER,'
                                 'v0 REAL,v1 REAL,v2 REAL,payload TEXT);')
        processed, eligible, rejected = 0, 0, 0
        with tempfile.TemporaryDirectory(prefix='ccd_parse_') as work:
            component_path = Path(work) / 'component.cif'
            def add(path):
                nonlocal processed, eligible, rejected
                _cancelled(cancel)
                processed += 1
                try:
                    record = component_record(path)
                except Exception:
                    record = None
                    rejected += 1
                if record is not None:
                    connection.execute('INSERT OR REPLACE INTO components VALUES(?,?,?,?,?,?,?)',
                                       (record['id'], record['name'], len(record['names']),
                                        *record['volumes'], json.dumps(record, separators=(',', ':'))))
                    eligible += 1
                if processed % 200 == 0:
                    connection.commit()
                    report(f'Indexing CCD: {processed:,} entries checked; {eligible:,} usable')
            if source.is_dir():
                for path in sorted(source.glob('*.cif')):
                    add(path)
            else:
                opener = gzip.open if source.suffix == '.gz' else open
                with opener(source, 'rt', encoding='utf-8') as stream:
                    for block in cif_blocks(stream):
                        component_path.write_text(block)
                        add(component_path)
        _cancelled(cancel)
        count = connection.execute('SELECT count(*) FROM components').fetchone()[0]
        if count == 0:
            raise RuntimeError('No usable CCD components were found; the previous library was kept.')
        metadata = dict(schema=SCHEMA, count=str(count), checked=str(processed), parse_errors=str(rejected),
                        created=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), source=str(source))
        connection.executemany('INSERT INTO metadata VALUES(?,?)', metadata.items())
        connection.commit()
        connection.close()
        connection = None
        os.replace(temporary, cache / 'ccd.sqlite')
        report(f'CCD library ready: {count:,} usable components')
        return metadata
    except (gzip.BadGzipFile, EOFError, zlib.error) as error:
        if managed_source:
            (cache / 'components.cif.gz').unlink(missing_ok=True)
            raise RuntimeError('The CCD download is damaged; retry Expand/Update CCD library to download it again. '
                               'The previous index was kept.') from error
        raise RuntimeError('The supplied CCD archive is damaged; use a complete archive and retry. '
                           'The previous index was kept.') from error
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)
        lock.unlink(missing_ok=True)


def library_info(cache):
    path = Path(cache).expanduser().resolve() / 'ccd.sqlite'
    if not path.is_file():
        return None
    with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as db:
        metadata = dict(db.execute('SELECT key,value FROM metadata'))
    if metadata.get('schema') != SCHEMA:
        raise ValueError('The CCD index needs rebuilding; click Expand/Update CCD library.')
    return metadata


def _file_digest(path, cancel):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk := stream.read(1024*1024):
            _cancelled(cancel)
            digest.update(chunk)
    return digest.hexdigest()


def _library_not_building(cache):
    for name in ('build.lock', 'chebi_build.lock'):
        path = Path(cache) / name
        if path.is_file():
            try:
                os.kill(int(path.read_text()), 0)
            except (ValueError, ProcessLookupError):
                continue
            raise ValueError('The ligand library is still being built. Wait for completion before packaging or installing.')


def _validate_library_snapshot(cache):
    """Validate a staged copy, never modify the source download/cache."""
    cache = Path(cache)
    if any(not (cache / name).is_file() for name in LIBRARY_FILES):
        raise ValueError('The combined library is incomplete: CCD, ChEBI, chemical groups and PDB counts are required.')
    for name in ('ccd.sqlite', 'chebi.sqlite', 'groups.sqlite'):
        with closing(sqlite3.connect((cache / name).as_uri() + '?mode=ro', uri=True)) as db:
            if db.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
                raise ValueError(f'The library database {name} failed its integrity check.')
    ccd = library_info(cache)
    chebi = chebi_api().info(cache)
    groups = group_info(cache)
    if ccd is None or chebi is None or not int(chebi['count']):
        raise ValueError('The combined library needs usable CCD and ChEBI entries.')
    if chebi.get('state') != 'complete' or int(chebi.get('errors', 0)):
        raise ValueError('The ChEBI expansion is not complete. Finish or resume its update before packaging.')
    if groups is None or not groups['current'] or frequency_info(cache) is None:
        raise ValueError('The library needs current chemical-group tags and valid PDB occurrence counts.')
    with closing(sqlite3.connect((cache / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        ids = {row[0] for row in db.execute('SELECT id FROM components')}
    with closing(sqlite3.connect((cache / 'chebi.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        extra_ids = [row[0] for row in db.execute('SELECT id FROM components')]
    if any(not ident.startswith('CHEBI:') and ident not in ids for ident in extra_ids):
        raise ValueError('ChEBI references an older CCD index; update the expansion before packaging.')
    new = sum(ident not in ids for ident in extra_ids)
    return dict(ccd=int(ccd['count']), chebi=int(chebi['count']), chebi_new=new,
                combined=len(ids)+new, ontology=chebi.get('ontology', 'unknown'))


def package_library(cache, output, asset_url, cancel=None, report=lambda message: None):
    """Create a distributable release from consistent SQLite snapshots only."""
    import gzip
    import tarfile
    import tempfile
    from urllib.parse import urlparse

    cache, output = Path(cache).expanduser().resolve(), Path(output).expanduser().resolve()
    cancel = cancel or threading.Event()
    if urlparse(asset_url).scheme != 'https':
        raise ValueError('asset_url must be the HTTPS URL of the future GitHub release archive.')
    _library_not_building(cache)
    state = chebi_api().info(cache)
    if state is None or state.get('state') != 'complete' or int(state.get('errors', 0)):
        raise ValueError('The ChEBI expansion is not complete. Wait for it to finish before packaging.')
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise ValueError('The output archive already exists; choose a new filename.')
    descriptor = output.with_name(output.name + '.release.json')
    if descriptor.exists():
        raise ValueError('The release descriptor already exists; choose a new archive filename.')
    temporary = output.with_name(output.name + '.part')
    lock = output.with_name(output.name + '.lock')
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise ValueError('A package operation or stale package lock exists; use another output filename.') from None
    with os.fdopen(fd, 'w') as stream:
        stream.write(str(os.getpid()))
    try:
        if output.exists() or descriptor.exists():
            raise ValueError('A package with this filename was created meanwhile; choose another output filename.')
        temporary.unlink(missing_ok=True)
        with tempfile.TemporaryDirectory(prefix='ligand_package_') as directory:
            staged = Path(directory)
            for name in LIBRARY_FILES:
                _cancelled(cancel)
                source = cache / name
                if not source.is_file():
                    raise ValueError(f'The library is missing {name}; complete its update first.')
                report('Ligand library: snapshotting ' + name)
                if name.endswith('.sqlite'):
                    with closing(sqlite3.connect(source.as_uri() + '?mode=ro', uri=True)) as first, \
                            closing(sqlite3.connect(staged / name)) as second:
                        first.backup(second, pages=256, progress=lambda *_: _cancelled(cancel))
                        if name == 'ccd.sqlite':
                            second.execute("UPDATE metadata SET value=? WHERE key='source'", (CCD_URL,))
                            second.commit()
                            # Do not retain the old local path in unused SQLite pages.
                            second.execute('VACUUM')
                else:
                    (staged / name).write_bytes(source.read_bytes())
            _library_not_building(cache)
            counts = _validate_library_snapshot(staged)
            notice = Path(__file__).with_name('DATA_NOTICE.md')
            if not notice.is_file():
                raise ValueError('Keep DATA_NOTICE.md beside the script when packaging library data.')
            (staged / 'DATA_NOTICE.md').write_bytes(notice.read_bytes())
            names = (*LIBRARY_FILES, 'DATA_NOTICE.md')
            manifest = dict(schema=1, created=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), counts=counts,
                            files={name: dict(bytes=(staged / name).stat().st_size,
                                              sha256=_file_digest(staged / name, cancel)) for name in names},
                            sources=dict(ccd=CCD_URL, chebi=CHEBI_URL, pubchem='https://pubchem.ncbi.nlm.nih.gov/docs/pubchem3d',
                                         pdb=PDB_DATA_URL, chemical_groups=CHEMONT_SERVICE))
            (staged / 'ligand-library.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
            report('Ligand library: compressing completed library')
            with temporary.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as compressed, \
                    tarfile.open(fileobj=compressed, mode='w') as archive:
                for name in ('ligand-library.json', *names):
                    _cancelled(cancel)
                    info = archive.gettarinfo(str(staged / name), arcname=name)
                    info.uid = info.gid = info.mtime = 0
                    info.uname = info.gname = ''
                    info.mode = 0o644
                    with (staged / name).open('rb') as stream:
                        archive.addfile(info, stream)
            release = dict(url=asset_url, bytes=temporary.stat().st_size, sha256=_file_digest(temporary, cancel))
            os.replace(temporary, output)
            descriptor.write_text(json.dumps(release, indent=2), encoding='utf-8')
            report(f'Ligand library: packaged {output.stat().st_size / 1e6:.1f} MB; release descriptor {descriptor.name}')
            return release
    finally:
        temporary.unlink(missing_ok=True)
        lock.unlink(missing_ok=True)


def install_library_archive(cache, archive_path, cancel=None, report=lambda message: None):
    """Validate and atomically install a flat bundle into a fresh cache only."""
    import tarfile
    import tempfile

    cache, archive_path = Path(cache).expanduser().resolve(), Path(archive_path)
    cancel = cancel or threading.Event()
    _library_not_building(cache)
    if cache.exists() and any(cache.iterdir()):
        raise ValueError('A local ligand cache already exists; it was not overwritten. Choose an empty cache folder.')
    cache.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.ligand_install_', dir=cache.parent) as directory:
        prepared = Path(directory) / 'library'
        prepared.mkdir()
        with tarfile.open(archive_path, 'r:gz') as archive:
            first = archive.next()
            if first is None or first.name != 'ligand-library.json' or not first.isfile() or first.size > 1024*1024:
                raise ValueError('The library archive has no valid manifest.')
            with archive.extractfile(first) as stream:
                manifest = json.load(stream)
            import re
            files = manifest.get('files', {}) if isinstance(manifest, dict) else {}
            if (not isinstance(manifest, dict) or manifest.get('schema') != 1 or not isinstance(files, dict) or
                    set(files) != {*LIBRARY_FILES, 'DATA_NOTICE.md'}):
                raise ValueError('The library archive has an unsupported file manifest.')
            if any(not isinstance(item, dict) or type(item.get('bytes')) is not int or not 0 < item['bytes'] <= 2_000_000_000 or
                   not isinstance(item.get('sha256'), str) or re.fullmatch(r'[a-f0-9]{64}', item['sha256']) is None
                   for item in files.values()):
                raise ValueError('The library manifest has invalid sizes or checksums.')
            if sum(item['bytes'] for item in files.values()) > 2_000_000_000:
                raise ValueError('The library archive exceeds the supported installed size.')
            seen, manifest_seen = set(), False
            for member in archive:
                if member.name == 'ligand-library.json':
                    if manifest_seen:
                        raise ValueError('The library archive contains duplicate manifests.')
                    manifest_seen = True
                    continue  # TarFile iteration includes the first member again.
                _cancelled(cancel)
                if not member.isfile() or member.name not in files or member.name in seen:
                    raise ValueError('The library archive contains unexpected files or links.')
                meta = files[member.name]
                if member.size != meta['bytes'] or not 0 < member.size <= 2_000_000_000:
                    raise ValueError('The library archive has an invalid file size.')
                target = prepared / member.name
                with archive.extractfile(member) as stream, target.open('xb') as output:
                    while chunk := stream.read(1024*1024):
                        _cancelled(cancel)
                        output.write(chunk)
                if _file_digest(target, cancel) != meta['sha256']:
                    raise ValueError(f'Library integrity check failed for {member.name}; download again.')
                seen.add(member.name)
                report('Ligand library: verified ' + member.name)
            if seen != set(files):
                raise ValueError('The library archive is incomplete; download it again.')
        _validate_library_snapshot(prepared)
        (prepared / 'ligand-library.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        _cancelled(cancel)
        # Recheck after extraction; never replace a cache created meanwhile.
        if cache.exists() and any(cache.iterdir()):
            raise ValueError('A local cache was created during installation; it was kept unchanged.')
        os.replace(prepared, cache)
    report('Ligand library: combined CCD + ChEBI library ready')


def install_prebuilt_library(cache, release=None, cancel=None, report=lambda message: None):
    """Serialize first-use downloads across ChimeraX instances sharing a cache."""
    release = LIBRARY_RELEASE if release is None else release
    if release is None:
        return False
    cache = Path(cache).expanduser().resolve()
    if cache.exists() and any(cache.iterdir()):
        return False
    cache.parent.mkdir(parents=True, exist_ok=True)
    # OS-managed locks are released on crash; leave the small lock file in place
    # so another process cannot lock a different inode during a rename/unlink.
    with (cache.parent / '.ligand-library-download.lock').open('a+b') as lock:
        if os.name == 'nt':
            import msvcrt
            lock.seek(0)
            if not lock.read(1):
                lock.write(b'0')
                lock.flush()
            lock.seek(0)
            try:
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError:
                raise ValueError('Another ligand-library download is running; wait for it to finish.') from None
        else:
            import fcntl
            try:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ValueError('Another ligand-library download is running; wait for it to finish.') from None
        try:
            return _install_prebuilt_library(cache, release, cancel, report)
        finally:
            if os.name == 'nt':
                lock.seek(0)
                msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _install_prebuilt_library(cache, release=None, cancel=None, report=lambda message: None):
    """Download the pinned first-use archive, reusing an interrupted transfer."""
    import urllib.request

    release = LIBRARY_RELEASE if release is None else release
    if release is None:
        return False
    cache = Path(cache).expanduser().resolve()
    cancel = cancel or threading.Event()
    if cache.exists() and any(cache.iterdir()):
        return False
    from urllib.parse import urlparse
    import re
    if (urlparse(release['url']).scheme != 'https' or not 0 < release['bytes'] <= 2_000_000_000 or
            re.fullmatch(r'[a-f0-9]{64}', release['sha256']) is None):
        raise ValueError('The published ligand-library descriptor is invalid; update Trimmings.')
    cache.parent.mkdir(parents=True, exist_ok=True)
    part = cache.parent / ('ligand-' + release['sha256'] + '.tar.gz.part')
    if part.exists() and part.stat().st_size > release['bytes']:
        part.unlink()
    start = part.stat().st_size if part.exists() else 0
    if start != release['bytes']:
        report(f'Ligand library: downloading combined CCD + ChEBI archive ({release["bytes"] / 1e6:.0f} MB)')
        headers = {'User-Agent': 'ChimeraX-ligandsearch/1'}
        if start:
            headers['Range'] = f'bytes={start}-'
        with urllib.request.urlopen(urllib.request.Request(release['url'], headers=headers), timeout=30) as response:
            append = start > 0 and response.status == 206
            if append and response.headers.get('Content-Range') != f'bytes {start}-{release["bytes"]-1}/{release["bytes"]}':
                raise ValueError('Unexpected ligand-library download range; retry the download.')
            done = start if append else 0
            with part.open('ab' if append else 'wb') as stream:
                while chunk := response.read(1024*1024):
                    _cancelled(cancel)
                    done += len(chunk)
                    if done > release['bytes']:
                        raise ValueError('The library download exceeds its published size.')
                    stream.write(chunk)
                    report(f'Ligand library: downloading {done / 1e6:.1f}/{release["bytes"] / 1e6:.1f} MB')
        if done != release['bytes']:
            raise ValueError('The library download was interrupted; retry to resume it.')
    if _file_digest(part, cancel) != release['sha256']:
        part.unlink(missing_ok=True)
        raise ValueError('The library checksum does not match the release; retry the download.')
    install_library_archive(cache, part, cancel, report)
    part.unlink(missing_ok=True)
    return True


def _public_json(url, cancel, payload=None, timeout=30, attempts=3):
    """Retry transient public-data errors, without adding a third-party client."""
    import urllib.error
    import urllib.request

    data = None if payload is None else json.dumps(payload).encode('utf-8')
    headers = {'User-Agent': 'ChimeraX-ligandsearch/1', 'Content-Type': 'application/json'}
    for attempt in range(attempts):
        _cancelled(cancel)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=headers), timeout=timeout) as response:
                return json.load(response)
        except (OSError, ValueError) as error:
            if isinstance(error, urllib.error.HTTPError) and error.code not in (429, 500, 502, 503, 504):
                raise
            if attempt == attempts - 1:
                raise
            cancel.wait(2 ** attempt)


def _entry_components(entry):
    """Count a component once per PDB entry, not once per symmetry/ligand copy."""
    components = set()
    for group, container, field, singular in (
        ('nonpolymer_entities', 'rcsb_nonpolymer_entity_container_identifiers', 'nonpolymer_comp_id', True),
        ('branched_entities', 'rcsb_branched_entity_container_identifiers', 'chem_comp_monomers', False),
        ('polymer_entities', 'rcsb_polymer_entity_container_identifiers', 'chem_comp_nstd_monomers', False),
    ):
        for entity in entry.get(group) or ():
            values = (entity.get(container) or {}).get(field)
            if values:
                components.update([values] if singular else values)
    return components


def update_frequencies(cache, cancel=None, report=lambda message: None):
    """Cache RCSB's bulk entry/component mapping; publish only complete counts."""
    from collections import Counter
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    import uuid

    cancel = cancel or threading.Event()
    report('Counting PDB occurrences: fetching entry list')
    ids = _public_json(PDB_IDS_URL, cancel)
    if not isinstance(ids, list) or not ids or any(not isinstance(value, str) for value in ids):
        raise ValueError('RCSB returned an invalid PDB entry list; previous frequency cache retained.')
    ids = sorted(set(ids))
    # Follow RCSB's ligand-composition recipe, including branched components and
    # nonstandard polymer monomers. Entry summary lists omit some bound ligands.
    fields = ('rcsb_id nonpolymer_entities{rcsb_nonpolymer_entity_container_identifiers{nonpolymer_comp_id}} '
              'branched_entities{rcsb_branched_entity_container_identifiers{chem_comp_monomers}} '
              'polymer_entities{rcsb_polymer_entity_container_identifiers{chem_comp_nstd_monomers}}')
    def fetch(batch):
        query = '{entries(entry_ids:' + json.dumps(batch) + '){' + fields + '}}'
        response = _public_json(PDB_DATA_URL, cancel, {'query': query})
        entries = (response.get('data') or {}).get('entries')
        if response.get('errors') or not isinstance(entries, list) or len(entries) != len(batch) or \
                {entry.get('rcsb_id') for entry in entries if entry} != set(batch):
            raise ValueError('RCSB returned incomplete PDB occurrence data; previous frequency cache retained.')
        return entries

    batches = iter([ids[start:start+1000] for start in range(0, len(ids), 1000)])
    pool, counts, done = ThreadPoolExecutor(max_workers=4), Counter(), 0
    pending = set()
    try:
        def dispatch():
            batch = next(batches, None)
            if batch:
                pending.add(pool.submit(fetch, batch))
        for _ in range(4):
            dispatch()
        while pending:
            _cancelled(cancel)
            ready, _ = wait(pending, timeout=.1, return_when=FIRST_COMPLETED)
            for future in ready:
                pending.remove(future)
                entries = future.result()
                for entry in entries:
                    counts.update(_entry_components(entry))
                done += len(entries)
                report(f'Counting PDB occurrences: {done:,}/{len(ids):,} entries')
                dispatch()
        _cancelled(cancel)
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
    record = dict(schema=1, source=PDB_DATA_URL, entries=len(ids), counts=dict(counts),
                  created=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    cache = Path(cache).expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)
    temporary = cache / ('frequency-' + uuid.uuid4().hex + '.building')
    try:
        temporary.write_text(json.dumps(record, separators=(',', ':')))
        _cancelled(cancel)
        os.replace(temporary, cache / 'frequencies.json')
    finally:
        temporary.unlink(missing_ok=True)
    report(f'PDB frequencies ready: {len(counts):,} components in {len(ids):,} entries')
    return record


def frequency_info(cache):
    path = Path(cache).expanduser().resolve() / 'frequencies.json'
    if not path.is_file():
        return None
    try:
        record = json.loads(path.read_text())
        if record['schema'] != 1 or not isinstance(record['counts'], dict) or \
                not isinstance(record['entries'], int) or record['entries'] <= 0:
            return None
        if any(not isinstance(code, str) or isinstance(count, bool) or not isinstance(count, int)
               or not 0 <= count <= record['entries'] for code, count in record['counts'].items()):
            return None
        return record
    except (OSError, ValueError, KeyError, TypeError):
        return None


def library_order(db, volume, fraction, frequencies=None):
    """Soft size match times a capped frequency boost; never exclude rare CCDs."""
    import numpy as np

    rows = db.execute('SELECT id,v0,v1,v2 FROM components').fetchall()
    priorities = []
    for ident, *volumes in rows:
        estimate = float(np.interp(fraction, VOLUME_LEVELS, volumes))
        priority = math.exp(-.5 * (math.log(max(estimate, 1e-8) / volume) / math.log(2)) ** 2)
        # Common ligands get at most 3x priority. A sufficiently better size
        # match still wins; zero-occurrence components retain their full prior.
        if frequencies is not None:
            priority *= 1 + min(2., .3 * math.log1p(frequencies.get(ident, 0)))
        priorities.append((priority, ident, estimate))
    priorities.sort(key=lambda row: (-row[0], row[1]))
    return priorities


def parse_groups(value='all'):
    import re
    values = value.split(',') if isinstance(value, str) else list(value)
    keys = tuple(dict.fromkeys(str(v).strip().lower().replace(' ', '_') for v in values))
    allowed = {'all', 'unclassified'} | {row[0] for row in CHEMICAL_GROUPS}
    unknown = {key for key in keys if key not in allowed and not re.fullmatch(r'chemont:[0-9]{7}', key)}
    if not keys or unknown or ('all' in keys and len(keys) > 1):
        raise ValueError('groups must be all, unclassified, or comma-separated chemical-group keys: '
                         + ', '.join(row[0] for row in CHEMICAL_GROUPS) + ', or chemont: followed by a seven-digit class ID. Use all by itself.')
    return keys


def parse_group_match(value):
    value = str(value).strip().lower()
    if value not in ('any', 'all'):
        raise ValueError('group_match must be any (OR) or all (AND).')
    return value


def hit_matches_groups(hit, groups, match='any'):
    """Filter a retained hit by its recorded tags, never by a guessed category."""
    if groups == ('all',):
        return True
    tags = hit.get('chemical_groups')
    if tags is None:
        return False  # Missing annotations do not establish group membership.
    matches = [not tags if key == 'unclassified' else key in tags for key in groups]
    return all(matches) if match == 'all' else any(matches)


def _component_signature(ids):
    import hashlib
    return hashlib.sha256('\n'.join(sorted(ids)).encode('ascii')).hexdigest()


def _group_definition():
    import hashlib
    return hashlib.sha256(repr((CHEMICAL_GROUPS, CHEMONT_PRESETS)).encode('ascii')).hexdigest()


def _group_ontology(cache, cancel, report, refresh=False):
    import urllib.error
    import urllib.request
    import uuid

    path = cache / 'chebi_core.obo.gz'
    if path.is_file() and not refresh:
        return path
    tag_path = cache / 'chebi_etag.json'
    try:
        etag = json.loads(tag_path.read_text()).get('etag')
    except (OSError, ValueError):
        etag = None
    headers = {'User-Agent': 'ChimeraX-ligandsearch/1'}
    if path.is_file() and etag:
        headers['If-None-Match'] = etag
    temporary = cache / ('chebi-' + uuid.uuid4().hex + '.part')
    try:
        try:
            response = urllib.request.urlopen(urllib.request.Request(CHEBI_URL, headers=headers), timeout=30)
        except urllib.error.HTTPError as error:
            if error.code == 304 and path.is_file():
                return path
            raise
        with response, temporary.open('wb') as stream:
            done = 0
            while True:
                _cancelled(cancel)
                chunk = response.read(1024*1024)
                if not chunk:
                    break
                stream.write(chunk)
                done += len(chunk)
                report(f'Chemical groups: downloading ChEBI, {done / 2**20:.1f} MB')
            expected = int(response.headers.get('Content-Length', 0))
            if expected and done != expected:
                raise RuntimeError('ChEBI download was incomplete; retry Update chemical groups.')
            etag = response.headers.get('ETag')
        os.replace(temporary, path)
        tag_path.write_text(json.dumps(dict(etag=etag)))
        return path
    finally:
        temporary.unlink(missing_ok=True)


def ontology_group_masks(stream, cancel=None, identities=None, molecules=None):
    """Propagate preset tags down is-a, selected charge/tautomer, and role edges.

    has-part and functional-parent relations are deliberately not class membership.
    A compound containing a steroid or sugar is not automatically that class.
    """
    from collections import defaultdict, deque

    cancel = cancel or threading.Event()
    nodes, aliases, node = {}, {}, None
    version = ''
    roles = {'RO:0000087', 'has_role'}
    charge = {'RO:0018033', 'RO:0018034', 'RO:0018036',
              'is_conjugate_acid_of', 'is_conjugate_base_of', 'is_tautomer_of'}
    def finish():
        if node is None or not node.get('id'):
            return
        ident = node['id']
        if node.get('obsolete'):
            if node.get('replacement'):
                aliases[ident] = node['replacement']
            return
        nodes[ident] = node
        if identities is not None and node.get('inchikey'):
            identities.setdefault(node['inchikey'], set()).add(ident)
        if molecules is not None and node.get('inchikey'):
            molecules[ident] = {key: node[key] for key in ('name', 'inchikey', 'formula', 'charge', 'smiles') if key in node}
        for alternate in node['aliases']:
            aliases[alternate] = ident
    for index, line in enumerate(stream):
        if index % 10000 == 0:
            _cancelled(cancel)
        line = line.strip()
        if line.startswith('data-version: '):
            version = line[14:]
        elif line.startswith('['):
            finish()
            node = dict(parents=[], roles=[], charge=[], aliases=[]) if line == '[Term]' else None
        elif node is not None:
            if line.startswith('id: '):
                node['id'] = line[4:]
            elif line.startswith('name: '):
                node['name'] = line[6:]
            elif line.startswith('alt_id: '):
                node['aliases'].append(line[8:])
            elif line.startswith('property_value: '):
                fields = line.split('"')
                if len(fields) >= 3:
                    property_name = line.split()[1].rsplit('/', 1)[-1].rsplit(':', 1)[-1].lower()
                    key = {'inchi_key_string': 'inchikey', 'inchikey': 'inchikey',
                           'smiles_string': 'smiles', 'smiles': 'smiles',
                           'generalized_empirical_formula': 'formula', 'formula': 'formula',
                           'formal_charge': 'charge', 'charge': 'charge'}.get(property_name)
                    if key:
                        node[key] = fields[1]
            elif line.startswith('is_a: '):
                node['parents'].append(line[6:].split()[0])
            elif line == 'is_obsolete: true':
                node['obsolete'] = True
            elif line.startswith('replaced_by: '):
                node['replacement'] = line[13:].split()[0]
            elif line.startswith('relationship: '):
                fields = line[14:].split()
                if len(fields) >= 2 and fields[0] in roles | charge:
                    node['roles' if fields[0] in roles else 'charge'].append(fields[1])
    finish()
    def canonical(ident):
        seen = set()
        while ident in aliases and ident not in seen:
            seen.add(ident)
            ident = aliases[ident]
        return ident
    aliases = {key: canonical(key) for key in aliases}
    role_mask = sum(1 << i for i, row in enumerate(CHEMICAL_GROUPS) if row[3] == 'role')
    all_mask = (1 << len(CHEMICAL_GROUPS)) - 1
    edges, masks, pending = defaultdict(list), {}, deque()
    for ident, data in nodes.items():
        for parent in data['parents']:
            edges[canonical(parent)].append((ident, all_mask))
        for parent in data['roles']:
            edges[canonical(parent)].append((ident, role_mask))
        for peer in data['charge']:
            peer = canonical(peer)
            if peer in nodes:
                edges[peer].append((ident, all_mask))
                edges[ident].append((peer, all_mask))
    for index, (key, label, roots, kind) in enumerate(CHEMICAL_GROUPS):
        for root in roots:
            root = canonical(root)
            if root not in nodes:
                raise ValueError(f'ChEBI is missing the required {label} root; previous group cache retained.')
            masks[root] = masks.get(root, 0) | (1 << index)
            pending.append(root)
    steps = 0
    while pending:
        ident = pending.popleft()
        steps += 1
        if steps % 4096 == 0:
            _cancelled(cancel)
        for child, allowed in edges[ident]:
            value = masks.get(child, 0) | (masks[ident] & allowed)
            if value != masks.get(child, 0):
                masks[child] = value
                pending.append(child)
    for ident in nodes:
        masks.setdefault(ident, 0)  # Distinguish known-but-unassigned terms from missing mappings.
    return masks, aliases, version


def identity_group_lookup(identities, masks):
    """Index exact identifiers and conservative connectivity-level tag agreement.

    Connectivity can join different stereochemical/charge forms. Retain only
    tags shared by every matching term, rather than borrowing one form's role.
    """
    all_mask = (1 << len(CHEMICAL_GROUPS)) - 1
    connectivity = {}
    for key, terms in identities.items():
        prefix = key.split('-')[0]
        value = all_mask
        for ident in terms:
            value &= masks.get(ident, 0)
        if prefix in connectivity:
            old_value, old_terms = connectivity[prefix]
            connectivity[prefix] = (old_value & value, old_terms | terms)
        else:
            connectivity[prefix] = (value, set(terms))
    return connectivity


def identity_group_evidence(key, identities, connectivity, masks):
    if not key:
        return []
    if key in identities:
        terms = identities[key]
        method = 'ChEBI exact InChIKey'
        value = 0
        for ident in terms:
            value |= masks.get(ident, 0)
    else:
        value, terms = connectivity.get(key.split('-')[0], (0, set()))
        method = 'ChEBI connectivity consensus'
    return [(group, method, ident) for i, (group, label, roots, kind) in enumerate(CHEMICAL_GROUPS)
            if value & (1 << i) for ident in sorted(terms) if masks.get(ident, 0) & (1 << i)]


def chemont_nodes(result, include_alternatives=False):
    """Primary lineage; speculative alternative matches are opt-in diagnostics."""
    if result.get('error'):
        return {}
    nodes = [result.get(key) for key in ('kingdom', 'superclass', 'class', 'subclass', 'direct_parent')]
    nodes += result.get('intermediate_nodes') or []
    if include_alternatives:
        nodes += result.get('alternative_parents') or []
    return {node['chemont_id']: node['name'] for node in nodes if isinstance(node, dict)
            and isinstance(node.get('name'), str) and isinstance(node.get('chemont_id'), str)
            and node['chemont_id'].startswith('CHEMONTID:')}


def classify_structures(cache, descriptors, cancel=None, report=lambda message: None):
    """Classify public CCD SMILES in batches, with a resumable versioned cache.

    This supplies structural families independently of ChEBI database coverage.
    It never predicts biological roles or changes the fitting candidate scores.
    """
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    import urllib.error

    cancel = cancel or threading.Event()
    version = _public_json(CHEMONT_SERVICE + '/version', cancel)
    stamp = version.get('rules_git_sha')
    if not isinstance(stamp, str) or not stamp or not version.get('engine_version'):
        raise ValueError('The chemical classifier did not report its rule version; retry Update chemical groups.')
    stamp = version['engine_version'] + ':' + stamp
    limits = version.get('limits') or {}
    batch_size = min(100, max(1, int(limits.get('max_batch', 100))))
    max_length = int(limits.get('max_smiles_length', 2000))
    chosen, results = {}, {}
    for code, row in descriptors.items():
        if row.get('elements') and 6 not in row['elements']:
            # Carbon-free CCD components are definitively inorganic. Do not ask
            # an organic-structure classifier to parse metal coordination cages.
            results[code] = dict(kingdom=dict(name='Inorganic compounds', chemont_id='CHEMONTID:0000001'),
                                 classification_source='CCD elemental composition')
            chosen[code] = None
            continue
        alternatives = [row.get('SMILES_stereo'), row.get('SMILES')]
        chosen[code] = next((s for s in alternatives if isinstance(s, str) and 0 < len(s) <= max_length), None)
        if chosen[code] is None:
            results[code] = dict(error={'code': 'missing_structure', 'message': 'No supported CCD SMILES available.'})
    with closing(sqlite3.connect(Path(cache) / 'chemont.sqlite')) as db:
        db.execute('CREATE TABLE IF NOT EXISTS classifications('
                   'smiles TEXT,version TEXT,payload TEXT,PRIMARY KEY(smiles,version))')
        stored = {s: json.loads(p) for s, p in db.execute('SELECT smiles,payload FROM classifications WHERE version=?', (stamp,))}
        pending_smiles = sorted({s for s in chosen.values() if s and s not in stored})
        completed = len({s for s in chosen.values() if s}) - len(pending_smiles)
        total = completed + len(pending_smiles)
        report(f'Chemical groups: structural classification {completed:,}/{total:,} structures cached')

        def fetch(batch):
            _cancelled(cancel)
            try:
                response = _public_json(CHEMONT_SERVICE + '/classify', cancel, {'smiles': batch}, timeout=35, attempts=1)
            except urllib.error.HTTPError as error:
                if error.code in (413, 504) and len(batch) > 1:
                    middle = len(batch) // 2
                    report(f'Chemical groups: isolating a slow {len(batch)}-structure batch')
                    return fetch(batch[:middle]) + fetch(batch[middle:])
                if error.code == 504 and len(batch) == 1:
                    return [dict(smiles=batch[0], error=dict(code='classification_timeout',
                                                          message='This structure exceeded the classifier time limit.'))]
                raise
            except TimeoutError:
                if len(batch) == 1:
                    return [dict(smiles=batch[0], error=dict(code='classification_timeout',
                                                          message='This structure exceeded the classifier time limit.'))]
                middle = len(batch) // 2
                report(f'Chemical groups: retrying a slow batch as {middle}-structure batches')
                return fetch(batch[:middle]) + fetch(batch[middle:])
            rows = response.get('results')
            if not isinstance(rows, list) or len(rows) != len(batch) or any(
                    not isinstance(row, dict) or row.get('smiles') != s for s, row in zip(batch, rows)):
                raise ValueError('The chemical classifier returned incomplete or mismatched structures; previous tags retained.')
            return rows

        # Two outstanding batches keep throughput high without flooding the service.
        # Isolate highly connected cages/very long structures: one expensive
        # rule evaluation must not hold an otherwise quick batch hostage.
        complex_smiles = {s for s in pending_smiles if len(s) > 1200 or sum(c.isdigit() for c in s) > 32}
        ordinary = [s for s in pending_smiles if s not in complex_smiles]
        batches = iter([ordinary[i:i+batch_size] for i in range(0, len(ordinary), batch_size)] +
                       [[s] for s in sorted(complex_smiles)])
        pool, pending = ThreadPoolExecutor(max_workers=2), set()
        def dispatch():
            batch = next(batches, None)
            if batch:
                pending.add(pool.submit(fetch, batch))
        try:
            for _ in range(2):
                dispatch()
            while pending:
                _cancelled(cancel)
                ready, _ = wait(pending, timeout=.1, return_when=FIRST_COMPLETED)
                for future in ready:
                    pending.remove(future)
                    for row in future.result():
                        s = row['smiles']
                        stored[s] = row
                        if (row.get('error') or {}).get('code') != 'internal_error':
                            db.execute('INSERT OR REPLACE INTO classifications VALUES(?,?,?)', (s, stamp, json.dumps(row)))
                        completed += 1
                    db.commit()
                    report(f'Chemical groups: structural classification {completed:,}/{total:,} structures checked')
                    dispatch()
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        for code, s in chosen.items():
            if s:
                results[code] = stored[s]
    # CCD supplies more than one independently generated representation. Retry
    # parser failures with its non-stereochemical descriptor, never edited text.
    fallback = {code: {'SMILES': row['SMILES']} for code, row in descriptors.items()
                if (results[code].get('error') or {}).get('code') == 'invalid_smiles'
                and isinstance(row.get('SMILES'), str) and row['SMILES'] != chosen[code]
                and 0 < len(row['SMILES']) <= max_length}
    if fallback:
        report(f'Chemical groups: retrying {len(fallback):,} alternate CCD descriptors')
        retried, other = classify_structures(cache, fallback, cancel, report)
        if (other.get('engine_version'), other.get('rules_git_sha')) != (version['engine_version'], version['rules_git_sha']):
            raise ValueError('Chemical classifier rules changed during the update; retry Update chemical groups.')
        for code, row in retried.items():
            if not row.get('error'):
                results[code] = dict(row, descriptor_fallback=True)
    return results, version


def update_groups(cache, cancel=None, report=lambda message: None, refresh=False, ontology=None):
    """Build an independent, atomic classification cache for the current CCD IDs."""
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    import gzip
    import uuid
    import zlib

    cache = Path(cache).expanduser().resolve()
    if library_info(cache) is None:
        raise ValueError('Build the CCD library before updating chemical groups.')
    cancel = cancel or threading.Event()
    report('Chemical groups: reading the CCD index')
    types, compositions, nitrogen = {}, {}, []
    with closing(sqlite3.connect((cache / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)) as source:
        for ident, payload in source.execute('SELECT id,payload FROM components'):
            record = json.loads(payload)
            types[ident] = record.get('type', '')
            compositions[ident] = sorted(set(record.get('elements') or ()))
            if has_carbon_nitrogen_bond(record):
                nitrogen.append(ident)
    ids = sorted(types)
    identities, retained = {}, {}
    try:
        ontology = Path(ontology) if ontology is not None else _group_ontology(cache, cancel, report, refresh)
        report('Chemical groups: preparing ChEBI class hierarchy')
        with gzip.open(ontology, 'rt', encoding='utf-8') as stream:
            masks, aliases, version = ontology_group_masks(stream, cancel, identities)
        if not identities:
            raise ValueError('ChEBI structure identifiers are missing.')
    except InterruptedError:
        raise
    except (OSError, ValueError, EOFError, zlib.error, RuntimeError) as error:
        if isinstance(error, (gzip.BadGzipFile, EOFError, zlib.error)) and ontology == cache / 'chebi_core.obo.gz':
            ontology.unlink(missing_ok=True)
        masks, aliases, identities, version = {}, {}, {}, 'unavailable'
        previous = cache / 'groups.sqlite'
        if previous.is_file():
            try:
                with closing(sqlite3.connect(previous.as_uri() + '?mode=ro', uri=True)) as old:
                    for ident, key, source, source_id in old.execute(
                            "SELECT comp_id,group_key,source,source_id FROM tags WHERE source LIKE 'ChEBI %'"):
                        if ident in types:
                            retained.setdefault(ident, []).append((ident, key, source, source_id))
                    meta = dict(old.execute('SELECT key,value FROM metadata'))
                    version = 'retained:' + meta.get('ontology', 'unknown').removeprefix('retained:')
            except sqlite3.Error:
                retained = {}
        report(f'Warning: ChEBI update unavailable ({error}); structural classification continues. '
               'Previously cached ChEBI tags are retained where available.')
    connectivity = identity_group_lookup(identities, masks)
    temporary = cache / ('groups-' + uuid.uuid4().hex + '.building')
    db, pool = sqlite3.connect(temporary), ThreadPoolExecutor(max_workers=4)
    try:
        db.executescript('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);'
                         'CREATE TABLE components(id TEXT PRIMARY KEY,mapped INTEGER,available INTEGER,chebi_ids TEXT,structural_parent TEXT);'
                         'CREATE TABLE group_labels(key TEXT PRIMARY KEY,label TEXT);'
                         'CREATE TABLE tags(comp_id TEXT,group_key TEXT,source TEXT,source_id TEXT,'
                         'PRIMARY KEY(comp_id,group_key,source,source_id));'
                         'CREATE INDEX tag_groups ON tags(group_key,comp_id);')
        def fetch(batch):
            query = ('{chem_comps(comp_ids:' + json.dumps(batch) + '){rcsb_id '
                     'rcsb_chem_comp_descriptor{InChIKey SMILES SMILES_stereo} '
                     'rcsb_chem_comp_related{resource_name resource_accession_code}}}')
            response = _public_json(PDB_DATA_URL, cancel, {'query': query})
            rows = (response.get('data') or {}).get('chem_comps')
            if response.get('errors') or not isinstance(rows, list):
                raise ValueError('RCSB chemical mappings are unavailable; previous tags retained.')
            found = [row['rcsb_id'] for row in rows if row is not None]
            if len(found) != len(set(found)) or not set(found) <= set(batch):
                raise ValueError('RCSB returned inconsistent chemical mappings; previous tags retained.')
            return batch, {row['rcsb_id']: row for row in rows if row is not None}
        batches = iter([ids[i:i+500] for i in range(0, len(ids), 500)])
        pending, done, available, descriptors = set(), 0, 0, {}
        def dispatch():
            batch = next(batches, None)
            if batch:
                pending.add(pool.submit(fetch, batch))
        for _ in range(4):
            dispatch()
        while pending:
            _cancelled(cancel)
            ready, _ = wait(pending, timeout=.1, return_when=FIRST_COMPLETED)
            for future in ready:
                pending.remove(future)
                batch, rows = future.result()
                for ident in batch:
                    row = rows.get(ident, {})
                    descriptors[ident] = dict(row.get('rcsb_chem_comp_descriptor') or {}, elements=compositions[ident])
                    chebi = {r['resource_accession_code'] for r in row.get('rcsb_chem_comp_related') or ()
                             if r.get('resource_name') == 'ChEBI' and r.get('resource_accession_code')}
                    resolved = {aliases.get(code, code) for code in chebi}
                    evidence = list(retained.get(ident, ()))
                    for code in chebi:
                        value = masks.get(aliases.get(code, code), 0)
                        evidence.extend((ident, key, 'ChEBI ' + kind, code)
                                        for i, (key, label, roots, kind) in enumerate(CHEMICAL_GROUPS) if value & (1 << i))
                    key = (row.get('rcsb_chem_comp_descriptor') or {}).get('InChIKey')
                    if not any(code in masks for code in resolved):
                        evidence.extend((ident, *item) for item in identity_group_evidence(key, identities, connectivity, masks))
                    matched = any(code in masks for code in resolved) or key in identities or bool(
                        key and key.split('-')[0] in connectivity)
                    db.execute('INSERT INTO components VALUES(?,?,?,?,?)',
                               (ident, int(matched), int(bool(row)), json.dumps(sorted(chebi)), None))
                    kind = types[ident].upper()
                    if 'PEPTIDE' in kind and 'LINKING' in kind:
                        evidence.append((ident, 'amino_acids', 'CCD monomer type', types[ident]))
                    if ('RNA' in kind or 'DNA' in kind) and ('LINKING' in kind or 'TERMINUS' in kind):
                        evidence.append((ident, 'nucleotides', 'CCD monomer type', types[ident]))
                    if 'SACCHARIDE' in kind:
                        evidence.append((ident, 'carbohydrates', 'CCD monomer type', types[ident]))
                    db.executemany('INSERT OR IGNORE INTO tags VALUES(?,?,?,?)', evidence)
                available += len(rows)
                done += len(batch)
                db.commit()
                report(f'Chemical groups: {done:,}/{len(ids):,} CCD entries checked')
                dispatch()
        _cancelled(cancel)
        if not available:
            raise ValueError('RCSB returned no chemical-component mappings; previous tags retained.')
        structural, classifier = classify_structures(cache, descriptors, cancel, report)
        for index, (ident, result) in enumerate(structural.items()):
            if index % 1000 == 0:
                _cancelled(cancel)
                report(f'Chemical groups: writing families {index:,}/{len(structural):,} CCD entries')
            nodes = chemont_nodes(result)
            origin = result.get('classification_source', 'ClassyFire2 primary lineage')
            parent = result.get('direct_parent') or result.get('class') or result.get('superclass') or result.get('kingdom')
            db.execute('UPDATE components SET structural_parent=? WHERE id=?', ((parent or {}).get('name'), ident))
            evidence = []
            for code, label in nodes.items():
                key = 'chemont:' + code.split(':')[1]
                db.execute('INSERT OR REPLACE INTO group_labels VALUES(?,?)', (key, label))
                evidence.append((ident, key, origin, code))
            for key, roots in CHEMONT_PRESETS.items():
                if key == 'organic_nitrogen':
                    continue  # Broad overlapping chemistry, not one primary taxonomy branch.
                evidence.extend((ident, key, origin, 'CHEMONTID:' + root)
                                for root in roots if 'CHEMONTID:' + root in nodes)
            db.executemany('INSERT OR IGNORE INTO tags VALUES(?,?,?,?)', evidence)
        db.executemany('INSERT OR IGNORE INTO tags VALUES(?,?,?,?)',
                       ((ident, 'organic_nitrogen', 'CCD bond connectivity', 'C-N') for ident in nitrogen))
        tagged = db.execute('SELECT count(DISTINCT comp_id) FROM tags').fetchone()[0]
        mapped = db.execute('SELECT sum(mapped) FROM components').fetchone()[0]
        metadata = dict(schema=GROUP_SCHEMA, nitrogen_rule='C-N-v1', definition=_group_definition(), signature=_component_signature(ids), total=str(len(ids)),
                        tagged=str(tagged), mapped=str(mapped), available=str(available), ontology=version, source=CHEBI_URL,
                        mapping_source=PDB_DATA_URL, classifier_source=CHEMONT_SERVICE, classifier=json.dumps(classifier),
                        structural=str(sum(bool(r.get('class')) for r in structural.values())),
                        structural_errors=str(sum(bool(r.get('error')) for r in structural.values())),
                        created=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        db.executemany('INSERT INTO metadata VALUES(?,?)', metadata.items())
        db.commit()
        db.close()
        db = None
        _cancelled(cancel)
        os.replace(temporary, cache / 'groups.sqlite')
        report(f'Chemical groups ready: {tagged:,}/{len(ids):,} entries tagged')
        return metadata
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
        if db is not None:
            db.close()
        temporary.unlink(missing_ok=True)


def group_info(cache, current_ids=None):
    path = Path(cache).expanduser().resolve() / 'groups.sqlite'
    if not path.is_file():
        return None
    try:
        repair_nitrogen_groups(cache)
        with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as db:
            metadata = dict(db.execute('SELECT key,value FROM metadata'))
            if metadata.get('schema') != GROUP_SCHEMA or metadata.get('definition') != _group_definition():
                return None
            counts = dict(db.execute('SELECT group_key,count(DISTINCT comp_id) FROM tags GROUP BY group_key'))
            labels = dict(db.execute('SELECT key,label FROM group_labels'))
            total = int(metadata['total'])
            if current_ids is None:
                source = Path(cache).expanduser().resolve() / 'ccd.sqlite'
                if not source.is_file():
                    return None
                with closing(sqlite3.connect(source.as_uri() + '?mode=ro', uri=True)) as library:
                    current_ids = [row[0] for row in library.execute('SELECT id FROM components')]
            current = metadata.get('signature') == _component_signature(current_ids)
            return dict(metadata, current=current, labels=labels,
                        counts=dict(counts, all=total, unclassified=total-int(metadata['tagged'])))
    except (OSError, sqlite3.Error, ValueError, KeyError):
        return None


def has_carbon_nitrogen_bond(record):
    elements = record['elements']
    return any({elements[a], elements[b]} == {6, 7} for a, b, _ in record['bonds'])


def repair_nitrogen_groups(cache):
    """One local, transactional upgrade; detailed ChemOnt classes are untouched."""
    cache = Path(cache).expanduser().resolve()
    with closing(sqlite3.connect(cache / 'groups.sqlite')) as db:
        metadata = dict(db.execute('SELECT key,value FROM metadata'))
        if metadata.get('nitrogen_rule') == 'C-N-v1':
            return
        if metadata.get('schema') != GROUP_SCHEMA or metadata.get('definition') != _group_definition():
            return
        with closing(sqlite3.connect((cache / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)) as source:
            ids = [ident for ident, payload in source.execute('SELECT id,payload FROM components')
                   if has_carbon_nitrogen_bond(json.loads(payload))]
        with db:
            db.execute("DELETE FROM tags WHERE group_key='organic_nitrogen'")
            db.executemany('INSERT OR IGNORE INTO tags VALUES(?,?,?,?)',
                           ((ident, 'organic_nitrogen', 'CCD bond connectivity', 'C-N') for ident in ids))
            tagged = db.execute('SELECT count(DISTINCT comp_id) FROM tags').fetchone()[0]
            db.executemany('INSERT OR REPLACE INTO metadata VALUES(?,?)',
                           [('nitrogen_rule', 'C-N-v1'), ('tagged', str(tagged))])


def group_members(db, groups, match='any'):
    groups = parse_groups(groups)
    match = parse_group_match(match)
    if groups == ('all',):
        return None
    selected = set()
    if match == 'all' and 'unclassified' in groups and len(groups) > 1:
        return selected  # Untagged compounds cannot also possess a selected tag.
    if 'unclassified' in groups:
        selected.update(row[0] for row in db.execute('SELECT id FROM components WHERE NOT EXISTS '
                                                    '(SELECT 1 FROM tags WHERE comp_id=components.id)'))
    keys = [key for key in groups if key != 'unclassified']
    if keys:
        query = 'SELECT comp_id FROM tags WHERE group_key IN (' + ','.join('?' for _ in keys) + ') GROUP BY comp_id'
        parameters = keys
        if match == 'all':
            query += ' HAVING count(DISTINCT group_key)=?'
            parameters = keys + [len(keys)]
        selected.update(row[0] for row in db.execute(query, parameters))
    return selected


def group_evidence(db, ident):
    component = db.execute('SELECT structural_parent FROM components WHERE id=?', (ident,)).fetchone()
    if component is None:
        return dict(chemical_groups=None, group_evidence=[], classification_status='unavailable')
    rows = db.execute('SELECT group_key,source,source_id FROM tags WHERE comp_id=? ORDER BY group_key,source,source_id', (ident,)).fetchall()
    return dict(chemical_groups=sorted({row[0] for row in rows}),
                structural_family=component[0],
                group_evidence=[dict(group=key, source=source, id=source_id) for key, source, source_id in rows],
                classification_status='tagged' if rows else 'unclassified')


def prepare_blob(data, center, resolution, level, radius=14.):
    """Isolate density near the seed, raising only an initially connected contour.

    The grid is read once. A monotonic search finds the lowest separating
    threshold while following one fixed local density peak, not another blob.
    The caller owns any displayed-contour changes.
    """
    import numpy as np
    from scipy.ndimage import label, distance_transform_edt, generate_binary_structure

    if not all(math.isfinite(v) and v > 0 for v in (resolution, level, radius)):
        raise ValueError('Resolution, contour level and region radius must be positive finite numbers.')
    if resolution < 2 * max(np.linalg.norm(data.ijk_to_xyz_transform.matrix[:, :3], axis=0)):
        raise ValueError('Resolution is finer than map sampling; supply the actual map resolution in Angstroms.')
    lower, upper = shape_api()._shape_cofr_crop_bounds(data, center, radius)
    if np.any(upper <= lower):
        raise ValueError('The picked point is outside this map; pick the density blob again.')
    if np.prod(upper-lower) > 4_000_000:
        raise ValueError('Map sampling is too fine for the search region; use a coarser-sampled map and pick the blob again.')
    matrix = data.matrix(tuple(lower), tuple(upper-lower)).astype(np.float32, copy=True)
    if not np.isfinite(matrix).all():
        raise ValueError('The selected map region contains nonfinite density values.')
    affine = data.ijk_to_xyz_transform.matrix[:, :3]
    metric = affine.T @ affine
    if not np.allclose(metric, np.diag(np.diag(metric)), atol=np.max(np.diag(metric)) * 1e-7):
        raise ValueError('This version needs an orthogonal map grid; resample the map onto an orthogonal grid first.')
    initial_level = level
    connectivity = generate_binary_structure(3, 1)
    labels, _ = label(matrix >= level, connectivity)
    k, j, i = np.indices(matrix.shape)
    points = data.ijk_to_xyz_transform.transform_points(
        np.column_stack((i.ravel(), j.ravel(), k.ravel())).astype(float) + lower)
    nearby = (np.linalg.norm(points-center, axis=1) <= 3) & (labels.ravel() > 0)
    if not np.any(nearby):
        raise ValueError('No blob within 3 A of the picked point at this contour. Adjust the displayed map contour and pick again.')
    ids, counts = np.unique(labels.ravel()[nearby], return_counts=True)
    blob = labels == ids[np.argmax(counts)]
    local = np.flatnonzero(nearby & blob.ravel())
    anchor = local[np.argmax(matrix.ravel()[local])]
    peak = float(matrix.ravel()[anchor])
    def touches_boundary(mask):
        return any(np.any(np.take(mask, edge, axis=axis)) for edge in (0, -1) for axis in range(3))
    if touches_boundary(blob):
        original = blob
        def at_level(value):
            components, _ = label(original & (matrix >= value), connectivity)
            return components == components.ravel()[anchor]
        values = np.unique(matrix[original & (matrix >= level) & (matrix <= peak)])
        lo, hi = 0, len(values)-1
        if touches_boundary(at_level(float(values[hi]))):
            raise ValueError('No contour isolates this blob from the region boundary before it disappears. '
                             'Pick farther inside the intended blob, or choose a more isolated density component.')
        while hi-lo > 1:
            middle = (lo+hi) // 2
            if touches_boundary(at_level(float(values[middle]))):
                lo = middle
            else:
                hi = middle
        # Just exclude the connecting voxels; jumping to the next distinct peak
        # value could unnecessarily shrink the displayed blob.
        level = float(np.nextafter(values[lo], np.float32(np.inf)))
        blob = at_level(level)
    steps = np.linalg.norm(affine, axis=0)
    outside = distance_transform_edt(~blob, sampling=steps[::-1])
    soft = np.exp(-.5 * (outside / (resolution / 3)) ** 2)
    soft[outside > resolution] = 0
    return dict(matrix=(matrix * soft).astype(np.float32), raw=matrix, mask=outside <= resolution,
                blob=blob, origin=np.array(data.ijk_to_xyz_transform * lower), step=np.array(data.step),
                rotation=np.array(data.rotation), cell_angles=np.array(data.cell_angles), resolution=resolution,
                level=level, volume=float(blob.sum() * abs(np.linalg.det(affine))),
                fraction=float(level / matrix[blob].max()), initial_level=initial_level,
                anchor=points[anchor].copy(), anchor_density=peak)


class BlobHighlight:
    """Reversible native connected-surface coloring, including recontouring."""
    def __init__(self):
        self.surface = None
        self.callback = self.recolor

    def select(self, surface, triangle=None, point=None, max_level=None, defer=False):
        import numpy as np

        self.clear()
        self.surface = surface
        self.base_recolor = surface.auto_recolor_vertices
        self.point = (np.mean(surface.vertices[surface.triangles[triangle]], axis=0)
                      if point is None else np.array(point, float))
        self.max_level = max_level
        self.vertices, self.indices = None, None
        if not defer:
            self.paint(triangle)
        surface.auto_recolor_vertices = self.callback

    def paint(self, triangle=None):
        import numpy as np
        from chimerax.pick_blobs.pickblobs import connected_component, color_blob

        s = self.surface
        self.vertices = s.vertices
        self.indices = None
        if s.vertices is None or s.triangles is None or not len(s.triangles):
            return
        if self.max_level is not None and s.level > self.max_level:
            return  # The selected peak disappeared; do not jump to a neighbor.
        if triangle is None:
            delta = s.vertices - self.point
            vertex = int(np.argmin(np.einsum('ij,ij->i', delta, delta)))
            incident = np.flatnonzero(np.any(s.triangles == vertex, axis=1))
            if not len(incident):
                return
            triangle = int(incident[0])
        self.indices, _ = connected_component(s.triangles, triangle)
        self.base_color = np.array(s.color, copy=True)
        self.base_colors = None if s.vertex_colors is None else s.vertex_colors.copy()
        alpha = (self.base_color[3] if self.base_colors is None else self.base_colors[self.indices, 3].copy())
        color_blob(s, self.indices, (*BLOB_HIGHLIGHT_RGB, int(self.base_color[3])))
        colors = s.vertex_colors
        colors[self.indices, 3] = alpha
        s.vertex_colors = colors
        self.applied = colors[self.indices].copy()

    def recolor(self):
        if self.base_recolor is not None:
            self.base_recolor()
        self.paint()
        self.surface.auto_recolor_vertices = self.callback

    def clear(self):
        import numpy as np

        s = self.surface
        self.surface = None
        if s is None or s.deleted:
            return
        current_recolor = s.auto_recolor_vertices
        owned = current_recolor is self.callback
        if s.vertices is self.vertices and self.indices is not None and s.vertex_colors is not None:
            colors = s.vertex_colors.copy()
            unchanged = np.all(colors[self.indices] == self.applied, axis=1)
            indices = self.indices[unchanged]
            colors[indices] = self.base_color if self.base_colors is None else self.base_colors[indices]
            # Preserve unrelated manual color edits made while the tool was open.
            uniform = self.base_colors is None and np.all(colors == s.color)
            s.vertex_colors = None if uniform else colors
        s.auto_recolor_vertices = self.base_recolor if owned else current_recolor


class TopHits:
    """A bounded, per-CCD leaderboard; completion order cannot break tied ranks."""
    def __init__(self, limit=20):
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
            raise ValueError('Top-hit buffer must contain 1 to 100 entries.')
        self.limit, self.rows = limit, []

    def add(self, row, replace_equal=False):
        if not math.isfinite(row['cc']) or not math.isfinite(row.get('score', row['cc'])):
            return
        previous = next((old for old in self.rows if old['id'] == row['id']), None)
        if previous is not None:
            old_score, new_score = previous.get('score', previous['cc']), row.get('score', row['cc'])
            if old_score > new_score or (old_score == new_score and not replace_equal):
                return
        self.rows = [old for old in self.rows if old['id'] != row['id']]
        self.rows.append(row)
        self.rows.sort(key=lambda item: (-item.get('score', item['cc']), item['id']))
        del self.rows[self.limit:]


class _Quiet:
    def status(self, *args, **kwargs):
        pass


def selected_conformers(record, observed=True):
    conformers = record['conformers'] if observed else [c for c in record['conformers'] if c['source'] == 'ideal']
    return conformers or record['conformers'][:1]


def validate_smiles(value):
    value = value.strip()
    if not value or len(value) > 2048 or not value.isascii() or any(c.isspace() or not c.isprintable() for c in value):
        raise ValueError('Enter one SMILES string without internal spaces (up to 2,048 characters).')
    return value


def smiles_component(smiles):
    """Use ChimeraX's built-in converter in the disposable fitting process."""
    from collections import Counter
    from types import SimpleNamespace
    import numpy as np
    from chimerax.atomic import initialize_atomic
    from chimerax.core import core_settings
    from chimerax.core.attributes import RegAttrManager
    from chimerax.core.nogui import UI
    from chimerax.core.session import Session
    from chimerax.smiles import fetch_smiles

    smiles = validate_smiles(smiles)
    session = Session(silent=True)
    session.ui = UI(session)
    core_settings.init(session)
    RegAttrManager(session)
    initialize_atomic(session)
    # The headless worker has no Toolshed providers. The native converter only
    # needs the SDF encoding before handing the file to ChimeraX's SDF reader.
    session.data_formats = {'sdf': SimpleNamespace(encoding='utf-8')}
    structures = []
    try:
        try:
            structures, _ = fetch_smiles(session, smiles)
        except Exception as error:
            raise ValueError('ChimeraX could not convert the SMILES to 3D. Check the string and network access, then retry.') from error
        if len(structures) != 1:
            raise ValueError('SMILES must convert to a single molecule.')
        model = structures[0]
        atoms = [atom for atom in model.atoms if atom.element.number > 1]
        xyz = np.array([atom.coord for atom in atoms])
        if len(atoms) < 4 or not np.isfinite(xyz).all() or np.linalg.matrix_rank(xyz-xyz.mean(axis=0)) < 2:
            raise ValueError('SMILES fitting needs a usable 3D structure with at least four heavy atoms.')
        index = {atom: i for i, atom in enumerate(atoms)}
        orders = {1: 'SING', 2: 'DOUB', 3: 'TRIP', 4: 'AROM', 1.5: 'AROM'}
        bonds, neighbors = [], [set() for _ in atoms]
        for bond in model.bonds:
            a, b = bond.atoms
            if a in index and b in index:
                i, j = index[a], index[b]
                bonds.append((i, j, orders.get(getattr(bond, 'order', None), 'UNKN')))
                neighbors[i].add(j)
                neighbors[j].add(i)
        seen, pending = {0}, [0]
        while pending:
            new = neighbors[pending.pop()] - seen
            pending.extend(new)
            seen.update(new)
        if len(seen) != len(atoms):
            raise ValueError('Enter a single connected molecule, without salt or solvent fragments.')
        counts = Counter(atom.element.name for atom in model.atoms)
        symbols = sorted(counts, key=lambda name: (name not in ('C', 'H'), name))
        formula = ' '.join(name + (str(counts[name]) if counts[name] > 1 else '') for name in symbols)
        return dict(id='LIG', name='Custom ligand', formula=formula, library_source='SMILES', smiles=smiles,
                    names=[atom.name for atom in atoms], elements=[atom.element.number for atom in atoms],
                    bonds=bonds, conformers=[dict(source='ChimeraX SMILES', xyz=xyz.tolist())])
    finally:
        for model in structures:
            model.delete()


def validate_conformer_limit(limit, maximum=MAX_CONFORMERS, name='conformers'):
    if isinstance(limit, bool) or not isinstance(limit, int) or not 2 <= limit <= maximum:
        raise ValueError(f'{name} must be an integer from 2 to {maximum} (total per ligand).')


def rotatable_bonds(record):
    """Conservative bridge-bond rotations using the existing heavy-atom CCD graph."""
    elements = record['elements']
    graph = [set() for _ in elements]
    orders = {}
    for a, b, order in record['bonds']:
        if not (0 <= a < len(graph) and 0 <= b < len(graph)) or a == b:
            raise ValueError('CCD bond indices are inconsistent with its atoms.')
        graph[a].add(b)
        graph[b].add(a)
        orders[tuple(sorted((a, b)))] = str(order).upper()

    def component(start, blocked):
        seen, pending = {start}, [start]
        while pending:
            a = pending.pop()
            for b in graph[a]:
                if tuple(sorted((a, b))) != blocked and b not in seen:
                    seen.add(b)
                    pending.append(b)
        return seen

    def pi_neighbors(a, exclude):
        return [b for b in graph[a] if b != exclude and
                orders[tuple(sorted((a, b)))] in ('DOUB', 'TRIP', 'AROM', 'DELO')]

    rotors = []
    for (a, b), order in sorted(orders.items()):
        if order != 'SING' or min(len(graph[a]), len(graph[b])) < 2:
            continue  # Terminal heavy-atom bonds cannot change the heavy-atom shape.
        if elements[a] not in (6, 7, 8, 15, 16) or elements[b] not in (6, 7, 8, 15, 16):
            continue  # Do not rotate metal coordination or unfamiliar bond chemistry.
        pa, pb = pi_neighbors(a, b), pi_neighbors(b, a)
        restricted = bool(pa and pb)
        for center, neighbor, pi in ((a, b, pa), (b, a, pb)):
            if elements[neighbor] == 7 and pi:
                restricted = True  # Amides, amidines, anilines and similar conjugated C-N bonds.
            if elements[neighbor] in (7, 8, 16) and any(elements[j] in (7, 8, 16) for j in pi):
                restricted = True  # Carbonyl/thiocarbonyl and sulfonyl/phosphoryl linkages.
        if restricted:
            continue
        left = component(a, (a, b))
        if b in left:
            continue  # Removing a ring bond does not disconnect its endpoints.
        right = component(b, (a, b))
        if len(left) < len(right):
            a, b, moving = b, a, left
        else:
            moving = right
        # Prefer changes that affect substantial fragments, not tiny tails.
        staggered = elements[a] == elements[b] == 6 and not pa and not pb
        angles = (120, -120) if staggered else (120, -120, 180, 60, -60)
        rotors.append(dict(a=a, b=b, moving=sorted(moving), angles=angles,
                           size=min(len(left), len(right))))
    rotors.sort(key=lambda r: (-r['size'], record['names'][r['a']], record['names'][r['b']]))
    return graph, rotors


def _rotate_fragment(xyz, rotor, angle):
    """Rotate one connected fragment, preserving bond lengths, angles and rings."""
    import numpy as np

    origin = xyz[rotor['a']]
    axis = xyz[rotor['b']] - origin
    length = np.linalg.norm(axis)
    if length < 1e-8:
        return False
    axis /= length
    indices = rotor['moving']
    offsets = xyz[indices] - origin
    theta = math.radians(angle)
    xyz[indices] = (origin + offsets * math.cos(theta) + np.cross(axis, offsets) * math.sin(theta)
                    + np.outer(offsets @ axis, axis) * (1 - math.cos(theta)))
    return True


def sparse_conformers(record, observed=True, limit=DEFAULT_CONFORMERS, maximum=MAX_CONFORMERS):
    """Bounded single/multiple torsion perturbations; no Cartesian product search."""
    import numpy as np
    from scipy.spatial import cKDTree
    from chimerax.geometry import align_points

    validate_conformer_limit(limit, maximum)
    seeds = selected_conformers(record, observed)
    result = [dict(seed) for seed in seeds]
    if not seeds or len(result) >= limit:
        return result[:limit]
    graph, rotors = rotatable_bonds(record)
    if not rotors:
        return result
    # Exclude bonded and 1-3 neighbors from the intraligand overlap test.
    excluded = set()
    for a, neighbors in enumerate(graph):
        for b in neighbors | {c for b in neighbors for c in graph[b]}:
            if a != b:
                excluded.add(tuple(sorted((a, b))))
    seed_xyz = [np.array(seed['xyz'], dtype=float) for seed in seeds]
    chosen = [np.array(seed['xyz'], dtype=float) for seed in result]
    rng, single_goal = np.random.default_rng(0), max(1, math.ceil((limit-len(seeds)) / 2))
    single_budget = len(seeds) * sum(len(rotor['angles']) for rotor in rotors)
    for attempt in range(limit * 32):
        if len(result) >= limit:
            break
        seed_index = attempt % len(seeds)
        original = seed_xyz[seed_index]
        xyz, changes = original.copy(), []
        if len(result)-len(seeds) < single_goal and attempt < single_budget:
            index = (attempt // len(seeds)) % len(rotors)
            rotor = rotors[index]
            angle = rotor['angles'][(attempt // (len(seeds)*len(rotors))) % len(rotor['angles'])]
            rotations = [(rotor, angle)]
        else:
            count = 1 if len(rotors) == 1 else int(rng.integers(2, min(3, len(rotors))+1))
            indices = sorted(rng.choice(len(rotors), count, replace=False))
            rotations = [(rotors[i], int(rng.choice(rotors[i]['angles']))) for i in indices]
        for rotor, angle in rotations:
            if _rotate_fragment(xyz, rotor, angle):
                changes.append(dict(bond=f'{record["names"][rotor["a"]]}-{record["names"][rotor["b"]]}', degrees=angle))
        if not changes:
            continue
        pairs = cKDTree(xyz).query_pairs(INTERNAL_OVERLAP, output_type='ndarray')
        pairs = np.array([pair for pair in pairs if tuple(pair) not in excluded], dtype=int).reshape((-1, 2))
        if len(pairs):
            distances = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
            starting = np.linalg.norm(original[pairs[:, 0]] - original[pairs[:, 1]], axis=1)
            if np.any(distances < starting - .05):
                continue  # Reject new/worsened overlaps, not pre-existing close contacts.
        if any(align_points(xyz, other)[1] < CONFORMER_RMSD for other in chosen):
            continue
        number = len(result)-len(seeds)+1
        result.append(dict(source=f'{seeds[seed_index]["source"]} torsion {number}', xyz=xyz.tolist(),
                           torsion_changes=changes))
        chosen.append(xyz)
    return result


class Fitter:
    """Per-process immutable target and shared shapealign preparation cache."""
    def __init__(self, config):
        import numpy as np
        from types import SimpleNamespace
        from chimerax.geometry import Place
        from chimerax.map_data import ArrayGridData

        self.config = config
        self.api = shape_api()
        self.owner = SimpleNamespace(scene_position=Place(), session=SimpleNamespace(logger=_Quiet()))
        data = ArrayGridData(config['matrix'], config['origin'], config['step'],
                             config['cell_angles'], config['rotation'])
        self.owner.data = data
        self.target = self.api._ShapeDensityInput(self.owner, None)
        self.target._shared_reference = True
        self.cache = self.api._ShapeBatchCache()
        self.cache.__enter__()

    def close(self):
        self.cache.__exit__(None, None, None)

    def source(self, xyz, elements):
        import numpy as np

        source = self.api._ShapeDensityInput.__new__(self.api._ShapeDensityInput)
        source.model, source.atoms, source.is_map = self.owner, None, False
        source.xyz, source.weights = np.array(xyz, float), np.array(elements, float)
        source.center = self.api._shape_center(source.xyz, source.weights)
        source.landmarks = source.xyz - source.center
        source.frame_landmarks = self.api._shape_geometry(source.landmarks, source.weights)
        source._cache, source._grid_cache_key, source._shared_reference = {}, object(), False
        return source

    def score(self, source, pose):
        import numpy as np
        from chimerax.map.molmap import add_gaussians
        from chimerax.map_data import ArrayGridData
        from chimerax.map_fit.fitmap import volume_values

        c = self.config
        matrix = np.zeros(c['matrix'].shape, np.float32)
        grid = ArrayGridData(matrix, c['origin'], c['step'], c['cell_angles'], c['rotation'])
        placed = pose.transform_points(source.xyz).astype(np.float32)
        add_gaussians(grid, placed, source.weights.astype(np.float32), c['resolution'] / (math.pi * math.sqrt(2)), 4)
        first, second = matrix[c['mask']].astype(float), c['raw'][c['mask']].astype(float)
        first -= first.mean()
        second -= second.mean()
        norm = np.linalg.norm(first) * np.linalg.norm(second)
        cc = float(np.dot(first, second) / norm) if norm > 1e-20 else -1.
        # Density completeness only: no receptor model or clash penalty is used.
        sampled = volume_values(placed, grid.xyz_to_ijk_transform, c['raw'])
        supported = float(np.mean(sampled >= c['level']))
        threshold = matrix.max() * c['fraction']
        predicted = matrix >= threshold if threshold > 0 else np.zeros(matrix.shape, bool)
        coverage = float(np.mean(predicted[c['blob']]))
        outside = float(np.count_nonzero(predicted & ~c['mask']) / max(1, np.count_nonzero(predicted)))
        return cc, supported, coverage, outside

    def fit(self, record, observed=True, fast=False, torsions=False, conformers=DEFAULT_CONFORMERS):
        import numpy as np
        from chimerax.geometry import Place
        from chimerax.core.errors import UserError

        options = dict(resolution=self.config['resolution'], radius=None, objective='density',
                       max_points=FAST_SAMPLES if fast else STANDARD_SAMPLES,
                       extra_starts=FAST_EXTRA_STARTS if fast else STANDARD_EXTRA_STARTS,
                       fallback_correlation=0., bandpass=False)
        best = first_fit_error = None
        tested, skipped = 0, 0
        variants = sparse_conformers(record, observed, conformers) if torsions else selected_conformers(record, observed)
        for conformer in variants:
            try:
                source = self.source(conformer['xyz'], record['elements'])
                _, stats = self.api._fit_shape_original_map_search(source, self.target, options, None, collect=True)
            except (UserError, self.api._ShapeEmptyMapSearch) as error:
                # One unusable conformer should not discard another's valid fit.
                if first_fit_error is None:
                    first_fit_error = error
                skipped += 1
                continue
            tested += 1
            for row in stats['candidates']:
                pose = Place(matrix=np.array(row['matrix']))
                cc, supported, coverage, outside = self.score(source, pose)
                # Penalize incomplete blob explanations and unsupported ligand
                # atoms softly. This is density support, never a receptor clash.
                score = cc * math.sqrt(supported * coverage) if cc > 0 else cc
                if best is None or score > best['score']:
                    best = dict(id=record['id'], name=record['name'], formula=record['formula'], cc=cc,
                                score=score, fast_search=bool(fast),
                                torsion_sampling=bool(torsions),
                                torsion_changes=conformer.get('torsion_changes', []),
                                pose=pose.matrix.tolist(), xyz=conformer['xyz'], source=conformer['source'],
                                names=record['names'], elements=record['elements'], bonds=record['bonds'],
                                supported=supported, coverage=coverage, outside=outside)
        if best is None:
            if first_fit_error is not None:
                raise first_fit_error
            raise ValueError('No usable coordinate conformer.')
        best.update(conformers_tested=tested, conformers_skipped=skipped)
        best.update({key: record[key] for key in ('library_source', 'chebi_ids', 'pubchem_cid', 'inchikey', 'formal_charge', 'smiles')
                     if key in record})
        return best

    def local_pose(self, source, initial):
        """Native rigid fit of simulated density, warm-started in map coordinates."""
        from chimerax.map_fit.fitmap import locate_maximum

        # Shape grids are centered; undo that centering before applying the
        # already-fitted pose. No PCA/FFT/global orientation search is repeated.
        points, weights = self.api._shape_density_samples(source.grid(self.config['resolution']), STANDARD_SAMPLES)
        points += source.center
        delta, _ = locate_maximum(initial.transform_points(points), weights, self.config['raw'],
                                  self.owner.data.xyz_to_ijk_transform, max_steps=240,
                                  ijk_step_size_min=.01, ijk_step_size_max=.5, metric='sum product')
        return delta * initial

    def refine(self, hit, conformers=FINAL_CONFORMERS):
        """Broader sparse sampling around the winning conformation, not small-angle descent."""
        import numpy as np
        from chimerax.geometry import Place
        from chimerax.core.errors import UserError

        validate_conformer_limit(conformers, MAX_FINAL_CONFORMERS, 'refine_conformers')
        record = dict(hit, conformers=[dict(source='ideal', xyz=hit['xyz'])])
        variants = sparse_conformers(record, limit=conformers, maximum=MAX_FINAL_CONFORMERS)
        best, initial = dict(hit), Place(matrix=np.array(hit['pose'], float))
        original_score = hit.get('score', hit['cc'])
        tested = skipped = 0
        for variant in variants:
            source = self.source(variant['xyz'], hit['elements'])
            try:
                fitted = self.local_pose(source, initial)
            except (UserError, self.api._ShapeEmptyMapSearch):
                skipped += 1
                continue
            tested += 1
            # The local sum-product optimizer and final completeness-weighted
            # score differ. Compare both poses; never discard a better seed.
            for pose in (initial, fitted):
                cc, supported, coverage, outside = self.score(source, pose)
                score = cc * math.sqrt(supported*coverage) if cc > 0 else cc
                if math.isfinite(score) and score > best.get('score', best['cc']) + 1e-10:
                    changes = variant.get('torsion_changes', [])
                    best.update(score=score, cc=cc, supported=supported, coverage=coverage, outside=outside,
                                xyz=variant['xyz'], pose=pose.matrix.tolist(),
                                source=hit['source'] + ('; final sparse torsions' if changes else '; final local fit'),
                                torsion_changes=list(hit.get('torsion_changes', [])) + changes,
                                final_torsion_changes=changes)
        best.update(final_refined=True, refinement_improved=best.get('score', best['cc']) > original_score + 1e-10,
                    pre_refine_score=original_score, pre_refine_cc=hit['cc'],
                    refine_conformers_tested=tested, refine_conformers_skipped=skipped,
                    refine_conformer_limit=conformers,
                    starting_pose=hit['pose'], starting_xyz=hit['xyz'], starting_source=hit['source'])
        return best


def worker_main(directory):
    import sys
    import numpy as np

    root = Path(directory)
    with np.load(root / 'blob.npz', allow_pickle=False) as values:
        config = {key: values[key] for key in values.files}
    for key in ('resolution', 'level', 'volume', 'fraction'):
        config[key] = float(config[key])
    fitter = Fitter(config)
    try:
        for line in sys.stdin:
            job = json.loads(line)
            try:
                if 'refine_hit' in job:
                    hit = fitter.refine(job['refine_hit'], job['refine_conformers'])
                else:
                    record = smiles_component(job['smiles']) if 'smiles' in job else job['component']
                    hit = fitter.fit(record, job.get('observed', True), fast=job.get('fast', False),
                                     torsions=job.get('torsions', False),
                                     conformers=job.get('conformers', DEFAULT_CONFORMERS))
                result = dict(result=hit)
            except Exception as error:
                result = dict(error=f'{type(error).__name__}: {error}')
            result['id'] = job['component']['id']
            print('LIGAND_RESULT:' + json.dumps(result), flush=True)
    finally:
        fitter.close()


class SearchJob:
    """Bounded worker dispatch. Stop terminates outstanding work, not saved hits."""
    def __init__(self, cache, config, workers=0, observed=True, ids=None, frequency=True, fast=False,
                 torsions=False, conformers=DEFAULT_CONFORMERS, groups='all', group_match='any', sources='ccd',
                 refine=False, refine_conformers=FINAL_CONFORMERS, top=20, smiles=None):
        validate_conformer_limit(conformers)
        validate_conformer_limit(refine_conformers, MAX_FINAL_CONFORMERS, 'refine_conformers')
        self.shortlist = TopHits(top)
        self.refine, self.refine_conformers = refine, refine_conformers
        self.cache, self.config, self.requested = Path(cache).expanduser().resolve(), config, workers
        self.observed, self.ids = observed, None if ids is None else set(ids)
        self.frequency = frequency
        self.fast = fast
        self.torsions, self.conformers = torsions, conformers
        self.groups = parse_groups(groups)
        self.group_match = parse_group_match(group_match)
        if sources not in ('ccd', 'ccd_chebi'):
            raise ValueError('sources must be ccd or ccd_chebi')
        self.sources = sources
        self.smiles = None if smiles is None else validate_smiles(smiles)
        self.events, self.cancel = queue.Queue(maxsize=128), threading.Event()
        self.thread = None

    def emit(self, event):
        while not self.cancel.is_set():
            try:
                self.events.put(event, timeout=.1)
                return
            except queue.Full:
                pass

    def start(self):
        self.thread = threading.Thread(target=self.run, name='Ligand search', daemon=True)
        self.thread.start()

    def run(self):
        import subprocess
        import tempfile
        from shutil import copyfile
        from contextlib import ExitStack
        import numpy as np

        processes, threads = [], []
        stopped = threading.Event()
        waiting, completed = queue.Queue(), queue.Queue()
        done, failures = 0, 0
        refined = refine_total = refine_failures = 0
        try:
            with tempfile.TemporaryDirectory(prefix='ligand_search_') as temporary, ExitStack() as stack:
                if self.smiles is not None:
                    ordered = [(None, 'LIG', None)]
                else:
                    if library_info(self.cache) is None:
                        raise ValueError('CCD library missing. Click Expand/Update CCD library first.')
                    db = stack.enter_context(closing(sqlite3.connect((self.cache / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)))
                    current_ids = [row[0] for row in db.execute('SELECT id FROM components')]
                    current_id_set = set(current_ids)
                    classification = group_info(self.cache, current_ids)
                    group_db = None
                    if classification is not None and classification['current']:
                        group_db = stack.enter_context(closing(sqlite3.connect((self.cache / 'groups.sqlite').as_uri() + '?mode=ro', uri=True)))
                    if self.groups != ('all',) and (classification is None or not classification['current']):
                        raise ValueError('Chemical-group tags are missing or out of date; click Update chemical groups, '
                                         'or use groups all to search without a chemical filter.')
                    members = group_members(group_db, self.groups, self.group_match) if group_db is not None else None
                    occurrence = frequency_info(self.cache) if self.frequency else None
                    counts = None if occurrence is None else occurrence['counts']
                    if self.frequency and occurrence is None:
                        self.emit(dict(kind='warning', text='PDB occurrence cache unavailable; using volume-only ordering. '
                                       'Expand/Update CCD library to fetch occurrence counts.'))
                    ordered = library_order(db, self.config['volume'], self.config['fraction'], counts)
                    if self.ids is not None:
                        ordered = [row for row in ordered if row[1] in self.ids]
                    if members is not None:
                        ordered = [row for row in ordered if row[1] in members]
                    expanded = {}
                    if self.sources == 'ccd_chebi':
                        expansion = chebi_api().info(self.cache)
                        if expansion is None or int(expansion['count']) == 0:
                            raise ValueError('ChEBI library missing; click Expand/Update ChEBI metabolites first.')
                        if expansion.get('state') != 'complete':
                            self.emit(dict(kind='warning', text='ChEBI expansion is incomplete; searching cached entries only. '
                                           'Use Expand/Update ChEBI metabolites to resume.'))
                        extra_db = stack.enter_context(closing(sqlite3.connect((self.cache / 'chebi.sqlite').as_uri() + '?mode=ro', uri=True)))
                        eligible = set()
                        for ident, payload in extra_db.execute('SELECT id,payload FROM components'):
                            record = json.loads(payload)
                            if self.ids is not None and not (self.ids & {ident, *record.get('chebi_ids', [])}):
                                continue
                            annotation = group_evidence(group_db, ident) if group_db is not None and ident in current_id_set else {}
                            record['chemical_groups'] = sorted(set(record.get('chemical_groups') or []) | set(annotation.get('chemical_groups') or []))
                            record['group_evidence'] += annotation.get('group_evidence', [])
                            if hit_matches_groups(record, self.groups, self.group_match):
                                expanded[ident] = record
                                eligible.add(ident)
                        # Apply the same volume/frequency ordering to both sources,
                        # not a cached-first or CCD-first queue.
                        ordered = [row for row in ordered if row[1] not in eligible]
                        extra_order = library_order(extra_db, self.config['volume'], self.config['fraction'], counts)
                        ordered.extend(row for row in extra_order if row[1] in eligible)
                        ordered.sort(key=lambda row: (-row[0], row[1]))
                total = len(ordered)
                if not total:
                    raise ValueError('No library entries match the selected groups and optional codes. '
                                     'Use groups all, choose another group, or clear the CCD codes filter.')
                api = shape_api()
                count = min(total, max(1, api._shape_available_cpus() - 1) if self.requested == 0 else self.requested,
                            api._shape_available_cpus())
                mode = 'fast' if self.fast else 'standard'
                detail = f'; sparse torsions, up to {self.conformers} conformers/ligand' if self.torsions else ''
                detail += f'; groups ({"AND" if self.group_match == "all" else "OR"}): ' + ', '.join(self.groups)
                self.emit(dict(kind='status', text=f'Searching {total:,} ligand entries with {count} workers ({mode} search{detail})',
                               total=total, done=0))
                if self.smiles is not None:
                    self.emit(dict(kind='status', text='Converting SMILES to 3D and fitting one ligand...'))
                root = Path(temporary)
                copyfile(Path(__file__), root / 'ligand_search.py')
                # Snapshot both files so a script reload cannot change a running job.
                copyfile(Path(__file__).with_name('chimerax_trimmings.py'), root / 'chimerax_trimmings.py')
                np.savez(root / 'blob.npz', **self.config)
                import sys
                bootstrap = ('import sys,runpy,json;sys.path[:]=json.loads(sys.argv[2]);'
                             'n=runpy.run_path(sys.argv[1]+"/ligand_search.py");n["worker_main"](sys.argv[1])')
                environment = os.environ.copy()
                for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                    environment[name] = '1'
                environment.pop('PYTHONHOME', None)
                environment['PYTHONDONTWRITEBYTECODE'] = '1'
                try:
                    from chimerax.core.python_utils import chimerax_python_executable
                    executable = chimerax_python_executable()
                except ImportError:
                    executable = sys.executable
                if not executable or not Path(executable).is_file():
                    raise ValueError('Cannot locate ChimeraX Python for fitting workers.')

                def exchange(process):
                    try:
                        while not stopped.is_set():
                            try:
                                job = waiting.get(timeout=.1)
                            except queue.Empty:
                                continue
                            process.stdin.write(json.dumps(job) + '\n')
                            process.stdin.flush()
                            while not stopped.is_set():
                                line = process.stdout.readline()
                                if not line:
                                    raise RuntimeError('A fitting worker exited unexpectedly.')
                                if line.startswith('LIGAND_RESULT:'):
                                    result = json.loads(line[len('LIGAND_RESULT:'):])
                                    if result['id'] != job['component']['id']:
                                        raise RuntimeError('A worker returned the wrong CCD component.')
                                    if 'result' in result and 'refine_hit' not in job and 'smiles' not in job:
                                        result['result'].update(estimated_volume=job['estimated_volume'],
                                                                search_priority=job['search_priority'],
                                                                pdb_entries=job['pdb_entries'], **job['classification'])
                                    completed.put(result)
                                    break
                    except Exception as error:
                        completed.put(dict(fatal=str(error)))

                try:
                    for _ in range(count):
                        process = subprocess.Popen([executable, '-u', '-c', bootstrap, str(root), json.dumps(sys.path)],
                                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                                   text=True, encoding='utf-8', env=environment,
                                                   **({'creationflags': subprocess.CREATE_NO_WINDOW} if sys.platform == 'win32' else {}))
                        processes.append(process)
                        thread = threading.Thread(target=exchange, args=(process,), daemon=True)
                        thread.start()
                        threads.append(thread)
                    cursor = iter(ordered)
                    def dispatch():
                        try:
                            prior, ident, estimated = next(cursor)
                        except StopIteration:
                            return
                        if self.smiles is not None:
                            waiting.put(dict(component={'id': ident}, smiles=self.smiles, observed=False,
                                             fast=self.fast, torsions=self.torsions, conformers=self.conformers))
                            return
                        record = expanded.get(ident)
                        if record is None:
                            record = json.loads(db.execute('SELECT payload FROM components WHERE id=?', (ident,)).fetchone()[0])
                        waiting.put(dict(component=record, observed=self.observed, fast=self.fast,
                                         torsions=self.torsions, conformers=self.conformers,
                                         estimated_volume=estimated, search_priority=prior,
                                         classification={key: record[key] for key in ('chemical_groups', 'group_evidence', 'classification_status')}
                                         if ident in expanded else group_evidence(group_db, ident) if group_db is not None else dict(
                                             chemical_groups=None, group_evidence=[], classification_status='unavailable'),
                                         pdb_entries=None if counts is None else counts.get(ident, 0)))
                    for _ in range(count):
                        dispatch()
                    while done < total:
                        _cancelled(self.cancel)
                        try:
                            message = completed.get(timeout=.1)
                        except queue.Empty:
                            continue
                        if 'fatal' in message:
                            raise RuntimeError(message['fatal'])
                        if self.smiles is not None and 'error' in message:
                            raise ValueError(message['error'])
                        done += 1
                        failures += 'error' in message
                        if self.refine and 'result' in message:
                            self.shortlist.add(message['result'])
                        self.emit(dict(kind='result', **message, done=done, total=total, failures=failures))
                        dispatch()
                    if self.refine and self.shortlist.rows:
                        # Freeze membership once screening completes. Reuse the
                        # same numerical workers; no library search is repeated.
                        refine_total = len(self.shortlist.rows)
                        self.emit(dict(kind='status', text=f'Final torsion search: 0/{refine_total} hits complete; '
                                       f'up to {self.refine_conformers} conformers per hit'))
                        for hit in self.shortlist.rows:
                            waiting.put(dict(component={'id': hit['id']}, refine_hit=hit,
                                             refine_conformers=self.refine_conformers))
                        while refined < refine_total:
                            _cancelled(self.cancel)
                            try:
                                message = completed.get(timeout=.1)
                            except queue.Empty:
                                continue
                            if 'fatal' in message:
                                raise RuntimeError(message['fatal'])
                            refined += 1
                            refine_failures += 'error' in message
                            self.emit(dict(kind='refined', **message, done=refined, total=refine_total,
                                           failures=refine_failures))
                finally:
                    stopped.set()
                    for process in processes:
                        if process.poll() is None:
                            process.terminate()
                    for process in processes:
                        try:
                            process.wait(timeout=2)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait()
                    for thread in threads:
                        thread.join(timeout=2)
                    for process in processes:
                        process.stdin.close()
                        process.stdout.close()
            self.emit(dict(kind='finished', done=done, total=total, failures=failures,
                           refined=refined, refine_total=refine_total, refine_failures=refine_failures))
        except InterruptedError:
            pass
        except Exception as error:
            self.emit(dict(kind='error', text=str(error)))


class LibraryJob:
    """Build in a separate native process so parsing/downloading cannot stall Qt."""
    def __init__(self, cache, source=None, groups_only=False, chebi=False, install=False):
        self.cache, self.source = Path(cache).expanduser(), source
        self.groups_only = groups_only
        self.chebi = chebi
        self.install = install
        self.cancel, self.events = threading.Event(), queue.Queue()
        self.thread = None

    def start(self):
        self.thread = threading.Thread(target=self.run, name='CCD library', daemon=True)
        self.thread.start()

    def run(self):
        import subprocess
        import sys

        process = None
        try:
            try:
                from chimerax.core.python_utils import chimerax_python_executable
                executable = chimerax_python_executable()
            except ImportError:
                executable = sys.executable
            action = '--build-chebi' if self.chebi else '--update-groups' if self.groups_only else '--build-library'
            command = [executable, '-u', str(Path(__file__).resolve()), action, str(self.cache), '--refresh']
            if self.install:
                command = [executable, '-u', str(Path(__file__).resolve()), '--install-library', str(self.source),
                           '--cache', str(self.cache)]
            elif self.source is not None:
                command.extend(('--ontology' if self.chebi else '--source', str(self.source)))
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       text=True, encoding='utf-8',
                                       **({'creationflags': subprocess.CREATE_NO_WINDOW} if sys.platform == 'win32' else {}))
            lines = queue.Queue()
            def read():
                for line in process.stdout:
                    lines.put(line.rstrip())
                lines.put(None)
            reader = threading.Thread(target=read, daemon=True)
            reader.start()
            tail = []
            while True:
                _cancelled(self.cancel)
                try:
                    line = lines.get(timeout=.1)
                except queue.Empty:
                    continue
                if line is None:
                    break
                tail = (tail + [line])[-5:]
                if line.startswith(('Downloading CCD:', 'Indexing CCD:', 'CCD library ready:',
                                    'Counting PDB occurrences:', 'PDB frequencies ready:', 'Chemical groups:', 'Chemical groups ready:',
                                    'ChEBI library:', 'Ligand library:')):
                    self.events.put(dict(kind='status', text=line))
                elif line.startswith('Warning:'):
                    self.events.put(dict(kind='warning', text=line[len('Warning:'):].strip()))
            if process.wait() != 0:
                action_name = 'Chemical-group update' if self.groups_only else 'ChEBI update' if self.chebi else 'Ligand library setup'
                raise RuntimeError(action_name + ' failed: ' + (tail[-1] if tail else 'no details from the worker.'))
            if self.chebi:
                self.events.put(dict(kind='chebi', info=chebi_api().info(self.cache)))
            elif self.groups_only:
                info = group_info(self.cache)
                if info is None or not info['current']:
                    raise RuntimeError('Chemical-group tags were not updated; retry Update chemical groups.')
                self.events.put(dict(kind='groups', info=info))
            else:
                info = library_info(self.cache)
                if info is None:
                    raise RuntimeError('No usable CCD index was installed. Retry the library setup.')
                extra = chebi_api().info(self.cache)
                self.events.put(dict(kind='library', info=info,
                                     combined=extra is not None and extra.get('state') == 'complete'))
        except InterruptedError:
            pass
        except Exception as error:
            self.events.put(dict(kind='error', text=str(error)))
        finally:
            if process is not None:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                process.stdout.close()


def make_model(session, row, volume, existing=None):
    """Materialize only a user-requested preview or kept hit, on the GUI thread."""
    import numpy as np
    from chimerax.atomic import AtomicStructure, Element
    from chimerax.atomic.colors import element_colors
    from chimerax.geometry import Place

    model = existing
    fresh = model is None or model.deleted
    old_residues = None
    if fresh:
        model = AtomicStructure(session, name=row['id'], auto_style=False, log_info=False)
    else:
        old_residues = model.residues
    residue = None
    try:
        # Deleting a structure's last residue can delete the native structure.
        # Build the replacement first, then discard only the old residues.
        residue = model.new_residue('LIG' if row['id'].startswith('CHEBI:') else row['id'], 'A', 1 if fresh else max(old_residues.numbers) + 1)
        xyz = Place(matrix=np.array(row['pose'])).transform_points(np.array(row['xyz'], dtype=float))
        atoms = []
        for name, element, point in zip(row['names'], row['elements'], xyz):
            atom = model.new_atom(name, Element.get_element(int(element)))
            residue.add_atom(atom)
            atom.coord = point
            atom.draw_mode = atom.STICK_STYLE
            atoms.append(atom)
        for a, b, order in row['bonds']:
            model.new_bond(atoms[a], atoms[b])
        if old_residues is not None:
            old_residues.delete()
            residue.number = 1
        model.atoms.colors = element_colors(np.array(row['elements'], dtype=np.uint8))
        label = ('SMILES' if row.get('library_source') == 'SMILES' else
                 row['id'] if row['id'].startswith('CHEBI:') else f'CCD {row["id"]}')
        model.name = f'{label}: {row["name"]}'
        model.position = volume.scene_position
        if fresh:
            session.models.add([model])
        model.display = True
        return model
    except BaseException:
        if fresh:
            model.delete()
        elif residue is not None and not residue.deleted:
            residue.delete()
        raise


def create_browser(session):
    from Qt.QtCore import QEvent, Qt, QTimer
    from Qt.QtGui import QColor, QFont, QPalette
    from Qt.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                              QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox,
                              QTableWidget, QTableWidgetItem, QVBoxLayout, QAbstractItemView,
                              QMenu, QScrollArea, QToolButton, QWidget, QGridLayout, QSizePolicy,
                              QHeaderView, QProgressBar)
    from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
    from chimerax.core.tools import ToolInstance
    from chimerax.ui import MainToolWindow
    from chimerax.map import Volume

    class NoWheel:
        # Moving the pointer over a control must not change its value.
        def wheelEvent(self, event):
            event.ignore()

    class EntrySpinBox(NoWheel, QSpinBox):
        pass

    class EntryDoubleSpinBox(NoWheel, QDoubleSpinBox):
        pass

    class EntryComboBox(NoWheel, QComboBox):
        pass

    class HitTable(QTableWidget):
        """Let names use spare space without crushing them in a narrow window."""
        def fit_name_column(self):
            if self.columnCount() != 8 or getattr(self, '_fitting_columns', False):
                return
            self._fitting_columns = True
            try:
                other = sum(self.columnWidth(i) for i in range(8) if i != 1)
                self.setColumnWidth(1, max(180, self.viewport().width()-other))
            finally:
                self._fitting_columns = False

        def viewportEvent(self, event):
            result = super().viewportEvent(event)
            if event.type() == QEvent.Type.Resize:
                self.fit_name_column()
            return result

    class LigandBrowser(ToolInstance):
        SESSION_SAVE = False
        SESSION_ENDURING = False

        def __init__(self):
            super().__init__(session, 'Ligand Search')
            self._closed, self.job, self.preview = False, None, None
            self._implementation = IMPLEMENTATION
            self.include_observed = True
            self.current, self.volume, self.watched_data = None, None, None
            self.blob_volume = 0.
            self.buffer, self.failures = TopHits(20), []
            self.result_groups, self.displayed_hits = None, []
            self.result_group_match = 'any'
            self.result_sources = None
            self.seed, self._pick_mode, self._previous_mode = None, None, None
            self.blob_highlight, self.picked_surface = BlobHighlight(), None
            self.stale, self.pending_start, self.started, self.library_completed = False, False, 0., False
            self._has_library = False
            self.group_labels, self.custom_groups = {}, set()
            self.group_counts = {}
            self._combined_group_counts = (None, {})
            self.tool_window = MainToolWindow(self)
            root = QVBoxLayout(self.tool_window.ui_area)
            root.setContentsMargins(10, 10, 10, 8)
            root.setSpacing(5)

            def emphasis(widget):
                font = widget.font()
                font.setWeight(QFont.Weight.DemiBold)
                widget.setFont(font)

            def secondary(widget):
                palette = widget.palette()
                foreground, background = palette.color(QPalette.ColorRole.WindowText), palette.color(QPalette.ColorRole.Window)
                # A subdued active label, not the near-invisible disabled color.
                color = QColor(*(round(.7*getattr(foreground, channel)() + .3*getattr(background, channel)())
                                 for channel in ('red', 'green', 'blue')))
                palette.setColor(QPalette.ColorRole.WindowText, color)
                widget.setPalette(palette)

            toolbar = QHBoxLayout()
            root.addLayout(toolbar)
            self.options_button = QToolButton()
            self.options_button.setText('Search options...')
            self.library_button = QToolButton()
            self.library_button.setText('Library')
            self.library_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
            toolbar.addWidget(self.options_button)
            toolbar.addWidget(self.library_button)
            toolbar.addStretch(1)
            setup = QGridLayout()
            setup.setHorizontalSpacing(8)
            setup.setVerticalSpacing(4)
            setup.setColumnStretch(1, 1)
            root.addLayout(setup)
            self.progress = QProgressBar()
            self.progress.setTextVisible(False)
            self.progress.setRange(0, 1)
            self.progress.setValue(0)
            self.progress.setFixedHeight(5)
            root.addWidget(self.progress)
            def auxiliary_dialog(title):
                dialog = QDialog(self.tool_window.ui_area)
                dialog.setWindowTitle(title)
                layout = QVBoxLayout(dialog)
                close = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
                close.rejected.connect(dialog.close)
                layout.addWidget(close)
                return dialog, layout

            def form_layout(layout):
                form = QFormLayout()
                form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
                form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
                form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
                form.setHorizontalSpacing(12)
                form.setVerticalSpacing(7)
                layout.insertLayout(0, form)
                return form
            self.options_dialog, options_layout = auxiliary_dialog('Ligand Search Options')
            self.options_button.clicked.connect(self.options_dialog.show)
            search_form = form_layout(options_layout)
            self.library_dialog, library_layout = auxiliary_dialog('Ligand Library Status and Cache')
            library_form = form_layout(library_layout)
            self.summary = QLabel('Density fit suggests candidates; it does not establish chemical identity. '
                                  'No receptor clash filtering.')
            self.summary.setWordWrap(True)
            notes = QScrollArea()
            notes.setFrameShape(QScrollArea.Shape.NoFrame)
            notes.setWidgetResizable(True)
            notes.setFixedHeight(72)
            notes.setWidget(self.summary)
            options_layout.insertWidget(1, notes)

            def setup_row(row, title, field, trailing=None):
                label = QLabel(title)
                label.setBuddy(field)
                setup.addWidget(label, row, 0)
                setup.addWidget(field, row, 1, 1, 2 if trailing is None else 1)
                if trailing is not None:
                    setup.addWidget(trailing, row, 2)
                return label

            self.maps = EntryComboBox()
            self.maps.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            self.maps.setMinimumContentsLength(10)
            self.maps.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.maps.currentIndexChanged.connect(self.map_changed)
            self.pick = QPushButton('Pick blob')
            self.pick.setToolTip('Then click the density surface in the graphics window.')
            self.pick.clicked.connect(self.begin_pick)
            setup_row(0, 'Density map', self.maps, self.pick)
            self.resolution = EntryDoubleSpinBox()
            self.resolution.setRange(.5, 20.)
            self.resolution.setDecimals(2)
            self.resolution.setValue(2.6)
            self.resolution.setSuffix(' A')
            self.resolution.setToolTip('Approximate local map resolution, in Angstroms.')
            self.location = QLabel('No blob picked')
            secondary(self.location)
            resolution_row = QHBoxLayout()
            resolution_row.setSpacing(8)
            resolution_row.addWidget(self.resolution)
            resolution_row.addWidget(self.location, 1)
            resolution_label = QLabel('Resolution')
            resolution_label.setBuddy(self.resolution)
            setup.addWidget(resolution_label, 1, 0)
            setup.addLayout(resolution_row, 1, 1, 1, 2)
            self.top = EntrySpinBox()
            self.top.setRange(1, 100)
            self.top.setValue(20)
            self.top.setToolTip('Maximum retained hits. One best pose per CCD code.')
            search_form.addRow('Top-hit buffer', self.top)
            self.workers = EntrySpinBox()
            self.workers.setRange(0, 256)
            self.workers.setSpecialValueText('Automatic')
            search_form.addRow('CPU workers', self.workers)
            self.codes = QLineEdit()
            self.codes.setPlaceholderText('All entries, or ATP,ADP,CHEBI:15422')
            search_form.addRow('Compound IDs', self.codes)
            self.sources = EntryComboBox()
            self.sources.addItem('CCD', 'ccd')
            self.sources.addItem('CCD + ChEBI metabolites', 'ccd_chebi')
            self.sources.addItem('Single ligand (SMILES)', 'smiles')
            self.sources.setCurrentIndex(self.sources.findData(DEFAULT_LIBRARY_SOURCE))
            self.sources.setToolTip('CCD + ChEBI adds metabolites that have usable PubChem3D conformers. '
                                    'Exact CCD matches reuse existing coordinates; duplicates are tested once.')
            setup_row(2, 'Compounds', self.sources)
            self.smiles = QLineEdit()
            self.smiles.setPlaceholderText('Enter a SMILES string')
            self.smiles.setMaxLength(2048)
            self.smiles.setToolTip('Uses ChimeraX\'s built-in online SMILES converter. Only the string is sent; '
                                   'do not enter confidential structures.')
            smiles_label = setup_row(3, 'SMILES (online)', self.smiles)
            smiles_label.hide()
            self.smiles.hide()
            self.library = QLineEdit(str(default_cache()))
            library_form.addRow('Cache folder', self.library)
            self.group_button = QPushButton()
            self.group_button.setAutoDefault(False)
            menu = QMenu(self.group_button)
            families = menu.addMenu('Chemical families')
            roles = menu.addMenu('Biological / usage roles')
            self.group_actions = {}
            for key, label in [('all', 'All compounds')] + [(r[0], r[1]) for r in CHEMICAL_GROUPS] + [('unclassified', 'Unclassified in these groups')]:
                kind = next((row[3] for row in CHEMICAL_GROUPS if row[0] == key), None)
                owner = roles if kind == 'role' else families if kind == 'class' else menu
                action = owner.addAction(label)
                action.setCheckable(True)
                action.setChecked(key == 'all')
                action.toggled.connect(lambda checked, key=key: self.choose_group(key, checked))
                self.group_actions[key] = action
            menu.addSeparator()
            menu.addAction('More structural classes...').triggered.connect(self.choose_structural_classes)
            self.group_button.setMenu(menu)
            self.group_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.group_button.setMinimumWidth(1)
            self.group_button.setStyleSheet('QPushButton { text-align: left; padding: 3px 22px 3px 7px; }')
            self.group_match = EntryComboBox()
            self.group_match.addItem('Any (OR)', 'any')
            self.group_match.addItem('All (AND)', 'all')
            self.group_match.setAccessibleName('Match groups')
            self.group_match.setToolTip('AND requires every selected chemical family/class/role; OR accepts any one.')
            self.group_match.currentIndexChanged.connect(lambda *_: self.set_groups(self.selected_groups()))
            groups_label = setup_row(4, 'Chemical groups', self.group_button, self.group_match)
            def source_changed():
                single = self.sources.currentData() == 'smiles'
                for widget in (smiles_label, self.smiles):
                    widget.setVisible(single)
                for widget in (groups_label, self.group_button, self.group_match):
                    widget.setVisible(not single)
                search_form.setRowVisible(self.codes, not single)
                self.set_busy(self.job is not None)
                self.refresh_group_counts()
            self.sources.currentIndexChanged.connect(source_changed)
            self.group_coverage = QLabel()
            self.group_coverage.setWordWrap(True)
            library_layout.insertWidget(1, self.group_coverage)
            self.library.editingFinished.connect(self.refresh_group_counts)
            self.fast = QCheckBox('Fast screening')
            self.fast.setToolTip('Tests the same CCD entries with fewer orientation starts and refinement samples. '
                                 'Final density scoring is unchanged, but poses and rankings can differ. '
                                 f'Uses {FAST_EXTRA_STARTS} extra starts and {FAST_SAMPLES:,} samples; '
                                 'unchecked uses the standard search settings.')
            search_form.addRow(self.fast)
            self.torsions = QCheckBox('Sample torsions during screening')
            self.torsions.setToolTip('Generate a small deterministic set of additional conformers around flexible '
                                     'single bonds. Slower; not exhaustive or energy-minimized. No receptor clash tests.')
            search_form.addRow(self.torsions)
            self.conformers = EntrySpinBox()
            self.conformers.setRange(2, MAX_CONFORMERS)
            self.conformers.setValue(DEFAULT_CONFORMERS)
            self.conformers.setToolTip('Maximum total per ligand, including its starting conformers. '
                                       'Rigid ligands or rejected/duplicate rotations may yield fewer.')
            self.conformers_label = QLabel(f'Search conformers (max {MAX_CONFORMERS})')
            search_form.addRow(self.conformers_label, self.conformers)
            self.torsions.toggled.connect(self.update_conformer_controls)
            self.refine = QCheckBox('Final torsion search on top hits')
            self.refine.setToolTip('After screening, sample larger torsion changes from each shortlisted fitted '
                                   'conformation and locally fit each variant. Keep the starting pose if none scores better. '
                                   'Independent of sparse torsion sampling during the full-library search.')
            search_form.addRow(self.refine)
            self.refine_conformers = EntrySpinBox()
            self.refine_conformers.setRange(2, MAX_FINAL_CONFORMERS)
            self.refine_conformers.setValue(FINAL_CONFORMERS)
            self.refine_conformers.setToolTip('Final-pass cap per shortlisted hit, including its starting conformation. '
                                              'Fewer variants may survive chemical restrictions and deduplication.')
            self.refine_conformers_label = QLabel(f'Final conformers (max {MAX_FINAL_CONFORMERS})')
            search_form.addRow(self.refine_conformers_label, self.refine_conformers)
            self.refine.toggled.connect(self.update_conformer_controls)
            self.frequency = QCheckBox('Use PDB frequency for search order')
            self.frequency.setChecked(True)
            self.frequency.setToolTip('Common components get a capped priority boost. All components remain eligible; '
                                      'frequency does not affect density scores or the final ranking.')
            search_form.addRow(self.frequency)
            self.start_button, self.stop_button = QPushButton('Start search'), QToolButton()
            self.stop_button.setText('Stop')
            emphasis(self.start_button)
            self.start_button.setAutoDefault(False)
            self.start_button.clicked.connect(self.start_search)
            self.stop_button.clicked.connect(self.stop)
            for button in (self.start_button, self.stop_button):
                toolbar.addWidget(button)
            library_menu = QMenu(self.library_button)
            self.library_button.setMenu(library_menu)
            self.library_button.setMinimumWidth(self.library_button.sizeHint().width() + 8)
            library_menu.addAction('Library status and cache...', self.library_dialog.show)
            library_menu.addSeparator()
            self.update_library_action = library_menu.addAction('Expand/Update CCD library', lambda: self.start_library())
            self.chebi_update = library_menu.addAction('Expand/Update ChEBI metabolites', self.start_chebi)
            self.group_update = library_menu.addAction('Update chemical groups', self.start_groups)
            self.install_action = library_menu.addAction('Install library archive...', self.install_archive)
            self.install_action.setToolTip('Install a prebuilt combined library into an empty cache folder; '
                                           'existing libraries and active downloads are never overwritten.')
            self.chebi_status = QLabel()
            self.chebi_status.setWordWrap(True)
            library_layout.insertWidget(2, self.chebi_status)
            self.library_dialog.resize(430, self.library_dialog.sizeHint().height())
            field_width = self.tool_window.ui_area.fontMetrics().horizontalAdvance('000.00 A') + 33
            for field in (self.resolution, self.top, self.workers,
                          self.conformers, self.refine_conformers):
                field.setFixedWidth(max(field.sizeHint().width(), field_width))
                field.setAlignment(Qt.AlignmentFlag.AlignRight)
                field.setKeyboardTracking(False)
            for widget in (self.codes, self.smiles, self.library, self.sources):
                widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            results = QWidget()
            results_layout = QVBoxLayout(results)
            results_layout.setContentsMargins(0, 0, 0, 0)
            root.addWidget(results, 1)
            self.table = HitTable(0, 8)
            self.table.setMinimumHeight(120)
            self.table.setHorizontalHeaderLabels(('ID', 'Compound', 'Score', 'CC', 'Blob %', 'Atom %', 'Size', 'PDB'))
            descriptions = ('CCD or ChEBI identifier', 'Compound name', 'Density score including blob/atom support',
                            'Density correlation coefficient', 'Percentage of the blob covered by the fitted ligand',
                            'Percentage of ligand atoms supported by density', 'Estimated ligand/blob volume ratio',
                            'Number of PDB entries containing this CCD component')
            for index, description in enumerate(descriptions):
                self.table.horizontalHeaderItem(index).setToolTip(description)
            self.table.setShowGrid(False)
            self.table.setAlternatingRowColors(True)
            self.table.setWordWrap(False)
            self.table.setTextElideMode(Qt.TextElideMode.ElideRight)
            self.table.verticalHeader().setDefaultSectionSize(self.table.fontMetrics().height() + 12)
            self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
            self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
            self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
            self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            self.table.itemSelectionChanged.connect(self.activate_row)
            self.table.horizontalHeader().setStretchLastSection(False)
            for index, width in enumerate((84, 200, 62, 62, 65, 65, 60, 64)):
                self.table.setColumnWidth(index, width)
            results_heading = QHBoxLayout()
            hits_title = QLabel('Top Hits')
            emphasis(hits_title)
            self.hit_count = QLabel('0 retained')
            secondary(self.hit_count)
            results_heading.addWidget(hits_title)
            results_heading.addWidget(self.hit_count)
            results_heading.addStretch(1)
            self.result_scope = QLabel('No search yet.')
            self.result_scope.setWordWrap(True)
            secondary(self.result_scope)
            results_layout.addLayout(results_heading)
            results_layout.addWidget(self.result_scope)
            results_layout.addWidget(self.table, 1)
            self.hit_title = QLabel('No hit selected')
            self.hit_title.setTextFormat(Qt.TextFormat.PlainText)
            self.hit_title.setWordWrap(True)
            emphasis(self.hit_title)
            self.details = QLabel('Choose a result to preview it in the density.')
            self.details.setTextFormat(Qt.TextFormat.PlainText)
            self.details.setWordWrap(True)
            self.details.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.details_scroll = QScrollArea()
            self.details_scroll.setWidgetResizable(True)
            self.details_scroll.setFrameShape(self.details_scroll.Shape.NoFrame)
            self.details_scroll.setMinimumHeight(48)
            self.details_scroll.setMaximumHeight(100)
            details_content = QWidget()
            details_layout = QVBoxLayout(details_content)
            details_layout.setContentsMargins(0, 3, 0, 3)
            details_layout.setSpacing(4)
            details_layout.addWidget(self.hit_title)
            details_layout.addWidget(self.details)
            details_layout.addStretch(1)
            self.details_scroll.setWidget(details_content)
            results_layout.addWidget(self.details_scroll)
            navigation = QGridLayout()
            self.follow_top = QCheckBox('Follow top hit')
            self.follow_top.setChecked(True)
            self.follow_top.setToolTip('Update the density preview when the best hit changes. '
                                       'Choosing another hit pauses following; check this to resume.')
            self.follow_top.toggled.connect(lambda checked: self.refresh_rows())
            results_heading.addWidget(self.follow_top)
            self.hit_buttons = []
            for text, function, row, column, span in (
                    ('Previous', lambda: self.step(-1), 0, 0, 1), ('Next', lambda: self.step(1), 0, 1, 1),
                    ('Focus', self.focus, 0, 2, 1), ('Keep hit', self.keep, 1, 0, 2), ('Export CSV', self.export, 1, 2, 1)):
                button = QPushButton(text)
                button.setAutoDefault(False)
                button.clicked.connect(function)
                self.hit_buttons.append(button)
                navigation.addWidget(button, row, column, 1, span)
            results_layout.addLayout(navigation)
            self.status = QLabel('Ready')
            self.status.setWordWrap(True)
            self.status.setMaximumHeight(3*self.status.fontMetrics().lineSpacing())
            root.addWidget(self.status)
            # Long filenames/chemical names must scroll or wrap, not impose a
            # window-wide minimum. Keep the result table's minimum height intact.
            for widget in (self.summary, self.details, self.hit_title, self.group_coverage,
                           self.chebi_status, self.result_scope, self.status, self.location):
                policy = widget.sizePolicy()
                policy.setHorizontalPolicy(QSizePolicy.Policy.Ignored)
                widget.setSizePolicy(policy)
                widget.setMinimumWidth(0)
            self.timer = QTimer(self.tool_window.ui_area)
            self.timer.setInterval(200)
            self.timer.timeout.connect(self.poll)
            self.handlers = [session.triggers.add_handler(ADD_MODELS, self.models_changed),
                             session.triggers.add_handler(REMOVE_MODELS, self.models_changed)]
            self.refresh_group_counts()
            self.refresh_maps()
            self.set_busy(False)
            self.tool_window.manage('side')

        def set_busy(self, busy):
            single = self.sources.currentData() == 'smiles'
            for widget in (self.update_library_action, self.install_action, self.resolution,
                           self.top, self.workers, self.fast, self.torsions, self.refine, self.frequency,
                           self.codes, self.smiles, self.library, self.maps, self.pick, self.sources):
                widget.setEnabled(not busy)
            for widget in (self.frequency, self.workers, self.top):
                widget.setEnabled(not busy and not single)
            self.group_button.setEnabled(not busy)
            self.group_match.setEnabled(not busy)
            self.group_update.setEnabled(not busy and self._has_library)
            self.chebi_update.setEnabled(not busy and self._has_library)
            self.update_conformer_controls()
            self.stop_button.setEnabled(busy)
            self.start_button.setText(('Searching...' if isinstance(self.job, SearchJob) else 'Updating...')
                                      if busy else 'Start search')
            self.start_button.setEnabled(not busy and self.has_picked_blob())
            if busy:
                self.progress.setRange(0, 0)
            elif self.progress.maximum() == 0:
                self.progress.setRange(0, 1)
                self.progress.setValue(0)

        def update_conformer_controls(self, *args):
            enabled = self.torsions.isEnabled() and self.torsions.isChecked()
            for widget in (self.conformers_label, self.conformers):
                widget.setEnabled(enabled)
            enabled = self.refine.isEnabled() and self.refine.isChecked()
            for widget in (self.refine_conformers_label, self.refine_conformers):
                widget.setEnabled(enabled)

        def selected_groups(self):
            keys = tuple(key for key, action in self.group_actions.items() if key != 'all' and action.isChecked())
            return keys + tuple(sorted(self.custom_groups)) or ('all',)

        def set_groups(self, keys):
            keys = parse_groups(keys)
            self.custom_groups = {key for key in keys if key.startswith('chemont:')}
            for key, action in self.group_actions.items():
                action.blockSignals(True)
                action.setChecked(key in keys)
                action.blockSignals(False)
            labels = dict((row[0], row[1]) for row in CHEMICAL_GROUPS)
            labels.update(all='All compounds', unclassified='Unclassified')
            labels.update(self.group_labels)
            summary = ', '.join(labels.get(key, key) for key in keys)
            self.group_button.setText(summary if len(keys) <= 2 else f'{len(keys)} groups selected')
            mode = 'All selected groups must match' if self.group_match.currentData() == 'all' else 'Any selected group matches'
            self.group_button.setToolTip(summary + '. ' + mode + '; compounds are not duplicated.')
            self.refresh_rows()

        def choose_group(self, key, checked):
            if key == 'all' and checked:
                self.set_groups(('all',))
            else:
                keys = [k for k in self.selected_groups() if k != 'all']
                self.set_groups(keys or ('all',))

        def choose_structural_classes(self):
            from Qt.QtWidgets import QDialog, QDialogButtonBox, QListWidget, QListWidgetItem

            if not self.group_labels:
                self.report_error('Choose Library > Update chemical groups first to download the structural classifications.')
                return
            dialog = QDialog(self.tool_window.ui_area)
            dialog.setWindowTitle('Structural chemical classes')
            layout = QVBoxLayout(dialog)
            layout.addWidget(QLabel('Select chemical families; Match groups controls OR versus AND.'))
            search = QLineEdit()
            search.setPlaceholderText('Find a family, e.g. pyrimidines, carbonates or peptides')
            layout.addWidget(search)
            listing = QListWidget()
            counts = self.group_counts
            for key, label in sorted(self.group_labels.items(), key=lambda pair: pair[1].casefold()):
                item = QListWidgetItem(f'{label} ({counts.get(key, 0):,})  [{key}]')
                item.setData(Qt.ItemDataRole.UserRole, key)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(Qt.CheckState.Checked if key in self.custom_groups else Qt.CheckState.Unchecked)
                listing.addItem(item)
            def filter_rows(text):
                for i in range(listing.count()):
                    item = listing.item(i)
                    item.setHidden(text.casefold() not in item.text().casefold())
            search.textChanged.connect(filter_rows)
            layout.addWidget(listing)
            buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
            buttons.accepted.connect(dialog.accept)
            buttons.rejected.connect(dialog.reject)
            layout.addWidget(buttons)
            dialog.resize(650, 550)
            if dialog.exec() == QDialog.DialogCode.Accepted:
                fixed = [key for key in self.selected_groups() if key != 'all' and not key.startswith('chemont:')]
                chosen = [listing.item(i).data(Qt.ItemDataRole.UserRole) for i in range(listing.count())
                          if listing.item(i).checkState() == Qt.CheckState.Checked]
                self.set_groups(fixed + chosen or ('all',))

        def refresh_group_counts(self):
            cache = Path(self.library.text()).expanduser().resolve()
            self.library.setToolTip(str(cache))
            try:
                library, info = library_info(cache), group_info(cache)
            except (ValueError, OSError, sqlite3.Error):
                library, info = None, None
            self._has_library = library is not None
            current = info is not None and info['current']
            self.group_labels = info['labels'] if current else {}
            counts = info['counts'] if current else {}
            total = int(library.get('count', 0)) if library else 0
            try:
                expansion = chebi_api().info(cache)
                if current and expansion is not None and self.sources.currentData() == 'ccd_chebi':
                    versions = tuple((entry.st_ino, entry.st_mtime_ns, entry.st_size)
                                     for entry in ((cache / name).stat() for name in ('groups.sqlite', 'chebi.sqlite')))
                    key = (cache, versions)
                    if self._combined_group_counts[0] != key:
                        self._combined_group_counts = (key, chebi_api().group_counts(cache))
                    counts = self._combined_group_counts[1]
                    total = counts['all']
            except (OSError, sqlite3.Error):
                expansion = None
            self.group_counts = counts
            labels = dict((r[0], r[1]) for r in CHEMICAL_GROUPS)
            labels.update(all='All compounds', unclassified='Unclassified in these groups')
            for key, action in self.group_actions.items():
                count = f'{total:,}' if key == 'all' else f'{counts.get(key, 0):,}' if current else '?'
                action.setText(f'{labels[key]} ({count})')
            if current:
                self.group_coverage.setText(f'{info["structural"]}/{info["total"]} CCD entries have class-level structural assignments; '
                                            f'{info["structural_errors"]} need structure review. '
                                            'Roles are separate curated annotations.')
            else:
                self.group_coverage.setText('Chemical-group tags are ' + ('out of date' if info else 'not downloaded') +
                                            '. Update chemical groups to filter; All compounds remains available.')
            self.group_update.setEnabled(self.job is None and self._has_library)
            self.chebi_status.setText('ChEBI metabolites: not downloaded.' if expansion is None else
                                      f'ChEBI metabolites: {expansion["count"]} cached entries; {expansion.get("state", "partial")}. '
                                      'Update resumes cached downloads. Group counts describe the selected library.')
            self.chebi_update.setEnabled(self.job is None and self._has_library)
            self.set_groups(self.selected_groups())

        def refresh_maps(self):
            selected = self.maps.currentData()
            self.maps.blockSignals(True)
            self.maps.clear()
            for model in session.models.list(type=Volume):
                self.maps.addItem(f'#{model.id_string} {model.name}', model)
            index = self.maps.findData(selected)
            if index >= 0:
                self.maps.setCurrentIndex(index)
            self.maps.blockSignals(False)
            self.maps.setToolTip(self.maps.currentText())
            if selected is not self.maps.currentData():
                self.map_changed()
            self.start_button.setEnabled(self.job is None and self.has_picked_blob())

        def has_picked_blob(self):
            model = self.maps.currentData()
            return (self.seed is not None and model is not None and not model.deleted and
                    self.picked_surface is not None and not self.picked_surface.deleted)

        def map_changed(self, *args):
            self.blob_highlight.clear()
            self.picked_surface = None
            self.seed = None
            self.unwatch()
            self.location.setText('No blob picked')
            self.start_button.setEnabled(False)

        def models_changed(self, name, models):
            if self.volume in models and name == REMOVE_MODELS:
                self.invalidate('The fitting map was closed; results are stale. Choose a map and restart.')
            self.refresh_maps()

        def invalidate(self, text):
            self.stale = True
            self.pending_start = False
            if isinstance(self.job, SearchJob):
                self.stop()
            self.remove_preview()
            self.hit_title.setText('Results need a new search')
            self.update_hit_actions()
            self.blob_highlight.clear()
            self.picked_surface = None
            self.seed = None
            self.location.setText('Pick the blob again')
            self.start_button.setEnabled(False)
            self.status.setText(text)

        def data_changed(self, reason):
            if reason in ('values changed', 'coordinates changed'):
                self.invalidate('Map data changed; results are stale. Restart the search.')

        def pick_blob(self, surface, triangle):
            import numpy as np

            owner = next((m for m in surface.drawing_lineage if isinstance(m, Volume)), None) if surface else None
            if owner is None or surface.level <= 0:
                self.status.setText('Click a positive-density map surface, not an atomic model.')
                return
            self.maps.setCurrentIndex(self.maps.findData(owner))
            point = np.mean(surface.vertices[surface.triangles[triangle]], axis=0)
            self.seed = owner.scene_position.inverse() * (surface.scene_position * point)
            self.picked_surface = surface
            self.blob_highlight.select(surface, triangle=triangle)
            self.unwatch()
            self.watched_data = owner.data
            self.watched_data.add_change_callback(self.data_changed)
            self.location.setText(f'Blob contour: {surface.level:.4g}')
            self.start_button.setEnabled(self.job is None)
            self.status.setText('Blob selected (green). Click Start search.')

        def show_blob_selection(self, model, config):
            surfaces = list(model.surfaces)
            if not surfaces:
                return False
            surface = self.picked_surface
            if surface not in surfaces:
                positive = [s for s in surfaces if s.level > 0]
                if not positive:
                    return False
                surface = min(positive, key=lambda s: abs(s.level-config.get('initial_level', config['level'])))
            changed = surface.level != config['level']
            self.blob_highlight.clear()
            if changed:
                levels = [s.level for s in surfaces]
                levels[surfaces.index(surface)] = config['level']
                model.set_parameters(surface_levels=levels)
            point = config.get('anchor', self.seed)
            if point is None:
                raise ValueError('Pick a blob before starting a search.')
            local = surface.scene_position.inverse() * (model.scene_position * point)
            self.blob_highlight.select(surface, point=local, max_level=config.get('anchor_density'), defer=changed)
            model.update_drawings()
            return changed

        def begin_pick(self):
            from chimerax.mouse_modes import MouseMode
            from chimerax.pick_blobs.pickblobs import _picked_triangle

            browser = self
            class Pick(MouseMode):
                name = 'ligand search blob'
                def mouse_down(mode, event):
                    try:
                        surface, triangle = _picked_triangle(session.main_view.picked_object(*event.position()))
                        browser.pick_blob(surface, triangle)
                    except Exception as error:
                        browser.report_error(error)
                    finally:
                        browser.restore_pick()
            self.restore_pick()
            self._previous_mode = session.ui.mouse_modes.mode('left', [], exact=True)
            self._pick_mode = Pick(session)
            session.ui.mouse_modes.bind_mouse_mode('left', [], self._pick_mode)
            self.status.setText('Click the desired density blob once; the previous mouse mode will be restored.')

        def restore_pick(self):
            if self._pick_mode is not None:
                manager = session.ui.mouse_modes
                if manager.mode('left', [], exact=True) is self._pick_mode:
                    manager.bind_mouse_mode('left', [], self._previous_mode)
                self._pick_mode = self._previous_mode = None

        def report_error(self, error):
            self.status.setText(str(error))
            self.status.setToolTip(str(error))
            self.progress.setRange(0, 1)
            self.progress.setValue(0)
            session.logger.warning('Ligand search: ' + str(error))

        def start_search(self):
            if self.job is not None:
                return
            try:
                model = self.maps.currentData()
                if model is None or model.deleted:
                    raise ValueError('Open a density map, then choose it in Density map.')
                if not self.has_picked_blob():
                    raise ValueError('Click Pick blob, then click a positive-density surface before starting a search.')
                sources = self.sources.currentData()
                smiles = validate_smiles(self.smiles.text()) if sources == 'smiles' else None
                cache = Path(self.library.text()).expanduser().resolve()
                groups = ('all',) if smiles else self.selected_groups()
                if smiles is None:
                    if library_info(cache) is None:
                        self.pending_start = True
                        self.start_library()
                        return
                    tags = group_info(cache) if groups != ('all',) else None
                    if groups != ('all',) and (tags is None or not tags['current']):
                        self.pending_start = True
                        self.start_groups()
                        return
                    if sources == 'ccd_chebi':
                        expansion = chebi_api().info(cache)
                        if expansion is None or int(expansion['count']) == 0:
                            self.pending_start = True
                            self.start_chebi()
                            return
                self.restore_pick()
                requested_level = float(self.picked_surface.level)
                config = prepare_blob(model.data, self.seed, self.resolution.value(), requested_level)
                displayed_changed = self.show_blob_selection(model, config)
                self.location.setText(f'Blob contour: {config["level"]:.4g}')
                adjustment = ''
                if config['level'] > requested_level:
                    adjustment = (f'Blob contour automatically raised from {requested_level:.6g} to '
                                  f'{config["level"]:.6g} to separate the selected blob from surrounding density.')
                    if displayed_changed:
                        adjustment += ' The displayed map contour was updated too.'
                    session.logger.info('Ligand search: ' + adjustment)
                    session.logger.status('Ligand search: ' + adjustment, blank_after=0)
                self.volume = model
                self.blob_volume = config['volume']
                self.unwatch()
                self.watched_data = model.data
                self.watched_data.add_change_callback(self.data_changed)
                self.remove_preview()
                top = 1 if smiles else self.top.value()
                self.buffer, self.current, self.failures = TopHits(top), None, []
                self.hit_title.setText('Searching for candidates')
                self.result_groups = groups
                self.result_group_match = self.group_match.currentData()
                self.result_sources = self.sources.currentData()
                self.stale = False
                self.follow_top.setChecked(True)
                self.details.setText('Searching; density scores are provisional. No clash filtering.')
                self.refresh_rows()
                mode = 'Fast search' if self.fast.isChecked() else 'Standard search'
                if self.torsions.isChecked():
                    mode += f' with sparse torsions (up to {self.conformers.value()} conformers/ligand)'
                if self.refine.isChecked():
                    mode += f'; final torsion search ({self.refine_conformers.value()} conformers/hit)'
                candidates = ('Fitting one SMILES ligand; no library search. ' if smiles else
                              f'Keeping top {top} hits. Volume and optional PDB frequency prioritize the queue. ')
                self.summary.setText(f'{mode}; blob volume {config["volume"]:.3g} A^3. {candidates}'
                                     'Ranking uses density fit. No clash tests. '
                                     + adjustment)
                if adjustment:
                    self.status.setText(adjustment)
                codes = [v.strip().upper() for v in self.codes.text().replace(' ', ',').split(',') if v.strip()] or None
                self.job = SearchJob(cache, config, 1 if smiles else self.workers.value(),
                                     observed=not smiles and self.include_observed, ids=None if smiles else codes,
                                     frequency=not smiles and self.frequency.isChecked(), fast=self.fast.isChecked(),
                                     torsions=self.torsions.isChecked(), conformers=self.conformers.value(), groups=groups,
                                     group_match=self.group_match.currentData(), sources='ccd' if smiles else sources,
                                     refine=self.refine.isChecked(), refine_conformers=self.refine_conformers.value(),
                                     top=top, smiles=smiles)
                self.begin_job()
            except Exception as error:
                self.report_error(error)

        def start_library(self, source=None):
            if self.job is not None:
                return
            self.job = LibraryJob(Path(self.library.text()).expanduser().resolve(), source)
            self.begin_job()

        def start_groups(self):
            if self.job is not None:
                return
            self.job = LibraryJob(Path(self.library.text()).expanduser().resolve(), groups_only=True)
            self.begin_job()

        def start_chebi(self):
            if self.job is not None:
                return
            self.job = LibraryJob(Path(self.library.text()).expanduser().resolve(), chebi=True)
            self.begin_job()

        def install_archive(self):
            if self.job is not None:
                return
            path, _ = QFileDialog.getOpenFileName(self.tool_window.ui_area, 'Install combined ligand library', '',
                                                  'Ligand library (*.tar.gz)')
            if not path:
                return
            self.job = LibraryJob(Path(self.library.text()).expanduser().resolve(), source=path, install=True)
            self.begin_job()

        def begin_job(self):
            self.started = time.monotonic()
            self.library_completed = False
            self.set_busy(True)
            self.job.start()
            self.timer.start()

        def stop(self):
            self.restore_pick()
            self.pending_start = False
            if self.job is not None:
                self.job.cancel.set()
                self.status.setText('Stopping; completed top hits are retained.')

        def poll(self):
            job = self.job
            if job is None or self._closed:
                return
            if isinstance(job, SearchJob) and self.volume is not None and not self.volume.deleted and \
                    self.volume.data is not self.watched_data:
                self.invalidate('Map data replaced; results are stale. Restart the search.')
            changed, last_status, finish_error, library_ready = False, None, None, None
            for _ in range(128):
                try:
                    event = job.events.get_nowait()
                except queue.Empty:
                    break
                kind = event['kind']
                if kind == 'result':
                    self.progress.setRange(0, max(1, event['total']))
                    self.progress.setValue(event['done'])
                    if 'result' in event:
                        self.buffer.add(event['result'])
                        changed = True
                    elif len(self.failures) < 10:
                        self.failures.append(f'{event["id"]}: {event["error"]}')
                    last_status = (f'{event["done"]:,}/{event["total"]:,} checked; '
                                   f'{event["failures"]:,} failed; {time.monotonic()-self.started:.0f}s elapsed')
                elif kind == 'refined':
                    self.progress.setRange(0, max(1, event['total']))
                    self.progress.setValue(event['done'])
                    if 'result' in event:
                        self.buffer.add(event['result'], replace_equal=True)
                        changed = True
                    elif len(self.failures) < 10:
                        self.failures.append(f'{event["id"]}: final torsion search failed; original retained: {event["error"]}')
                    last_status = (f'Final torsion search: {event["done"]}/{event["total"]} hits complete; '
                                   f'{event["failures"]} failed (originals retained); {time.monotonic()-self.started:.0f}s elapsed')
                elif kind == 'status':
                    last_status = event['text']
                    if 'total' in event:
                        self.progress.setRange(0, max(1, event['total']))
                        self.progress.setValue(event.get('done', 0))
                    elif last_status.startswith('Final torsion search:'):
                        self.progress.setRange(0, 0)
                elif kind == 'warning':
                    session.logger.warning('Ligand search: ' + event['text'])
                    self.summary.setText(self.summary.text() + ' ' + event['text'])
                elif kind == 'library':
                    self.progress.setRange(0, 1)
                    self.progress.setValue(1)
                    library_ready = event['info']
                    self.library_completed = True
                    self.refresh_group_counts()
                    if event.get('combined'):
                        last_status = 'Combined CCD + ChEBI library ready'
                    else:
                        last_status = f'CCD library ready: {library_ready["count"]} usable components'
                elif kind == 'groups':
                    self.progress.setRange(0, 1)
                    self.progress.setValue(1)
                    self.library_completed = True
                    self.refresh_group_counts()
                    last_status = f'Chemical groups ready: {event["info"]["tagged"]}/{event["info"]["total"]} entries tagged'
                elif kind == 'chebi':
                    self.progress.setRange(0, 1)
                    self.progress.setValue(1)
                    self.library_completed = True
                    self.refresh_group_counts()
                    last_status = f'ChEBI library ready: {event["info"]["count"]} cached entries'
                    if int(event['info']['count']) == 0:
                        self.pending_start = False
                        finish_error = 'No usable ChEBI 3D conformers found; use CCD or retry the ChEBI update.'
                elif kind == 'error':
                    finish_error = event['text']
                    self.pending_start = False
                elif kind == 'finished':
                    self.progress.setRange(0, 1)
                    self.progress.setValue(1)
                    last_status = (f'Finished: {event["done"]:,}/{event["total"]:,} checked; '
                                   f'{event["failures"]:,} failed; {time.monotonic()-self.started:.0f}s')
                    if event.get('refine_total'):
                        last_status += f'; final torsion search {event["refined"]}/{event["refine_total"]} '
                        last_status += f'({event["refine_failures"]} failed; originals retained)'
            if changed:
                self.refresh_rows()
            if last_status and not self.stale:
                self.status.setText(last_status)
                self.status.setToolTip(last_status)
                session.logger.status('Ligand search: ' + last_status, blank_after=0)
            if finish_error:
                self.report_error(finish_error)
            if not job.thread.is_alive() and job.events.empty():
                self.job = None
                self.timer.stop()
                self.set_busy(False)
                if job.cancel.is_set() and not self.stale:
                    self.progress.setRange(0, 1)
                    self.progress.setValue(0)
                    self.status.setText(f'Stopped; {len(self.buffer.rows)} top hits retained.')
                    if isinstance(job, LibraryJob) and job.chebi:
                        self.refresh_group_counts()
                if self.failures:
                    session.logger.warning('Ligand search skipped some candidates. Examples:\n' + '\n'.join(self.failures))
                if self.library_completed and self.pending_start:
                    self.pending_start = False
                    self.start_search()

        def visible_hit(self, hit, groups, match):
            source = self.sources.currentData()
            custom = hit.get('library_source') == 'SMILES'
            if custom or source == 'smiles':
                return custom and source == 'smiles'
            return ((source != 'ccd' or not hit['id'].startswith('CHEBI:')) and
                    hit_matches_groups(hit, groups, match))

        def filtered_hits(self):
            groups, match = self.selected_groups(), self.group_match.currentData()
            return [hit for hit in self.buffer.rows if self.visible_hit(hit, groups, match)]

        def refresh_rows(self):
            groups = self.selected_groups()
            match = self.group_match.currentData()
            labels = dict((row[0], row[1]) for row in CHEMICAL_GROUPS)
            labels.update(all='All compounds', unclassified='Unclassified')
            labels.update(self.group_labels)
            def describe(keys, mode):
                joiner = ' AND ' if mode == 'all' else ' OR '
                if len(keys) > 2:
                    return f'{len(keys)} groups ({"AND" if mode == "all" else "OR"})'
                return joiner.join(labels.get(key, key) for key in keys)
            if self.sources.currentData() == 'smiles':
                self.result_scope.setText('Single ligand (SMILES); no library search.')
            elif self.result_groups is not None and set(groups) == set(self.result_groups) and match == self.result_group_match:
                self.result_scope.setText('Search groups: ' + describe(self.result_groups, self.result_group_match))
            elif self.result_groups is not None or self.buffer.rows:
                source = describe(self.result_groups, self.result_group_match) if self.result_groups is not None else 'the existing shortlist'
                self.result_scope.setText(f'Showing {describe(groups, match)} from {source}. '
                                          'Start search to screen the selected groups fully.')
            else:
                self.result_scope.setText('No search yet.')
            if self.result_sources is not None and self.result_sources != self.sources.currentData():
                self.result_scope.setText(self.result_scope.text() + ' Library changed; start a new search for complete coverage.')
            self.displayed_hits = self.filtered_hits()
            self.hit_count.setText(f'{len(self.displayed_hits)} retained')
            if self.current is not None and not self.visible_hit(self.current, groups, match):
                # A manually pinned preview must not masquerade as a filtered hit.
                self.current = None
                self.remove_preview()
                self.hit_title.setText('No matching hit selected')
                self.details.setText('No matching retained hit selected. Start search to screen these groups fully.')
                self.details.setToolTip('')
            selected = self.current['id'] if self.current else None
            if (self.follow_top.isChecked() or self.current is None) and self.displayed_hits:
                selected = self.displayed_hits[0]['id']
            self.table.blockSignals(True)
            self.table.setRowCount(len(self.displayed_hits))
            index = -1
            for row_index, hit in enumerate(self.displayed_hits):
                values = (hit['id'], hit['name'], f'{hit.get("score", hit["cc"]):.3f}', f'{hit["cc"]:.3f}', f'{100*hit["coverage"]:.0f}%',
                          f'{100*hit["supported"]:.0f}%',
                          f'{hit["estimated_volume"] / self.blob_volume:.2f}'
                          if self.blob_volume > 0 and 'estimated_volume' in hit else '',
                          f'{hit["pdb_entries"]:,}' if hit.get('pdb_entries') is not None else '')
                for column, value in enumerate(values):
                    item = QTableWidgetItem(value)
                    item.setToolTip(value)
                    item.setTextAlignment((Qt.AlignmentFlag.AlignLeft if column < 2 else Qt.AlignmentFlag.AlignRight) |
                                          Qt.AlignmentFlag.AlignVCenter)
                    if column == 0:
                        font = item.font()
                        font.setWeight(QFont.Weight.DemiBold)
                        item.setFont(font)
                    self.table.setItem(row_index, column, item)
                if hit['id'] == selected:
                    index = row_index
            self.table.clearSelection()
            if index >= 0:
                self.table.selectRow(index)
            self.table.blockSignals(False)
            self.table.fit_name_column()
            # Do not rebuild the preview for worse incoming hits. Object identity
            # also catches an improved pose for the same CCD component.
            if index >= 0 and (self.current is not self.displayed_hits[index] or
                               self.preview is None or self.preview.deleted):
                self.activate_row(manual=False)
            elif self.current is not None and index < 0:
                self.show_details(pinned=True)
            self.update_hit_actions()

        def update_hit_actions(self):
            available = self.current is not None and not self.stale and self.preview is not None and not self.preview.deleted
            for button in self.hit_buttons:
                if button.text() in ('Previous', 'Next'):
                    button.setEnabled(not self.stale and len(self.displayed_hits) > 1)
                elif button.text() == 'Export CSV':
                    button.setEnabled(bool(self.displayed_hits))
                else:
                    button.setEnabled(available)

        def activate_row(self, manual=True):
            if self.stale:
                return
            rows = self.table.selectionModel().selectedRows()
            if not rows:
                return
            if manual:
                self.follow_top.blockSignals(True)
                self.follow_top.setChecked(False)
                self.follow_top.blockSignals(False)
            self.current = self.displayed_hits[rows[0].row()]
            try:
                if self.volume is None or self.volume.deleted:
                    raise ValueError('The fitting map is no longer available.')
                self.preview = make_model(session, self.current, self.volume, self.preview)
                self.preview.SESSION_SAVE = False
                self.show_details()
                self.update_hit_actions()
            except Exception as error:
                self.report_error(error)

        def show_details(self, pinned=False):
            if self.current is not None:
                hit = self.current
                keys = hit.get('chemical_groups')
                labels = {r[0]: r[1] for r in CHEMICAL_GROUPS}
                labels.update(self.group_labels)
                family = 'Tags unavailable' if keys is None else ', '.join(
                    labels.get(key, key) for key in keys if not key.startswith('chemont:')) or 'No preset family'
                if hit.get('structural_family'):
                    family += '; structural class: ' + hit['structural_family']
                count = hit.get('conformers_tested', 1)
                trials = f'{count} conformer{"s" if count != 1 else ""} tested'
                if hit.get('conformers_skipped'):
                    trials += f'; {hit["conformers_skipped"]} skipped'
                if hit.get('final_refined'):
                    trials += (f'\nFinal search: {hit["refine_conformers_tested"]} conformers; '
                               f'score {hit["pre_refine_score"]:.3f} to {hit.get("score", hit["cc"]):.3f}; '
                               + ('improved. ' if hit['refinement_improved'] else 'starting pose retained. '))
                self.hit_title.setText(f'{hit["id"]}  {hit["name"]}')
                self.details.setText(f'{hit["formula"]} | {hit["source"]} coordinates | {trials}\n'
                                     f'Groups: {family}\n'
                                     f'CC {hit["cc"]:.3f} | Blob {100*hit["coverage"]:.0f}% | '
                                     f'Atom support {100*hit["supported"]:.0f}% | Outside region {100*hit["outside"]:.0f}%' +
                                     ('\nPreview pinned outside the current top-hit buffer.' if pinned else ''))
                self.details.setToolTip('\n'.join(f'{labels.get(e["group"], e["group"])}: {e["source"]}, {e["id"]}'
                                                   for e in hit.get('group_evidence', [])))

        def step(self, direction):
            if self.displayed_hits and not self.stale:
                index = next((i for i, row in enumerate(self.displayed_hits)
                              if self.current is not None and row['id'] == self.current['id']), -1)
                self.follow_top.blockSignals(True)
                self.follow_top.setChecked(False)
                self.follow_top.blockSignals(False)
                self.table.selectRow((index + direction) % len(self.displayed_hits))

        def focus(self):
            if self.preview is not None and not self.preview.deleted:
                from chimerax.core.commands import run
                run(session, f'view #{self.preview.id_string}')

        def keep(self):
            if self.current is None or self.stale or self.volume is None or self.volume.deleted:
                return
            from chimerax.core.undo import UndoState
            try:
                kept = make_model(session, self.current, self.volume)
                undo = UndoState('keep ligand search hit')
                presence = shape_api()._ShapeCopyPresence(session, kept)
                undo.add(presence, 'restore', False, True, option='M')
                session.undo.register(undo)
                session.logger.info(f'Kept {self.current["id"]} as #{kept.id_string}; target models unchanged.')
            except Exception as error:
                self.report_error(error)

        def export(self):
            import csv
            path, _ = QFileDialog.getSaveFileName(self.tool_window.ui_area, 'Export ligand hits', 'ligand_hits.csv', 'CSV (*.csv)')
            if not path:
                return
            try:
                with open(path, 'w', newline='') as stream:
                    writer = csv.writer(stream)
                    writer.writerow(('CCD', 'name', 'formula', 'score', 'density_cc', 'blob_coverage', 'atoms_supported',
                                     'outside_region', 'conformer', 'pdb_entries', 'fast_search',
                                     'torsion_sampling', 'conformers_tested', 'conformers_skipped', 'torsion_changes',
                                     'chemical_groups', 'structural_family', 'group_evidence', 'classification_status', 'stale',
                                     'library_source', 'chebi_ids', 'pubchem_cid', 'inchikey', 'formal_charge',
                                     'final_refined', 'refinement_improved', 'pre_refine_score', 'pre_refine_cc',
                                     'refine_conformers_tested', 'refine_conformers_skipped', 'refine_conformer_limit',
                                     'final_torsion_changes', 'smiles'))
                    for row in self.filtered_hits():
                        writer.writerow((row['id'], row['name'], row['formula'], f'{row.get("score", row["cc"]):.3g}', f'{row["cc"]:.3g}',
                                         f'{row["coverage"]:.3g}', f'{row["supported"]:.3g}',
                                         f'{row["outside"]:.3g}', row['source'], row.get('pdb_entries'),
                                         row.get('fast_search', False), row.get('torsion_sampling', False),
                                         row.get('conformers_tested', 1), row.get('conformers_skipped', 0),
                                         '; '.join(f'{change["bond"]}:{change["degrees"]:+g}'
                                                   for change in row.get('torsion_changes', [])),
                                         '; '.join(row.get('chemical_groups') or ()), row.get('structural_family', ''),
                                         json.dumps(row.get('group_evidence', [])),
                                         row.get('classification_status', 'unavailable'), self.stale,
                                         row.get('library_source', 'CCD'), '; '.join(row.get('chebi_ids', [])),
                                         row.get('pubchem_cid'), row.get('inchikey', ''), row.get('formal_charge'),
                                         row.get('final_refined', False), row.get('refinement_improved', False),
                                         f'{row["pre_refine_score"]:.3g}' if 'pre_refine_score' in row else '',
                                         f'{row["pre_refine_cc"]:.3g}' if 'pre_refine_cc' in row else '',
                                         row.get('refine_conformers_tested', ''), row.get('refine_conformers_skipped', ''),
                                         row.get('refine_conformer_limit', ''), json.dumps(row.get('final_torsion_changes', [])),
                                         row.get('smiles', '')))
            except OSError as error:
                self.report_error(error)

        def remove_preview(self):
            if self.preview is not None and not self.preview.deleted:
                session.models.close([self.preview])
            self.preview = None

        def unwatch(self):
            if self.watched_data is not None and self.data_changed in self.watched_data.change_callbacks:
                self.watched_data.remove_change_callback(self.data_changed)
            self.watched_data = None

        def delete(self):
            if self._closed:
                return
            self._closed = True
            self.stop()
            self.options_dialog.close()
            self.library_dialog.close()
            self.timer.stop()
            self.unwatch()
            self.blob_highlight.clear()
            for handler in self.handlers:
                handler.remove()
            self.remove_preview()
            if getattr(session, '_trimmings_ligand_browser', None) is self:
                session._trimmings_ligand_browser = None
            super().delete()

    return LigandBrowser()


def ligand_search(session, map_model=None, resolution=2.6, top=20,
                  library=None, workers=0, observed=True, codes=None, frequency=True, fast=False,
                  torsions=False, conformers=DEFAULT_CONFORMERS, groups='all', group_match='any', sources=None,
                  refine=False, refine_conformers=FINAL_CONFORMERS, smiles=None):
    from chimerax.core.errors import UserError

    if not session.ui.is_gui:
        raise UserError('ligandsearch is a GUI tool; run it in a graphical ChimeraX session.')
    try:
        if sources not in (None, 'ccd', 'ccd_chebi'):
            raise ValueError('sources must be ccd or ccd_chebi')
        validate_conformer_limit(conformers)
        validate_conformer_limit(refine_conformers, MAX_FINAL_CONFORMERS, 'refine_conformers')
        groups = parse_groups(groups)
        group_match = parse_group_match(group_match)
        if smiles is not None:
            smiles = validate_smiles(smiles)
            if codes is not None or groups != ('all',) or sources not in (None, 'ccd'):
                raise ValueError('smiles fits one ligand; omit codes, groups and sources.')
    except ValueError as error:
        raise UserError(str(error)) from error
    if not math.isfinite(resolution) or not .5 <= resolution <= 20:
        raise UserError('resolution must be 0.5-20 Angstroms.')
    TopHits(top)
    if isinstance(workers, bool) or not isinstance(workers, int) or not 0 <= workers <= 256:
        raise UserError('workers must be 0 (automatic) or an integer from 1 to 256.')
    browser = getattr(session, '_trimmings_ligand_browser', None)
    if browser is not None and not browser._closed and browser.job is not None:
        browser.tool_window.shown = True
        raise UserError('Ligand search is already running; use its Stop button first.')
    if browser is not None and not browser._closed and getattr(browser, '_implementation', None) is not IMPLEMENTATION:
        browser.delete()
        browser = None
    if browser is None or browser._closed:
        browser = session._trimmings_ligand_browser = create_browser(session)
    browser.resolution.setValue(resolution)
    browser.top.setValue(top)
    browser.workers.setValue(workers)
    browser.include_observed = observed
    browser.frequency.setChecked(frequency)
    browser.fast.setChecked(fast)
    browser.torsions.setChecked(torsions)
    browser.conformers.setValue(conformers)
    browser.refine.setChecked(refine)
    browser.refine_conformers.setValue(refine_conformers)
    if library is not None:
        browser.library.setText(str(Path(library).expanduser().resolve()))
    browser.group_match.setCurrentIndex(browser.group_match.findData(group_match))
    source = 'smiles' if smiles else sources or DEFAULT_LIBRARY_SOURCE
    source_index = browser.sources.findData(source)
    if browser.sources.currentIndex() == source_index:
        browser.refresh_group_counts()
    else:
        browser.sources.setCurrentIndex(source_index)
    if smiles is not None:
        browser.smiles.setText(smiles)
    browser.set_groups(groups)
    if codes is not None:
        browser.codes.setText(codes)
    if map_model is not None:
        browser.maps.setCurrentIndex(browser.maps.findData(map_model))
    browser.tool_window.shown = True
    if map_model is not None:
        browser.begin_pick()
    return browser


def ligand_stop(session):
    browser = getattr(session, '_trimmings_ligand_browser', None)
    if browser is not None:
        browser.stop()


def cli():
    import argparse
    import signal

    parser = argparse.ArgumentParser(description='Build the cached CCD ligand-search index.')
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--build-library', type=Path)
    action.add_argument('--update-frequency', type=Path)
    action.add_argument('--update-groups', type=Path)
    action.add_argument('--build-chebi', type=Path)
    action.add_argument('--package-library', type=Path, help='Package a completed combined library cache for a release')
    action.add_argument('--install-library', type=Path, help='Install a prebuilt archive into an empty cache folder')
    parser.add_argument('--output', type=Path, help='Output archive for --package-library')
    parser.add_argument('--asset-url', help='Future HTTPS release archive URL for --package-library')
    parser.add_argument('--cache', type=Path, help='Empty target cache for --install-library')
    parser.add_argument('--ontology', type=Path)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--refresh', action='store_true')
    args = parser.parse_args()
    def stopped(*args):
        raise InterruptedError('Stopped by user.')
    signal.signal(signal.SIGTERM, stopped)
    try:
        report = lambda message: print(message, flush=True)
        if args.package_library:
            if args.output is None or args.asset_url is None:
                parser.error('--package-library requires --output and --asset-url')
            print(json.dumps(package_library(args.package_library, args.output, args.asset_url, report=report)))
            return
        if args.install_library:
            if args.cache is None:
                parser.error('--install-library requires --cache pointing to an empty folder')
            install_library_archive(args.cache, args.install_library, report=report)
            return
        if args.build_chebi:
            chebi_api().build(args.build_chebi, report=report, ontology=args.ontology, refresh=args.refresh)
            return
        if args.update_groups:
            update_groups(args.update_groups, report=report, refresh=args.refresh, ontology=args.ontology)
            return
        if args.build_library:
            if args.source is None and library_info(args.build_library) is None:
                from urllib.error import URLError

                try:
                    if install_prebuilt_library(args.build_library, report=report):
                        return
                except (URLError, TimeoutError) as error:
                    raise SystemExit(f'Could not download the combined library ({error}). Retry to resume, or '
                                     'download the archive manually and choose Library > Install library archive.') from error
            print(json.dumps(build_library(args.build_library, args.source, report=report, refresh=args.refresh)))
        if args.update_frequency or args.source is None:
            try:
                update_frequencies(args.update_frequency or args.build_library, report=report)
            except InterruptedError:
                raise
            except Exception as error:
                if args.update_frequency:
                    raise
                report(f'Warning: PDB frequency update failed: {error}. Existing counts retained, or volume-only ordering used.')
        if args.build_library and args.source is None:
            try:
                update_groups(args.build_library, report=report, refresh=args.refresh)
            except InterruptedError:
                raise
            except Exception as error:
                report(f'Warning: Chemical-group update failed: {error}. Existing tags retained; groups all remains available.')
    except InterruptedError:
        print('CCD library update stopped; previous index retained.', flush=True)


if __name__ == '__main__':
    cli()
