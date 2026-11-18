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
SCHEMA = '1'
VOLUME_LEVELS = (.15, .30, .50)
VOLUME_RESOLUTION = 2.6
STANDARD_EXTRA_STARTS, STANDARD_SAMPLES = 24, 3000
FAST_EXTRA_STARTS, FAST_SAMPLES = 8, 1500
DEFAULT_CONFORMERS, MAX_CONFORMERS = 8, 32
CONFORMER_RMSD = .25  # Ignore torsion variants indistinguishable after rigid alignment.
INTERNAL_OVERLAP = 1.6  # Reject newly introduced gross heavy-atom overlaps, in Angstroms.
SHAPE = None
IMPLEMENTATION = object()


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
    # A small, reusable descriptor permits immediate queueing of a full library.
    # It is a broad scheduling prior, not a molecular identity or fit score.
    xyz = np.array(conformers[0]['xyz'], dtype=np.float32)
    xyz -= xyz.mean(axis=0)
    try:
        grid = shape_api()._shape_grid(xyz, np.array(elements, np.float32), VOLUME_RESOLUTION)
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


def _public_json(url, cancel, payload=None):
    """Retry transient public-data errors, without adding a third-party client."""
    import urllib.error
    import urllib.request

    data = None if payload is None else json.dumps(payload).encode('utf-8')
    headers = {'User-Agent': 'ChimeraX-ligandsearch/1', 'Content-Type': 'application/json'}
    for attempt in range(3):
        _cancelled(cancel)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=headers), timeout=30) as response:
                return json.load(response)
        except (OSError, ValueError) as error:
            if isinstance(error, urllib.error.HTTPError) and error.code not in (429, 500, 502, 503, 504):
                raise
            if attempt == 2:
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


def prepare_blob(data, center, resolution, level, radius=14.):
    """Pick a connected component near CoFR without changing the Volume/surface."""
    import numpy as np
    from scipy.ndimage import label, distance_transform_edt, generate_binary_structure

    if not all(math.isfinite(v) and v > 0 for v in (resolution, level, radius)):
        raise ValueError('Resolution, contour level and region radius must be positive finite numbers.')
    if resolution < 2 * max(np.linalg.norm(data.ijk_to_xyz_transform.matrix[:, :3], axis=0)):
        raise ValueError('Resolution is finer than map sampling; supply the actual map resolution in Angstroms.')
    lower, upper = shape_api()._shape_cofr_crop_bounds(data, center, radius)
    if np.any(upper <= lower):
        raise ValueError('CoFR is outside this map; center on the unknown blob and try again.')
    if np.prod(upper-lower) > 4_000_000:
        raise ValueError('This region is too large; reduce Region radius, or use a coarser-sampled input map.')
    matrix = data.matrix(tuple(lower), tuple(upper-lower)).astype(np.float32, copy=True)
    if not np.isfinite(matrix).all():
        raise ValueError('The selected map region contains nonfinite density values.')
    affine = data.ijk_to_xyz_transform.matrix[:, :3]
    metric = affine.T @ affine
    if not np.allclose(metric, np.diag(np.diag(metric)), atol=np.max(np.diag(metric)) * 1e-7):
        raise ValueError('This version needs an orthogonal map grid; resample the map onto an orthogonal grid first.')
    labels, _ = label(matrix >= level, generate_binary_structure(3, 1))
    k, j, i = np.indices(matrix.shape)
    points = data.ijk_to_xyz_transform.transform_points(
        np.column_stack((i.ravel(), j.ravel(), k.ravel())).astype(float) + lower)
    nearby = (np.linalg.norm(points-center, axis=1) <= 3) & (labels.ravel() > 0)
    if not np.any(nearby):
        raise ValueError('No blob within 3 A of CoFR at this contour. Recenter on the blob or lower Blob contour.')
    ids, counts = np.unique(labels.ravel()[nearby], return_counts=True)
    blob = labels == ids[np.argmax(counts)]
    if any(np.any(np.take(blob, edge, axis=axis)) for edge in (0, -1) for axis in range(3)):
        raise ValueError('The blob reaches the region boundary. Raise Blob contour to separate it from '
                         'nearby density, or increase Region radius if the isolated blob is larger.')
    steps = np.linalg.norm(affine, axis=0)
    outside = distance_transform_edt(~blob, sampling=steps[::-1])
    soft = np.exp(-.5 * (outside / (resolution / 3)) ** 2)
    soft[outside > resolution] = 0
    return dict(matrix=(matrix * soft).astype(np.float32), raw=matrix, mask=outside <= resolution,
                blob=blob, origin=np.array(data.ijk_to_xyz_transform * lower), step=np.array(data.step),
                rotation=np.array(data.rotation), cell_angles=np.array(data.cell_angles), resolution=resolution,
                level=level, volume=float(blob.sum() * abs(np.linalg.det(affine))),
                fraction=float(level / matrix[blob].max()))


class TopHits:
    """A bounded, per-CCD leaderboard; completion order cannot break tied ranks."""
    def __init__(self, limit=20):
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
            raise ValueError('Top-hit buffer must contain 1 to 100 entries.')
        self.limit, self.rows = limit, []

    def add(self, row):
        if not math.isfinite(row['cc']) or not math.isfinite(row.get('score', row['cc'])):
            return
        previous = next((old for old in self.rows if old['id'] == row['id']), None)
        if previous is not None and previous.get('score', previous['cc']) >= row.get('score', row['cc']):
            return
        self.rows = [old for old in self.rows if old['id'] != row['id']]
        self.rows.append(row)
        self.rows.sort(key=lambda item: (-item.get('score', item['cc']), item['id']))
        del self.rows[self.limit:]


class _Quiet:
    def status(self, *args, **kwargs):
        pass


def selected_conformers(record, observed=False):
    conformers = record['conformers'] if observed else [c for c in record['conformers'] if c['source'] == 'ideal']
    return conformers or record['conformers'][:1]


def validate_conformer_limit(limit):
    if isinstance(limit, bool) or not isinstance(limit, int) or not 2 <= limit <= MAX_CONFORMERS:
        raise ValueError(f'conformers must be an integer from 2 to {MAX_CONFORMERS} (total per ligand).')


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


def sparse_conformers(record, observed=False, limit=DEFAULT_CONFORMERS):
    """Bounded single/multiple torsion perturbations; no Cartesian product search."""
    import numpy as np
    from scipy.spatial import cKDTree
    from chimerax.geometry import align_points

    validate_conformer_limit(limit)
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

    def fit(self, record, observed=False, fast=False, torsions=False, conformers=DEFAULT_CONFORMERS):
        import numpy as np
        from chimerax.geometry import Place
        from chimerax.core.errors import UserError

        options = dict(resolution=self.config['resolution'], radius=None, objective='density',
                       max_points=FAST_SAMPLES if fast else STANDARD_SAMPLES,
                       extra_starts=FAST_EXTRA_STARTS if fast else STANDARD_EXTRA_STARTS,
                       fallback_correlation=0., bandpass=False)
        best = None
        tested, skipped = 0, 0
        variants = sparse_conformers(record, observed, conformers) if torsions else selected_conformers(record, observed)
        for conformer in variants:
            source = self.source(conformer['xyz'], record['elements'])
            try:
                _, stats = self.api._fit_shape_original_map_search(source, self.target, options, None, collect=True)
            except (UserError, self.api._ShapeEmptyMapSearch):
                if not conformer.get('torsion_changes'):
                    raise
                # An expanded trial can exceed a grid limit or have no fit.
                # Keep the successful original conformer rather than losing its CCD.
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
            raise ValueError('No usable coordinate conformer.')
        best.update(conformers_tested=tested, conformers_skipped=skipped)
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
                result = dict(result=fitter.fit(job['component'], job['observed'], fast=job.get('fast', False),
                                               torsions=job.get('torsions', False),
                                               conformers=job.get('conformers', DEFAULT_CONFORMERS)))
            except Exception as error:
                result = dict(error=f'{type(error).__name__}: {error}')
            result['id'] = job['component']['id']
            print('LIGAND_RESULT:' + json.dumps(result), flush=True)
    finally:
        fitter.close()


class SearchJob:
    """Bounded worker dispatch. Stop terminates outstanding work, not saved hits."""
    def __init__(self, cache, config, workers=0, observed=False, ids=None, frequency=True, fast=False,
                 torsions=False, conformers=DEFAULT_CONFORMERS):
        validate_conformer_limit(conformers)
        self.cache, self.config, self.requested = Path(cache).expanduser().resolve(), config, workers
        self.observed, self.ids = observed, None if ids is None else set(ids)
        self.frequency = frequency
        self.fast = fast
        self.torsions, self.conformers = torsions, conformers
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
        import numpy as np

        processes, threads = [], []
        stopped = threading.Event()
        waiting, completed = queue.Queue(), queue.Queue()
        done, failures = 0, 0
        try:
            info = library_info(self.cache)
            if info is None:
                raise ValueError('CCD library missing. Click Expand/Update CCD library first.')
            with closing(sqlite3.connect((self.cache / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)) as db, \
                    tempfile.TemporaryDirectory(prefix='ligand_search_') as temporary:
                occurrence = frequency_info(self.cache) if self.frequency else None
                counts = None if occurrence is None else occurrence['counts']
                if self.frequency and occurrence is None:
                    self.emit(dict(kind='warning', text='PDB occurrence cache unavailable; using volume-only ordering. '
                                   'Expand/Update CCD library to fetch occurrence counts.'))
                ordered = library_order(db, self.config['volume'], self.config['fraction'], counts)
                if self.ids is not None:
                    ordered = [row for row in ordered if row[1] in self.ids]
                total = len(ordered)
                if not total:
                    raise ValueError('No requested components are present in the CCD library.')
                api = shape_api()
                count = min(total, max(1, api._shape_available_cpus() - 1) if self.requested == 0 else self.requested,
                            api._shape_available_cpus())
                mode = 'fast' if self.fast else 'standard'
                detail = f'; sparse torsions, up to {self.conformers} conformers/ligand' if self.torsions else ''
                self.emit(dict(kind='status', text=f'Searching {total:,} CCD entries with {count} workers ({mode} search{detail})',
                               total=total, done=0))
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
                                    if 'result' in result:
                                        result['result'].update(estimated_volume=job['estimated_volume'],
                                                                search_priority=job['search_priority'],
                                                                pdb_entries=job['pdb_entries'])
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
                        record = json.loads(db.execute('SELECT payload FROM components WHERE id=?', (ident,)).fetchone()[0])
                        waiting.put(dict(component=record, observed=self.observed, fast=self.fast,
                                         torsions=self.torsions, conformers=self.conformers,
                                         estimated_volume=estimated, search_priority=prior,
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
                        done += 1
                        failures += 'error' in message
                        self.emit(dict(kind='result', **message, done=done, total=total, failures=failures))
                        dispatch()
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
            self.emit(dict(kind='finished', done=done, total=total, failures=failures))
        except InterruptedError:
            pass
        except Exception as error:
            self.emit(dict(kind='error', text=str(error)))


class LibraryJob:
    """Build in a separate native process so parsing/downloading cannot stall Qt."""
    def __init__(self, cache, source=None):
        self.cache, self.source = Path(cache).expanduser(), source
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
            command = [executable, '-u', str(Path(__file__).resolve()), '--build-library', str(self.cache), '--refresh']
            if self.source is not None:
                command.extend(('--source', str(self.source)))
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
                                    'Counting PDB occurrences:', 'PDB frequencies ready:')):
                    self.events.put(dict(kind='status', text=line))
                elif line.startswith('Warning:'):
                    self.events.put(dict(kind='warning', text=line[len('Warning:'):].strip()))
            if process.wait() != 0:
                raise RuntimeError('CCD library update failed: ' + ' '.join(tail))
            self.events.put(dict(kind='library', info=library_info(self.cache)))
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
        model = AtomicStructure(session, name=f'CCD {row["id"]}', auto_style=False, log_info=False)
    else:
        old_residues = model.residues
    residue = None
    try:
        # Deleting a structure's last residue can delete the native structure.
        # Build the replacement first, then discard only the old residues.
        residue = model.new_residue(row['id'], 'A', 1 if fresh else max(old_residues.numbers) + 1)
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
        model.name = f'CCD {row["id"]}: {row["name"]}'
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
    from Qt.QtCore import Qt, QTimer
    from Qt.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                              QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox,
                              QTableWidget, QTableWidgetItem, QVBoxLayout, QAbstractItemView)
    from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
    from chimerax.core.tools import ToolInstance
    from chimerax.ui import MainToolWindow
    from chimerax.map import Volume

    class LigandBrowser(ToolInstance):
        SESSION_SAVE = False
        SESSION_ENDURING = False

        def __init__(self):
            super().__init__(session, 'Ligand Search')
            self._closed, self.job, self.preview = False, None, None
            self._implementation = IMPLEMENTATION
            self.current, self.volume, self.watched_data = None, None, None
            self.blob_volume = 0.
            self.buffer, self.failures = TopHits(20), []
            self.seed, self._pick_mode, self._previous_mode = None, None, None
            self.stale, self.pending_start, self.started, self.library_completed = False, False, 0., False
            self.tool_window = MainToolWindow(self)
            root = QVBoxLayout(self.tool_window.ui_area)
            self.summary = QLabel('Identify candidate CCD components for an isolated density blob.')
            self.summary.setWordWrap(True)
            root.addWidget(self.summary)
            form = QFormLayout()
            self.maps = QComboBox()
            self.maps.currentIndexChanged.connect(self.map_changed)
            form.addRow('Density map', self.maps)
            self.resolution = QDoubleSpinBox()
            self.resolution.setRange(.5, 20.)
            self.resolution.setDecimals(2)
            self.resolution.setValue(2.6)
            self.resolution.setSuffix(' A')
            form.addRow('Local resolution', self.resolution)
            self.level = QDoubleSpinBox()
            self.level.setRange(.00000001, 1e6)
            self.level.setDecimals(8)
            self.level.setValue(1.)
            form.addRow('Blob contour', self.level)
            self.radius = QDoubleSpinBox()
            self.radius.setRange(3., 100.)
            self.radius.setValue(14.)
            self.radius.setSuffix(' A')
            form.addRow('Region radius', self.radius)
            self.top = QSpinBox()
            self.top.setRange(1, 100)
            self.top.setValue(20)
            self.top.setToolTip('Maximum retained hits. One best pose per CCD code.')
            form.addRow('Top-hit buffer', self.top)
            self.workers = QSpinBox()
            self.workers.setRange(0, 256)
            self.workers.setSpecialValueText('Automatic')
            form.addRow('CPU workers', self.workers)
            self.codes = QLineEdit()
            self.codes.setPlaceholderText('All indexed CCD entries, or comma-separated codes')
            form.addRow('CCD codes (optional)', self.codes)
            self.library = QLineEdit(str(default_cache()))
            form.addRow('Library cache folder', self.library)
            self.fast = QCheckBox(f'Fast search ({FAST_EXTRA_STARTS} extra starts, {FAST_SAMPLES:,} samples)')
            self.fast.setToolTip('Tests the same CCD entries with fewer orientation starts and refinement samples. '
                                 'Final density scoring is unchanged, but poses and rankings can differ. '
                                 'Unchecked uses the standard search settings.')
            form.addRow(self.fast)
            self.torsions = QCheckBox('Sparse torsion sampling')
            self.torsions.setToolTip('Generate a small deterministic set of additional conformers around flexible '
                                     'single bonds. Slower; not exhaustive or energy-minimized. No receptor clash tests.')
            form.addRow(self.torsions)
            self.conformers = QSpinBox()
            self.conformers.setRange(2, MAX_CONFORMERS)
            self.conformers.setValue(DEFAULT_CONFORMERS)
            self.conformers.setToolTip('Maximum total per CCD entry, including selected ideal/observed conformers. '
                                       'Rigid ligands or rejected/duplicate rotations may yield fewer.')
            self.conformers_label = QLabel('Max conformers per ligand')
            form.addRow(self.conformers_label, self.conformers)
            self.observed = QCheckBox('Also fit the CCD observed-coordinate conformer')
            self.observed.setToolTip('Adds the CCD example structure, not your bound ligand pose. Slower.')
            form.addRow(self.observed)
            self.torsions.toggled.connect(self.update_conformer_controls)
            self.frequency = QCheckBox('Prioritize volume match and PDB occurrence frequency')
            self.frequency.setChecked(True)
            self.frequency.setToolTip('Common components get a capped priority boost. All components remain eligible; '
                                      'frequency does not affect density scores or the final ranking.')
            form.addRow(self.frequency)
            root.addLayout(form)
            location = QHBoxLayout()
            self.pick = QPushButton('Pick blob')
            self.pick.clicked.connect(self.begin_pick)
            location.addWidget(self.pick)
            cofr = QPushButton('Use CoFR')
            cofr.clicked.connect(self.use_cofr)
            location.addWidget(cofr)
            self.location = QLabel('Seed: current center of rotation')
            location.addWidget(self.location)
            root.addLayout(location)
            actions = QHBoxLayout()
            self.start_button, self.stop_button, self.update_button = (
                QPushButton('Start search'), QPushButton('Stop'), QPushButton('Expand/Update CCD library'))
            self.start_button.clicked.connect(self.start_search)
            self.stop_button.clicked.connect(self.stop)
            self.update_button.clicked.connect(lambda: self.start_library())
            for button in (self.start_button, self.stop_button, self.update_button):
                actions.addWidget(button)
            root.addLayout(actions)
            self.table = QTableWidget(0, 8)
            self.table.setHorizontalHeaderLabels(('CCD', 'Score', 'Density CC', 'Blob covered', 'Atom support', 'Volume ratio', 'PDB entries', 'Name'))
            self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
            self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
            self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            self.table.itemSelectionChanged.connect(self.activate_row)
            self.table.horizontalHeader().setStretchLastSection(True)
            root.addWidget(self.table, 1)
            self.details = QLabel('Density ranking is a hypothesis, not chemical identification. No clash filtering.')
            self.details.setWordWrap(True)
            self.details.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            root.addWidget(self.details)
            navigation = QHBoxLayout()
            self.follow_top = QCheckBox('Follow top hit')
            self.follow_top.setChecked(True)
            self.follow_top.setToolTip('Update the density preview when the best hit changes. '
                                       'Choosing another hit pauses following; check this to resume.')
            self.follow_top.toggled.connect(lambda checked: self.refresh_rows())
            navigation.addWidget(self.follow_top)
            for text, function in (('Previous', lambda: self.step(-1)), ('Next', lambda: self.step(1)),
                                   ('Focus', self.focus), ('Keep hit', self.keep), ('Export CSV', self.export)):
                button = QPushButton(text)
                button.clicked.connect(function)
                navigation.addWidget(button)
            root.addLayout(navigation)
            self.status = QLabel('Ready')
            self.status.setWordWrap(True)
            root.addWidget(self.status)
            self.timer = QTimer(self.tool_window.ui_area)
            self.timer.setInterval(200)
            self.timer.timeout.connect(self.poll)
            self.handlers = [session.triggers.add_handler(ADD_MODELS, self.models_changed),
                             session.triggers.add_handler(REMOVE_MODELS, self.models_changed)]
            self.refresh_maps()
            self.set_busy(False)
            self.tool_window.manage('side')

        def set_busy(self, busy):
            for widget in (self.start_button, self.update_button, self.resolution, self.radius,
                           self.level, self.top, self.workers, self.fast, self.torsions, self.frequency, self.codes, self.library, self.maps, self.pick):
                widget.setEnabled(not busy)
            self.update_conformer_controls()
            self.stop_button.setEnabled(busy)

        def update_conformer_controls(self, *args):
            enabled = self.torsions.isEnabled() and self.torsions.isChecked()
            for widget in (self.conformers_label, self.conformers, self.observed):
                widget.setEnabled(enabled)

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
            if selected is not self.maps.currentData():
                self.map_changed()

        def map_changed(self, *args):
            self.seed = None
            model = self.maps.currentData()
            if model is not None and not model.deleted:
                levels = [surface.level for surface in model.surfaces if surface.level > 0]
                if levels:
                    self.level.setValue(min(levels))
            self.location.setText('Seed: current center of rotation')

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
            self.status.setText(text)

        def data_changed(self, reason):
            if reason in ('values changed', 'coordinates changed'):
                self.invalidate('Map data changed; results are stale. Restart the search.')

        def use_cofr(self):
            self.restore_pick()
            self.seed = None
            self.location.setText('Seed: current center of rotation')

        def begin_pick(self):
            from chimerax.mouse_modes import MouseMode
            from chimerax.pick_blobs.pickblobs import _picked_triangle
            import numpy as np

            browser = self
            class Pick(MouseMode):
                name = 'ligand search blob'
                def mouse_down(mode, event):
                    try:
                        surface, triangle = _picked_triangle(session.main_view.picked_object(*event.position()))
                        owner = next((m for m in surface.drawing_lineage if isinstance(m, Volume)), None) if surface else None
                        if owner is None:
                            browser.status.setText('Click a density-map surface, not an atomic model.')
                            return
                        browser.maps.setCurrentIndex(browser.maps.findData(owner))
                        point = np.mean(surface.vertices[surface.triangles[triangle]], axis=0)
                        scene = surface.scene_position * point
                        browser.seed = owner.scene_position.inverse() * scene
                        browser.location.setText('Seed: picked density surface point')
                        browser.status.setText('Blob selected. Click Start search.')
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
            session.logger.warning('Ligand search: ' + str(error))

        def start_search(self):
            if self.job is not None:
                return
            try:
                model = self.maps.currentData()
                if model is None or model.deleted:
                    raise ValueError('Open a density map, then choose it in Density map.')
                info = library_info(Path(self.library.text()).expanduser().resolve())
                if info is None:
                    self.pending_start = True
                    self.start_library()
                    return
                self.restore_pick()
                center = model.scene_position.inverse() * session.main_view.center_of_rotation if self.seed is None else self.seed
                config = prepare_blob(model.data, center, self.resolution.value(), self.level.value(), self.radius.value())
                self.volume = model
                self.blob_volume = config['volume']
                self.unwatch()
                self.watched_data = model.data
                self.watched_data.add_change_callback(self.data_changed)
                self.remove_preview()
                self.buffer, self.current, self.failures = TopHits(self.top.value()), None, []
                self.stale = False
                self.follow_top.setChecked(True)
                self.details.setText('Searching; density scores are provisional. No clash filtering.')
                self.refresh_rows()
                mode = 'Fast search' if self.fast.isChecked() else 'Standard search'
                if self.torsions.isChecked():
                    mode += f' with sparse torsions (up to {self.conformers.value()} conformers/ligand)'
                self.summary.setText(f'{mode}; blob volume {config["volume"]:.3g} A^3; keeping top {self.top.value()} CCD hits. '
                                     'Volume and optional PDB frequency prioritize the queue; ranking uses density fit. No clash tests.')
                codes = [v.strip().upper() for v in self.codes.text().replace(' ', ',').split(',') if v.strip()] or None
                self.job = SearchJob(Path(self.library.text()).expanduser().resolve(), config, self.workers.value(),
                                     self.observed.isChecked(), codes, self.frequency.isChecked(), fast=self.fast.isChecked(),
                                     torsions=self.torsions.isChecked(), conformers=self.conformers.value())
                self.begin_job()
            except Exception as error:
                self.report_error(error)

        def start_library(self, source=None):
            if self.job is not None:
                return
            self.job = LibraryJob(Path(self.library.text()).expanduser().resolve(), source)
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
                    if 'result' in event:
                        self.buffer.add(event['result'])
                        changed = True
                    elif len(self.failures) < 10:
                        self.failures.append(f'{event["id"]}: {event["error"]}')
                    last_status = (f'{event["done"]:,}/{event["total"]:,} checked; '
                                   f'{event["failures"]:,} failed; {time.monotonic()-self.started:.0f}s elapsed')
                elif kind == 'status':
                    last_status = event['text']
                elif kind == 'warning':
                    session.logger.warning('Ligand search: ' + event['text'])
                    self.summary.setText(self.summary.text() + ' ' + event['text'])
                elif kind == 'library':
                    library_ready = event['info']
                    self.library_completed = True
                    last_status = f'CCD library ready: {library_ready["count"]} usable components'
                elif kind == 'error':
                    finish_error = event['text']
                    self.pending_start = False
                elif kind == 'finished':
                    last_status = (f'Finished: {event["done"]:,}/{event["total"]:,} checked; '
                                   f'{event["failures"]:,} failed; {time.monotonic()-self.started:.0f}s')
            if changed:
                self.refresh_rows()
            if last_status and not self.stale:
                self.status.setText(last_status)
                session.logger.status('Ligand search: ' + last_status, blank_after=0)
            if finish_error:
                self.report_error(finish_error)
            if not job.thread.is_alive() and job.events.empty():
                self.job = None
                self.timer.stop()
                self.set_busy(False)
                if job.cancel.is_set() and not self.stale:
                    self.status.setText(f'Stopped; {len(self.buffer.rows)} top hits retained.')
                if self.failures:
                    session.logger.warning('Ligand search skipped some candidates. Examples:\n' + '\n'.join(self.failures))
                if self.library_completed and self.pending_start:
                    self.pending_start = False
                    self.start_search()

        def refresh_rows(self):
            selected = self.current['id'] if self.current else None
            if (self.follow_top.isChecked() or self.current is None) and self.buffer.rows:
                selected = self.buffer.rows[0]['id']
            self.table.blockSignals(True)
            self.table.setRowCount(len(self.buffer.rows))
            index = -1
            for row_index, hit in enumerate(self.buffer.rows):
                values = (hit['id'], f'{hit.get("score", hit["cc"]):.3f}', f'{hit["cc"]:.3f}', f'{100*hit["coverage"]:.0f}%',
                          f'{100*hit["supported"]:.0f}%',
                          f'{hit["estimated_volume"] / self.blob_volume:.2f}'
                          if self.blob_volume > 0 and 'estimated_volume' in hit else '',
                          str(hit['pdb_entries']) if hit.get('pdb_entries') is not None else '', hit['name'])
                for column, value in enumerate(values):
                    self.table.setItem(row_index, column, QTableWidgetItem(value))
                if hit['id'] == selected:
                    index = row_index
            self.table.clearSelection()
            if index >= 0:
                self.table.selectRow(index)
            self.table.blockSignals(False)
            # Do not rebuild the preview for worse incoming hits. Object identity
            # also catches an improved pose for the same CCD component.
            if index >= 0 and (self.current is not self.buffer.rows[index] or
                               self.preview is None or self.preview.deleted):
                self.activate_row(manual=False)
            elif self.current is not None and index < 0:
                self.show_details(pinned=True)

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
            self.current = self.buffer.rows[rows[0].row()]
            try:
                if self.volume is None or self.volume.deleted:
                    raise ValueError('The fitting map is no longer available.')
                self.preview = make_model(session, self.current, self.volume, self.preview)
                self.preview.SESSION_SAVE = False
                self.show_details()
            except Exception as error:
                self.report_error(error)

        def show_details(self, pinned=False):
            if self.current is not None:
                hit = self.current
                trials = f'Conformers tested: {hit.get("conformers_tested", 1)}. '
                if hit.get('conformers_skipped'):
                    trials += f'Trials skipped: {hit["conformers_skipped"]}. '
                self.details.setText(f'{hit["id"]}: {hit["name"]}\n{hit["formula"]}; {hit["source"]} coordinates. '
                                     f'{trials}'
                                     f'CC {hit["cc"]:.3f}; {100*hit["outside"]:.0f}% predicted contour outside scoring region. '
                                     'No clash filtering. Rankings are provisional while the search runs.' +
                                     (' Preview pinned outside the current top-hit buffer.' if pinned else ''))

        def step(self, direction):
            if self.buffer.rows and not self.stale:
                index = next((i for i, row in enumerate(self.buffer.rows)
                              if self.current is not None and row['id'] == self.current['id']), -1)
                self.follow_top.blockSignals(True)
                self.follow_top.setChecked(False)
                self.follow_top.blockSignals(False)
                self.table.selectRow((index + direction) % len(self.buffer.rows))

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
                session.logger.info(f'Kept CCD {self.current["id"]} as #{kept.id_string}; target models unchanged.')
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
                                     'torsion_sampling', 'conformers_tested', 'conformers_skipped', 'torsion_changes', 'stale'))
                    for row in self.buffer.rows:
                        writer.writerow((row['id'], row['name'], row['formula'], f'{row.get("score", row["cc"]):.3g}', f'{row["cc"]:.3g}',
                                         f'{row["coverage"]:.3g}', f'{row["supported"]:.3g}',
                                         f'{row["outside"]:.3g}', row['source'], row.get('pdb_entries'),
                                         row.get('fast_search', False), row.get('torsion_sampling', False),
                                         row.get('conformers_tested', 1), row.get('conformers_skipped', 0),
                                         '; '.join(f'{change["bond"]}:{change["degrees"]:+g}'
                                                   for change in row.get('torsion_changes', [])), self.stale))
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
            self.timer.stop()
            self.unwatch()
            for handler in self.handlers:
                handler.remove()
            self.remove_preview()
            if getattr(session, '_trimmings_ligand_browser', None) is self:
                session._trimmings_ligand_browser = None
            super().delete()

    return LigandBrowser()


def ligand_search(session, map_model=None, resolution=2.6, level=None, radius=14., top=20,
                  library=None, workers=0, observed=False, codes=None, frequency=True, fast=False,
                  torsions=False, conformers=DEFAULT_CONFORMERS):
    from chimerax.core.errors import UserError

    if not session.ui.is_gui:
        raise UserError('ligandsearch is a GUI tool; run it in a graphical ChimeraX session.')
    try:
        validate_conformer_limit(conformers)
    except ValueError as error:
        raise UserError(str(error)) from error
    if not math.isfinite(resolution) or not .5 <= resolution <= 20 or not math.isfinite(radius) or not 3 <= radius <= 100:
        raise UserError('resolution must be 0.5-20 Angstroms; radius must be 3-100 Angstroms.')
    if level is not None and (not math.isfinite(level) or not 1e-8 <= level <= 1e6):
        raise UserError('level must be a positive finite map contour value from 1e-8 to 1e6.')
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
    browser.radius.setValue(radius)
    browser.top.setValue(top)
    browser.workers.setValue(workers)
    browser.observed.setChecked(observed)
    browser.frequency.setChecked(frequency)
    browser.fast.setChecked(fast)
    browser.torsions.setChecked(torsions)
    browser.conformers.setValue(conformers)
    if library is not None:
        browser.library.setText(str(Path(library).expanduser().resolve()))
    if codes is not None:
        browser.codes.setText(codes)
    if map_model is not None:
        browser.maps.setCurrentIndex(browser.maps.findData(map_model))
    if level is not None:
        browser.level.setValue(level)
    browser.tool_window.shown = True
    if map_model is not None:
        browser.start_search()
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
    parser.add_argument('--source', type=Path)
    parser.add_argument('--refresh', action='store_true')
    args = parser.parse_args()
    def stopped(*args):
        raise InterruptedError('Stopped by user.')
    signal.signal(signal.SIGTERM, stopped)
    try:
        report = lambda message: print(message, flush=True)
        if args.build_library:
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
    except InterruptedError:
        print('CCD library update stopped; previous index retained.', flush=True)


if __name__ == '__main__':
    cli()
