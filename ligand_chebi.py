"""Optional ChEBI metabolite library; only validated 3D coordinates are fitted.

This runs in the library subprocess, never on the ChimeraX GUI thread. Exact
InChIKeys preserve stereochemistry and protonation when reusing coordinates.
"""

from contextlib import closing
from pathlib import Path
import gzip
import json
import os
import re
import sqlite3
import threading
import time

# Injected by ligand_search.chebi_api so runscript installs need no sys.path edits.
ls = None

PUBCHEM = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
SCHEMA = '1'


def info(cache):
    path = Path(cache) / 'chebi.sqlite'
    if not path.is_file():
        return None
    with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as db:
        metadata = dict(db.execute('SELECT key,value FROM metadata'))
        if metadata.get('schema') != SCHEMA:
            return None
        return dict(metadata, count=str(db.execute('SELECT count(*) FROM components').fetchone()[0]))


def group_counts(cache):
    """Union memberships so a reused CCD identity is never counted twice."""
    cache = Path(cache).resolve()
    with closing(sqlite3.connect((cache / 'groups.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        db.execute('ATTACH DATABASE ? AS expansion', ((cache / 'chebi.sqlite').as_uri() + '?mode=ro',))
        combined = ('WITH combined AS (SELECT comp_id,group_key FROM tags UNION '
                    "SELECT c.id,j.value FROM expansion.components c,json_each(c.payload,'$.chemical_groups') j), "
                    'ids AS (SELECT id FROM components UNION SELECT id FROM expansion.components) ')
        counts = dict(db.execute(combined + 'SELECT group_key,count(*) FROM combined GROUP BY group_key'))
        counts['all'] = db.execute(combined + 'SELECT count(*) FROM ids').fetchone()[0]
        counts['unclassified'] = db.execute(combined + 'SELECT count(*) FROM ids WHERE NOT EXISTS '
                                           '(SELECT 1 FROM combined WHERE comp_id=ids.id)').fetchone()[0]
        return counts


def candidates(ontology, cancel):
    molecules = {}
    with gzip.open(ontology, 'rt', encoding='utf-8') as stream:
        masks, _, version = ls.ontology_group_masks(stream, cancel, molecules=molecules)
    bit = 1 << next(i for i, row in enumerate(ls.CHEMICAL_GROUPS) if row[0] == 'metabolites')
    by_identity = {}
    for ident, molecule in molecules.items():
        if not masks[ident] & bit:
            continue
        key = molecule['inchikey']
        if re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]', key) is None:
            continue
        entry = by_identity.setdefault(key, dict(molecule, chebi_ids=[], chemical_groups=set(), group_evidence=[]))
        entry['chebi_ids'].append(ident)
        groups = [row[0] for i, row in enumerate(ls.CHEMICAL_GROUPS) if masks[ident] & (1 << i)]
        entry['chemical_groups'].update(groups)
        entry['group_evidence'].extend(dict(group=group, source='ChEBI ontology', id=ident) for group in groups)
    for entry in by_identity.values():
        entry['chebi_ids'].sort(key=lambda ident: int(ident.split(':')[1]))
        entry['chemical_groups'] = sorted(entry['chemical_groups'])
    return sorted(by_identity.values(), key=lambda row: int(row['chebi_ids'][0].split(':')[1])), version


def ccd_identities(cache, cancel, report):
    """Cache exact identities in batches; no connectivity-only coordinate reuse."""
    path = Path(cache) / 'identities.sqlite'
    with closing(sqlite3.connect(path)) as db:
        db.execute('CREATE TABLE IF NOT EXISTS identities(id TEXT PRIMARY KEY,inchikey TEXT)')
        with closing(sqlite3.connect((Path(cache) / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)) as source:
            ids = [row[0] for row in source.execute('SELECT id FROM components ORDER BY id')]
        known = dict(db.execute('SELECT id,inchikey FROM identities'))
        missing = [ident for ident in ids if ident not in known or not known[ident]]
        for offset in range(0, len(missing), 500):
            ls._cancelled(cancel)
            batch = missing[offset:offset+500]
            query = '{chem_comps(comp_ids:' + json.dumps(batch) + '){rcsb_id rcsb_chem_comp_descriptor{InChIKey}}}'
            result = ls._public_json(ls.PDB_DATA_URL, cancel, {'query': query})
            rows = (result.get('data') or {}).get('chem_comps')
            if result.get('errors') or not isinstance(rows, list):
                raise ValueError('CCD identity lookup failed; retry the ChEBI library update.')
            for row in rows:
                if row and row['rcsb_id'] in batch:
                    key = (row.get('rcsb_chem_comp_descriptor') or {}).get('InChIKey')
                    if key:
                        known[row['rcsb_id']] = key
                        db.execute('INSERT OR REPLACE INTO identities VALUES(?,?)', (row['rcsb_id'], key))
            db.commit()
            report(f'ChEBI library: CCD identities {min(offset+500, len(missing)):,}/{len(missing):,}')
        reverse = {}
        for ident in ids:
            if known.get(ident):
                reverse.setdefault(known[ident], ident)
        return reverse


def pubchem_record(entry, properties, compound):
    """Reject 2D records, changed identities, disconnected salts and bad geometry."""
    import numpy as np

    if properties['InChIKey'] != entry['inchikey']:
        raise ValueError('PubChem returned a different chemical identity')
    if compound['id']['id']['cid'] != properties['CID']:
        raise ValueError('PubChem returned coordinates for a different CID')
    if entry.get('charge') is not None and int(entry['charge']) != properties['Charge']:
        raise ValueError('PubChem returned a different charge state')
    coordinates = next((c for c in compound['coords'] if 2 in c['type']), None)
    if coordinates is None:
        raise ValueError('No 3D conformer')
    atoms = compound['atoms']
    elements = dict(zip(atoms['aid'], atoms['element'], strict=True))
    heavy = [aid for aid in atoms['aid'] if elements[aid] > 1]
    if len(heavy) < 4 or any(not 1 <= element <= 118 for element in elements.values()):
        raise ValueError('Too few heavy atoms or unsupported elements')
    index = {aid: i for i, aid in enumerate(heavy)}
    conf = coordinates['conformers'][0]
    points = dict(zip(coordinates['aid'], zip(conf['x'], conf['y'], conf['z'], strict=True), strict=True))
    xyz = np.array([points[aid] for aid in heavy], float)
    if not np.isfinite(xyz).all() or np.linalg.matrix_rank(xyz-xyz.mean(axis=0)) < 2:
        raise ValueError('Invalid conformer geometry')
    bonds, adjacency = [], [set() for _ in heavy]
    source = compound.get('bonds', {})
    for a, b, order in zip(source.get('aid1', []), source.get('aid2', []), source.get('order', []), strict=True):
        if a in index and b in index:
            if order not in (1, 2, 3, 4):
                raise ValueError('Unsupported bond order')
            i, j = index[a], index[b]
            bonds.append((i, j, {1: 'SING', 2: 'DOUB', 3: 'TRIP', 4: 'AROM'}[order]))
            adjacency[i].add(j)
            adjacency[j].add(i)
    visited, pending = set(), [0]
    while pending:
        i = pending.pop()
        if i not in visited:
            visited.add(i)
            pending.extend(adjacency[i]-visited)
    if len(visited) != len(heavy):
        raise ValueError('Disconnected chemical components')
    charge = sum(atom['value'] for atom in atoms.get('charge', []))
    if charge != properties['Charge']:
        raise ValueError('3D conformer charge does not match the identity record')
    return dict(id=entry['chebi_ids'][0], name=entry['name'], formula=properties['MolecularFormula'],
                type='NON-POLYMER', names=[f'A{i+1}' for i in range(len(heavy))],
                elements=[elements[aid] for aid in heavy], bonds=bonds,
                conformers=[dict(source='PubChem3D', xyz=xyz.tolist())],
                library_source='ChEBI', pubchem_cid=properties['CID'], inchikey=entry['inchikey'],
                formal_charge=charge)


def fetch_record(entry, cache, cancel):
    """Resume cached successes/definite absences; retry transient failures later."""
    import urllib.error

    folder = Path(cache) / 'pubchem3d'
    folder.mkdir(exist_ok=True)
    path = folder / (entry['inchikey'] + '.json')
    if path.is_file():
        try:
            saved = json.loads(path.read_text())
            record = saved.get('record')
            if record and record['inchikey'] == entry['inchikey'] and record.get('conformers'):
                return record, ''
            # Retry definite absences monthly, not on every resume.
            if record is None and saved.get('reason') and time.time()-path.stat().st_mtime < 30*86400:
                return None, saved['reason']
        except (OSError, ValueError, KeyError, AttributeError, TypeError):
            pass  # A damaged cache entry is fetched again, not deferred forever.
    try:
        properties = ls._public_json(PUBCHEM + '/compound/inchikey/' + entry['inchikey'] +
                                     '/property/InChIKey,MolecularFormula,Charge/JSON', cancel)
    except urllib.error.HTTPError as error:
        if error.code != 404:
            raise
        properties = {}
    records = properties.get('PropertyTable', {}).get('Properties', [])
    records = sorted((row for row in records if row.get('InChIKey') == entry['inchikey']), key=lambda row: row['CID'])
    reason, record = 'No exact PubChem identity with usable 3D coordinates', None
    for props in records:
        ls._cancelled(cancel)
        try:
            result = ls._public_json(PUBCHEM + f'/compound/cid/{props["CID"]}/JSON?record_type=3d', cancel)
            record = pubchem_record(entry, props, result['PC_Compounds'][0])
            break
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
        except (KeyError, ValueError, IndexError) as error:
            reason = str(error)
        cancel.wait(.3)
    # Coordinate payloads retain provenance. Scheduling volumes are computed by
    # the native library process, not by a web service.
    temporary = path.with_suffix('.part')
    temporary.write_text(json.dumps(dict(record=record, reason='' if record else reason)))
    os.replace(temporary, path)
    return record, reason


def build(cache, cancel=None, report=lambda message: None, ontology=None, refresh=False):
    cache = Path(cache).expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)
    lock = cache / 'chebi_build.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        try:
            os.kill(int(lock.read_text()), 0)
        except (ValueError, ProcessLookupError):
            lock.unlink(missing_ok=True)
            return build(cache, cancel, report, ontology, refresh)
        raise RuntimeError('Another ChEBI library update is running; wait for it or stop that update.') from None
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(str(os.getpid()))
        return _build(cache, cancel, report, ontology, refresh)
    finally:
        lock.unlink(missing_ok=True)


def _build(cache, cancel=None, report=lambda message: None, ontology=None, refresh=False):
    """Incremental, resumable expansion. Completed entries survive Stop/errors."""
    cache = Path(cache).expanduser().resolve()
    cancel = cancel or threading.Event()
    if ls.library_info(cache) is None:
        raise ValueError('Build the CCD library first.')
    ontology = ontology or ls._group_ontology(cache, cancel, report, refresh)
    entries, version = candidates(ontology, cancel)
    if not entries:
        raise ValueError('No defined ChEBI metabolites found; previous expansion retained.')
    report(f'ChEBI library: {len(entries):,} distinct metabolite identities; preparing CCD matches')
    identities = ccd_identities(cache, cancel, report)
    with closing(sqlite3.connect(cache / 'chebi.sqlite')) as db, \
            closing(sqlite3.connect((cache / 'ccd.sqlite').as_uri() + '?mode=ro', uri=True)) as ccd:
        db.executescript('CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT);'
                        'CREATE TABLE IF NOT EXISTS components(id TEXT PRIMARY KEY,name TEXT,atoms INTEGER,'
                        'v0 REAL,v1 REAL,v2 REAL,payload TEXT);')
        db.executemany('INSERT OR REPLACE INTO metadata VALUES(?,?)',
                       dict(schema=SCHEMA, state='building', ontology=version, total=str(len(entries))).items())
        db.commit()
        added = reused = unavailable = errors = 0
        retained_ids = set()
        for index, entry in enumerate(entries):
            ls._cancelled(cancel)
            try:
                ident = identities.get(entry['inchikey'])
                if ident:
                    record = json.loads(ccd.execute('SELECT payload FROM components WHERE id=?', (ident,)).fetchone()[0])
                    record.update(library_source='CCD', inchikey=entry['inchikey'])
                    reused += 1
                else:
                    # One request stream, deliberately below PubChem's public
                    # rate limit; it is independent of the fitting worker count.
                    record, reason = fetch_record(entry, cache, cancel)
                    cancel.wait(.4)
                    if record is None:
                        unavailable += 1
                        db.execute('DELETE FROM components WHERE id=?', (entry['chebi_ids'][0],))
                        db.commit()
                        continue
                    record = ls.describe_component(record)
                    added += 1
                groups = set(entry['chemical_groups'])
                evidence = list(entry.get('group_evidence') or [dict(group=key, source='ChEBI ontology', id=entry['chebi_ids'][0])
                                                              for key in groups])
                if ls.has_carbon_nitrogen_bond(record):
                    groups.add('organic_nitrogen')
                    evidence.append(dict(group='organic_nitrogen', source='Bond connectivity', id='C-N'))
                record.update(chebi_ids=entry['chebi_ids'], chemical_groups=sorted(groups),
                              group_evidence=evidence,
                              classification_status='tagged')
                retained_ids.add(record['id'])
                db.execute('DELETE FROM components WHERE id=?', (entry['chebi_ids'][0],))
                db.execute('INSERT OR REPLACE INTO components VALUES(?,?,?,?,?,?,?)',
                           (record['id'], record['name'], len(record['names']), *record['volumes'], json.dumps(record)))
                db.commit()
            except InterruptedError:
                raise
            except Exception as error:
                errors += 1
                if errors <= 5:
                    report(f'Warning: {entry["chebi_ids"][0]} deferred: {error}')
                if errors >= 10 and added+reused == 0:
                    raise RuntimeError('ChEBI/PubChem services unavailable; retry the expansion later. Cached entries retained.') from error
            finally:
                report(f'ChEBI library: {index+1:,}/{len(entries):,}; {added:,} new 3D, {reused:,} CCD matches, '
                       f'{unavailable:,} unavailable, {errors:,} deferred')
        if errors == 0:
            # A successful ontology refresh must not retain withdrawn entries or
            # metabolites that no longer have a supported coordinate identity.
            obsolete = [row[0] for row in db.execute('SELECT id FROM components') if row[0] not in retained_ids]
            db.executemany('DELETE FROM components WHERE id=?', ((ident,) for ident in obsolete))
        db.executemany('INSERT OR REPLACE INTO metadata VALUES(?,?)',
                       dict(state='complete' if errors == 0 else 'partial', added=str(added), reused=str(reused),
                            unavailable=str(unavailable), errors=str(errors),
                            updated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())).items())
        db.commit()
    return info(cache)
