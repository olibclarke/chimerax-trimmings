"""Fetch and superimpose mapped domain exemplars without changing the query."""
import json
import math
from pathlib import Path
import re
import tempfile
import time
from types import SimpleNamespace
from urllib.request import Request, urlopen


def source_spec(record, structure):
    """Return only explicit source identities and declared domain boundaries."""
    match = structure.get('domain_match', {})
    if match.get('locations'):
        return {'provider': 'pdb', 'id': structure['accession'].lower(),
                'database': match['database'], 'accession': match['accession'],
                'regions': match.get('regions', []), 'mapped': True}
    reference = record.get('structural_match', {}).get('reference', {})
    same_url = bool(structure.get('url')) and structure['url'] == reference.get('structure_url')
    if not same_url and structure.get('accession') not in (
            reference.get('pdb'), reference.get('ecod_domain'), reference.get('id'), reference.get('uniprot')):
        return None
    def source(provider, identifier, **extra):
        spec = {'provider': provider, 'id': identifier, 'regions': regions, 'mapped': False, **extra}
        uid = reference.get('ecod_uid')
        if reference.get('catalog') == 'ECOD' and isinstance(uid, int) and uid > 0:
            spec['ecod_uid'] = uid
        return spec
    def domain_only():
        uid = reference.get('ecod_uid')
        return source('ecod', str(uid), domain_only=True) if (
            reference.get('catalog') == 'ECOD' and isinstance(uid,int) and uid > 0) else None
    # ECOD's predicted-domain rows omit PDB/chain/UniProt fields. Its AlphaFold
    # identifiers carry the source accession and fragment explicitly.
    predicted = None
    if reference.get('catalog') == 'ECOD' and reference.get('structure_source') == 'predicted_or_other':
        predicted = re.fullmatch(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})_F1_nD[1-9][0-9]*',
                                 reference.get('ecod_domain') or '')
    chain_id = reference.get('chain') or ('A' if predicted else None)
    regions = []
    text = reference.get('ecod_seqid_range', '')
    if text:
        pieces = {}
        for item in text.split(','):
            parsed = re.fullmatch(r'(?:(.+):)?(\d+)-(\d+)', item.strip())
            if not parsed or int(parsed[2]) > int(parsed[3]):
                return domain_only()
            chain = parsed[1] or chain_id
            if not chain:
                return domain_only()
            pieces.setdefault(chain, []).append({'start': int(parsed[2]), 'end': int(parsed[3])})
        regions = [{'chain': chain, 'segments': segments} for chain, segments in pieces.items()]
    elif reference.get('segments') and reference.get('chain') and not reference.get('parent_target'):
        regions = [{'chain': reference['chain'], 'segments': [{'start': a, 'end': b}
                   for a, b in reference['segments']]}]
    if not regions:
        return domain_only()
    if reference.get('pdb'):
        return source('pdb', reference['pdb'].lower())
    if reference.get('structure_source') == 'alphafold' and reference.get('uniprot'):
        return source('alphafold', reference['uniprot'])
    if predicted:
        return source('alphafold', predicted[1], fragment=1)
    return domain_only()


def _coordinate_file_valid(path):
    with Path(path).open('rb') as stream:
        data = stream.read(65536).lstrip()
    if Path(path).suffix == '.pdb':
        return data.startswith(b'ATOM  ') or b'\nATOM  ' in data
    return data.startswith(b'data_')


def _download(url, path, progress):
    request = Request(url, headers={'User-Agent': 'ChimeraX-Trimmings/1.0'})
    temporary = None
    try:
        with urlopen(request, timeout=45) as response, tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            total, size, updated = int(response.headers.get('Content-Length', 0)), 0, 0
            while block := response.read(1024 * 1024):
                size += len(block)
                if size > 1024 * 1024 * 1024:
                    raise ValueError('Structure download exceeds 1 GB.')
                stream.write(block)
                if progress and time.monotonic() - updated > 1:
                    progress(f'Downloading {size / 1048576:.1f}' + (f'/{total / 1048576:.1f}' if total else '') + ' MB')
                    updated = time.monotonic()
        if total and size != total:
            raise ValueError('Incomplete structure download.')
        # Validate according to the destination format before publishing it.
        with temporary.open('rb') as stream:
            data = stream.read(65536).lstrip()
        valid = (data.startswith(b'ATOM  ') or b'\nATOM  ' in data) if path.suffix == '.pdb' else data.startswith(b'data_')
        if not valid:
            raise ValueError('The source did not return a PDB structure.' if path.suffix == '.pdb'
                             else 'The source did not return an mmCIF structure.')
        temporary.replace(path)
    finally:
        if temporary:
            temporary.unlink(missing_ok=True)


def fetch_source(spec, cache, progress=None):
    """Network/cache work only; safe to run outside the GUI thread."""
    spec = dict(spec)
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    provider, identifier = spec['provider'], spec['id']
    if provider == 'ecod':
        uid = spec.get('ecod_uid')
        if not isinstance(uid,int) or uid <= 0:
            raise ValueError('Invalid ECOD domain UID.')
        path = cache / f'ecod-{uid}.pdb'
        if not (path.is_file() and 0 <= time.time()-path.stat().st_mtime < 14*86400 and _coordinate_file_valid(path)):
            _download(f'http://prodata.swmed.edu/ecod/api/domain/{uid}/pdb',path,progress)
        return path,spec
    if provider == 'pdb':
        if not re.fullmatch(r'[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]+', identifier):
            raise ValueError('Invalid PDB accession.')
        url = f'https://files.rcsb.org/download/{identifier.upper()}.cif'
    elif provider == 'alphafold':
        if not re.fullmatch(r'[A-Z0-9]{6,10}', identifier):
            raise ValueError('Invalid AlphaFold accession.')
        with urlopen(Request(f'https://alphafold.ebi.ac.uk/api/prediction/{identifier}',
                             headers={'User-Agent': 'ChimeraX-Trimmings/1.0'}), timeout=30) as response:
            predictions = json.load(response)
        first = min(p['start'] for r in spec['regions'] for p in r['segments'])
        last = max(p['end'] for r in spec['regions'] for p in r['segments'])
        expected = f"AF-{identifier}-F{spec['fragment']}" if spec.get('fragment') else None
        model = next((m for m in predictions if m.get('uniprotAccession') == identifier and m.get('cifUrl')
                      and not m.get('isComplex')
                      and (expected is None or m.get('modelEntityId',m.get('entryId')) == expected)
                      and int(m.get('uniprotStart',m.get('sequenceStart',1))) <= first
                      and int(m.get('uniprotEnd',m.get('sequenceEnd',0))) >= last), None)
        if model is None:
            raise ValueError('No full AlphaFold model is available for this reference.')
        url = model['cifUrl']
        if not url.startswith('https://alphafold.ebi.ac.uk/'):
            raise ValueError('Unexpected AlphaFold coordinate source.')
        spec['sequence_offset'] = int(model.get('uniprotStart', model.get('sequenceStart', 1))) - 1
        # The public prediction service is authoritative for the full-model chain.
        if model.get('chainId'):
            spec['regions'] = [{**region, 'chain': model['chainId']} for region in spec['regions']]
    else:
        raise ValueError('Unsupported structure source.')
    suffix = f"-F{spec['fragment']}" if spec.get('fragment') else ''
    if provider == 'alphafold':
        suffix += '-v' + str(model.get('latestVersion', 'current'))
    path = cache / f'{provider}-{identifier}{suffix}.cif'
    usable = path.is_file() and 0 <= time.time() - path.stat().st_mtime < 14 * 86400
    if usable:
        usable = _coordinate_file_valid(path)
    if not usable:
        _download(url, path, progress)
    if spec.get('ecod_uid'):
        uid = spec['ecod_uid']
        domain_path = cache / f'ecod-{uid}.pdb'
        usable = (domain_path.is_file() and 0 <= time.time()-domain_path.stat().st_mtime < 14*86400
                  and _coordinate_file_valid(domain_path))
        if not usable:
            _download(f'http://prodata.swmed.edu/ecod/api/domain/{uid}/pdb', domain_path, progress)
        spec['domain_path'] = str(domain_path)
    return path, spec


def region_residues(model, region, offset=0):
    from chimerax.atomic import Residues

    candidates = [c for c in model.chains if c.chain_id == region['chain']]
    if not candidates:
        candidates = [c for c in model.chains if any(r.mmcif_chain_id == region['chain'] for r in c.existing_residues)]
    selected = []
    for chain in candidates:
        for segment in region['segments']:
            a, b = segment['start'] - offset, segment['end'] - offset
            if not 1 <= a <= b <= len(chain.residues):
                continue
            members = [r for r in chain.residues[a-1:b] if r is not None]
            # Explicit author endpoints disambiguate structures whose sequence
            # numbering differs from the API mapping. Never guess a chain.
            for index, field in [(a-1, 'author_start'), (b-1, 'author_end')]:
                expected = segment.get(field)
                residue = chain.residues[index]
                if expected is not None and residue is not None and str(residue.number) + residue.insertion_code.strip() != str(expected):
                    raise ValueError('The exemplar residue numbering disagrees with its domain mapping.')
            selected.extend(members)
    return Residues(list(dict.fromkeys(selected)))


def copy_region(session, model, residues, name):
    """Copy only selected atoms/bonds, avoiding a whole-complex temporary copy."""
    from chimerax.atomic import AtomicStructure

    result = AtomicStructure(session, name=name, auto_style=False, log_info=False)
    atom_map = {}
    try:
        for residue in residues:
            new = result.new_residue(residue.name, residue.chain_id, residue.number, residue.insertion_code)
            new.ss_type, new.ss_id = residue.ss_type, residue.ss_id
            new.ribbon_color = residue.ribbon_color
            for atom in residue.atoms:
                copied = result.new_atom(atom.name, atom.element)
                new.add_atom(copied)
                copied.coord, copied.bfactor, copied.occupancy = atom.coord, atom.bfactor, atom.occupancy
                copied.color = atom.color
                atom_map[atom] = copied
        for bond in residues.atoms.intra_bonds:
            a, b = bond.atoms
            result.new_bond(atom_map[a], atom_map[b])
        result.scene_position = model.scene_position
        result.atoms.displays = False
        result.residues.ribbon_displays = True
        return result
    except BaseException:
        result.delete()
        raise


def align_full(session, record, full, spec):
    from chimerax.match_maker import match, defaults, CP_SPECIFIC_SPECIFIC

    baseline, best = full.scene_position, None
    failures = []
    for region in spec['regions']:
        full.scene_position = baseline
        try:
            selected = region_residues(full, region, spec.get('sequence_offset', 0))
            if len(selected) < 3:
                continue
            chains = selected.unique_chains
            if len(chains) != 1:
                raise ValueError('A domain occurrence must map to one polymer chain.')
            fits = match(session, CP_SPECIFIC_SPECIFIC, [(record['chain'], chains[0])],
                         defaults['matrix'], defaults['alignment_algorithm'], defaults['gap_open'], defaults['gap_extend'],
                         domain_residues=(record['residues'], selected), cutoff_distance=defaults['iter_cutoff'],
                         always_raise_errors=True, log_parameters=False, verbose=None, compute_ss=False)
            if not fits:
                continue
            fit = fits[0]
            count, rmsd = len(fit['final match atoms']), fit['final RMSD']
            if count < 3 or not math.isfinite(rmsd):
                continue
            rank = (-count, rmsd)
            if best is None or rank < best[0]:
                best = rank, full.scene_position, selected, count, rmsd
        except Exception as error:
            failures.append(str(error))
    if best is None:
        full.scene_position = baseline
        raise ValueError('No mapped occurrence could be aligned: ' + (failures[0] if failures else 'fewer than three matching residues.'))
    _, full.scene_position, selected, count, rmsd = best
    return selected, count, rmsd


def open_aligned(session, record, spec, path):
    from chimerax.core.models import Model
    from chimerax.mmcif import open_mmcif

    group = Model(f"{spec['id'].upper()} aligned to {record['chain_spec']} {record['name']}", session)
    models = []
    try:
        full = None
        if not spec.get('domain_only'):
            models, _ = open_mmcif(session, str(path), auto_style=False, log_info=False, max_models=1)
            if not models or models[0].num_atoms == 0:
                raise ValueError('The coordinate file contains no atoms.')
            full = models[0]
            full.name = f"{spec['id'].upper()} full context"
            full.display = False
        session.models.add([group])
        if models:
            session.models.add(models, parent=group)
        if full is not None and not full.ss_assigned:
            from chimerax.dssp import compute_ss
            compute_ss(full)
        if spec.get('domain_path') or spec.get('domain_only'):
            from chimerax.pdb import open_pdb
            domains, _ = open_pdb(session, str(path) if spec.get('domain_only') else spec['domain_path'],
                                  auto_style=False, log_info=False, max_models=1)
            if not domains or not domains[0].num_atoms:
                raise ValueError('ECOD returned an empty domain structure.')
            domain = domains[0]
            # Parent immediately so a failed fit closes every new model.
            domain.display = False
            session.models.add(domains, parent=group)
            if not domain.ss_assigned:
                from chimerax.dssp import compute_ss
                compute_ss(domain)
            domain_spec = {'regions': [{'chain': c.chain_id, 'segments': [{'start':1,'end':len(c.residues)}]}
                                       for c in domain.chains if len(c.residues)]}
            members, count, rmsd = align_full(session, record, domain, domain_spec)
            # Fit context to the exact domain. This also handles source files
            # whose coordinate frame differs, without moving the query.
            context_rmsd = 0
            if full is not None:
                selected, _, context_rmsd = align_full(session, {'chain':members.unique_chains[0], 'residues':members}, full, spec)
                style_reference(full, selected)
            style_reference(domain, members)
            crop = copy_region(session, domain, members, f"ECOD {spec['ecod_uid']} aligned domain")
            session.models.add([crop], parent=group)
            crop.scene_position = domain.scene_position
            session.models.close(domains)
            relation = full.scene_position.inverse() * crop.scene_position if full is not None else None
        else:
            selected, count, rmsd = align_full(session, record, full, spec)
            style_reference(full, selected)
            crop = copy_region(session, full, selected, f"{spec['id'].upper()} aligned domain")
            session.models.add([crop], parent=group)
            crop.scene_position = full.scene_position
            relation, context_rmsd = None, 0
        return SimpleNamespace(group=group, domain=crop, full=full, count=count, rmsd=rmsd,
                               source=spec, context_relation=relation, context_rmsd=context_rmsd)
    except BaseException:
        if group.id is not None:
            session.models.close([group])
        else:
            for model in models:
                model.delete()
            group.delete()
        raise


def style_reference(full, selected):
    import numpy as np
    full.atoms.displays = False
    full.residues.ribbon_displays = True
    full.atoms.colors = np.tile(np.array([180,180,180,255], dtype=np.uint8), (len(full.atoms),1))
    full.residues.ribbon_colors = np.tile(np.array([180,180,180,255], dtype=np.uint8), (len(full.residues),1))
    selected.atoms.colors = np.tile(np.array([255,190,70,255], dtype=np.uint8), (len(selected.atoms),1))
    selected.ribbon_colors = np.tile(np.array([255,190,70,255], dtype=np.uint8), (len(selected),1))


def realign_pair(session, record, pair):
    relation = getattr(pair, 'context_relation', None)
    if relation is not None or pair.full is None:
        original = pair.domain.scene_position
        try:
            regions = [{'chain': c.chain_id, 'segments':[{'start':1,'end':len(c.residues)}]}
                       for c in pair.domain.chains if len(c.residues)]
            _, count, rmsd = align_full(session, record, pair.domain, {'regions':regions})
            if pair.full is not None:
                pair.full.scene_position = pair.domain.scene_position * relation.inverse()
            pair.count, pair.rmsd = count, rmsd
        except BaseException:
            pair.domain.scene_position = original
            raise
        return
    showing, original = pair.full.display, pair.full.scene_position
    crop = None
    pair.full.display = False
    try:
        selected, count, rmsd = align_full(session, record, pair.full, pair.source)
        style_reference(pair.full, selected)
        crop = copy_region(session, pair.full, selected, pair.domain.name)
        session.models.add([crop], parent=pair.group)
        crop.scene_position = pair.full.scene_position
        session.models.close([pair.domain])
        pair.domain, pair.count, pair.rmsd = crop, count, rmsd
        pair.full.display, pair.domain.display = showing, not showing
    except BaseException:
        if crop is not None and crop is not pair.domain and not crop.deleted:
            if crop.id is not None:
                session.models.close([crop])
            else:
                crop.delete()
        pair.full.scene_position, pair.full.display = original, showing
        raise


def show_context(pair, show):
    if pair.full is None:
        raise ValueError('ECOD provides this domain, but its full source structure is not identified.')
    relation = getattr(pair, 'context_relation', None)
    if show:
        pair.full.scene_position = (pair.domain.scene_position * relation.inverse()
                                    if relation is not None else pair.domain.scene_position)
    else:
        pair.domain.scene_position = (pair.full.scene_position * relation
                                      if relation is not None else pair.full.scene_position)
    pair.full.display, pair.domain.display = bool(show), not show


def focus_pair(session, pair):
    """Frame the visible exemplar, updating clipping after a context switch."""
    model = pair.full if pair.full is not None and pair.full.display else pair.domain
    bounds = model.atoms.scene_bounds
    if bounds is None:
        return
    view, center = session.main_view, bounds.center()
    view.view_all(bounds, pad=0.15)
    normal, depth = view.camera.view_direction(), max(bounds.radius(), 1.0)
    view.clip_planes.set_clip_position('near', center - depth * normal, view)
    view.clip_planes.set_clip_position('far', center + depth * normal, view)
    view.center_of_rotation = center
