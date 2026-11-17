#!/usr/bin/env python3
"""Complete a folded-domain Pfam checkpoint using InterPro-linked structures.

Inventory and candidate attempts persist separately from verified templates.
Visits are breadth-first: one candidate per unresolved family in each round.
Downloads are shared through the builder cache; no incomplete release is published.
"""

from collections import Counter
from contextlib import closing
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import random
import re
import sqlite3
from threading import Event
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse

import build_pfam_library as builder

VERSION = 1
FOLDED_TYPES = builder.FOLDED_PROFILE_TYPES


class CandidateRejected(ValueError):
    def __init__(self, reason, message):
        self.reason = reason
        super().__init__(message)


def eligible_family(family):
    # A short profile may have longer natural insertions, so its length alone
    # does not establish that no usable >=20-residue template can exist.
    return family['type'] in FOLDED_TYPES


def inventory_database(path):
    state = sqlite3.connect(path)
    state.executescript('''
        CREATE TABLE IF NOT EXISTS candidates(family TEXT,key TEXT,data TEXT,
            PRIMARY KEY(family,key));
        CREATE TABLE IF NOT EXISTS pages(family TEXT,source TEXT,url TEXT,next TEXT,total INTEGER,
            checked TEXT,PRIMARY KEY(family,source,url));
        CREATE TABLE IF NOT EXISTS attempts(family TEXT,key TEXT,version INTEGER,status TEXT,
            detail TEXT,checked TEXT,PRIMARY KEY(family,key,version));
    ''')
    return state


def candidate_key(candidate):
    return hashlib.sha256(json.dumps(candidate, sort_keys=True).encode()).hexdigest()


def source_url(accession, source, page_size):
    if source == 'alphafold':
        return f'{builder.fs.INTERPRO}/protein/uniprot/entry/pfam/{accession}/?has_model=true&page_size={page_size}&extra_fields=sequence'
    return f'{builder.fs.INTERPRO}/structure/pdb/entry/pfam/{accession}/?page_size={page_size}'


def safe_page(url, accession, source):
    parsed = urlparse(url)
    expected = '/interpro/api/' + ('protein/uniprot' if source == 'alphafold' else 'structure/pdb') + '/entry/pfam/' + accession + '/'
    if parsed.scheme != 'https' or parsed.hostname != 'www.ebi.ac.uk' or parsed.path.lower() != expected.lower():
        raise ValueError('Unsafe InterPro pagination URL.')
    return url


def parse_candidates(payload, accession, source):
    if payload is None:
        return []
    if not isinstance(payload, dict) or not isinstance(payload.get('results'), list):
        raise ValueError('Unrecognized InterPro candidate response.')
    candidates = []
    for item in payload['results']:
        meta = item['metadata']
        identifier = meta['accession']
        if source == 'alphafold':
            if not builder.UNIPROT.fullmatch(identifier) or meta.get('in_alphafold') is False:
                continue
        elif not re.fullmatch(r'[A-Za-z0-9]{4}', identifier):
            raise ValueError('Invalid PDB accession in InterPro response.')
        for entry in item.get('entries', []):
            if entry.get('accession', '').upper() != accession:
                continue
            sequence = (item.get('extra_fields', {}).get('sequence') if source == 'alphafold' else entry.get('sequence'))
            if isinstance(sequence, dict):
                sequence = sequence.get('sequence')
            if not isinstance(sequence, str) or not re.fullmatch(r'[A-Za-z]+', sequence):
                continue
            sequence = sequence.upper()
            locations = entry.get('entry_protein_locations', []) if source == 'alphafold' else entry.get('entry_structure_locations', [])
            for location in locations:
                pieces = sorted((fragment['start'], fragment['end']) for fragment in location.get('fragments', []))
                if not pieces or any(type(a) is not int or type(b) is not int or not 1 <= a <= b <= len(sequence)
                                     for a, b in pieces):
                    raise ValueError('Invalid InterPro domain boundaries.')
                if any(b >= c for (_, b), (c, _) in zip(pieces, pieces[1:])):
                    raise ValueError('Overlapping InterPro domain fragments.')
                score = location.get('score')
                candidate = {'source': source, 'id': identifier, 'sequence': sequence, 'segments': pieces,
                             'chain': entry.get('chain', '') if source == 'pdb' else '',
                             'representative': bool(location.get('representative')), 'annotation_score': score,
                             'title': meta.get('name', identifier), 'resolution': meta.get('resolution')}
                if source == 'pdb' and not candidate['chain']:
                    continue
                candidates.append(candidate)
    return candidates


def rank(candidate, family):
    length = sum(b - a + 1 for a, b in candidate['segments'])
    unknown = sum(candidate['sequence'][a - 1:b].count('X') for a, b in candidate['segments']) / max(1, length)
    return (candidate['source'] != 'pdb', unknown, not candidate.get('representative', False),
            abs(length - family['length']) / max(1, family['length']),
            float(candidate.get('resolution') or 99), candidate['id'])


def retryable(error):
    if isinstance(error, HTTPError):
        return error.code in (429, 500, 502, 503, 504)
    return isinstance(error, (URLError, TimeoutError, ConnectionError))


def load_candidates(family, database, cache, page_size, advance=False):
    accession = family['accession']
    with closing(inventory_database(database)) as state:
        for source in ['pdb', 'alphafold']:
            pages = state.execute('SELECT url,next FROM pages WHERE family=? AND source=? ORDER BY rowid',
                                  (accession, source)).fetchall()
            url = source_url(accession, source, page_size) if not pages else pages[-1][1] if advance else None
            if url is None:
                continue
            safe_page(url, accession, source)
            if any(url == previous for previous, _ in pages):
                raise ValueError('Cyclic InterPro candidate pagination.')
            payload = cache.json(url)
            if payload is None:
                candidates, following, total = [], None, 0
            else:
                candidates = parse_candidates(payload, accession, source)
                following, total = payload.get('next'), payload.get('count')
                if following is not None:
                    safe_page(following, accession, source)
                if type(total) is not int or total < 0:
                    raise ValueError('Invalid InterPro candidate count.')
            for candidate in candidates:
                candidate['annotation_url'] = url
                key = candidate_key(candidate)
                state.execute('INSERT OR IGNORE INTO candidates VALUES(?,?,?)', (accession, key, json.dumps(candidate)))
            state.execute('INSERT INTO pages VALUES(?,?,?,?,?,?)',
                          (accession, source, url, following, total, datetime.now(timezone.utc).isoformat()))
            state.commit()
        rows = state.execute('SELECT key,data FROM candidates WHERE family=? AND key NOT IN '
                             '(SELECT key FROM attempts WHERE family=? AND version=? AND status!="retryable")',
                             (accession, accession, VERSION)).fetchall()
    return sorted(((key, json.loads(data)) for key, data in rows), key=lambda row: rank(row[1], family))


def occurrence_hits(candidate, family, profile, folder, args, cancel):
    """Validate the annotated natural sequence before any coordinate download."""
    sequence, pieces = candidate['sequence'], candidate['segments']
    a, b = max(1, pieces[0][0] - 12), min(len(sequence), pieces[-1][1] + 12)
    accepted = []
    for _, _, first, last, score, evalue, coverage in builder.profile_hits(
            sequence[a - 1:b], family, profile, folder, args, cancel):
        first, last = first + a - 1, last + a - 1
        shared = max(0, min(last, pieces[-1][1]) - max(first, pieces[0][0]) + 1)
        if shared >= .5 * min(last - first + 1, pieces[-1][1] - pieces[0][0] + 1):
            accepted.append((first, last, score, evalue, coverage))
    if not accepted:
        raise CandidateRejected('profile_or_coverage', 'No current gathering-threshold match with sufficient profile coverage at this occurrence.')
    return accepted


@lru_cache(maxsize=4)
def coordinate_chains(data, chain):
    """Reuse parsing for proteins supplying several families; keep memory bounded."""
    return tuple(builder.polymer_chains(data, chain))


def validate_coordinates(data, candidate, family, profile, folder, work, args, cancel, matches=None):
    expected = candidate['sequence']
    pieces = candidate['segments']
    saw_sequence, saw_profile, reasons = False, False, []
    for sequence, residues, chain in coordinate_chains(data, candidate['chain']):
        if sequence != expected:
            continue
        saw_sequence = True
        # Bound sequence revalidation to this annotated occurrence. An unrelated
        # occurrence of the same family elsewhere must not validate this crop.
        hits = matches if matches is not None else occurrence_hits(candidate, family, profile, folder, args, cancel)
        for first, last, score, evalue, coverage in hits:
            saw_profile = True
            segments = [(first, last)] if len(pieces) == 1 else [
                (max(first, x), min(last, y)) for x, y in pieces if max(first, x) <= min(last, y)]
            try:
                # Check modeled coverage of each natural piece; the HMM saw the
                # natural intervening sequence, never an artificial concatenation.
                selected = []
                for x, y in segments:
                    crop, _ = builder.crop_positions(sequence, residues, x, y, args, plddt=False)
                    selected.extend(crop)
                if not args.min_residues <= len(selected) <= 9999:
                    raise CandidateRejected('insufficient_size', f'{len(selected)} modeled residues; need {args.min_residues}–9999.')
                mean = sum(r.find_atom('CA', '*').b_iso for r in selected) / len(selected)
                if candidate['source'] == 'alphafold' and mean < args.min_plddt:
                    raise CandidateRejected('low_confidence', f'Mean domain pLDDT {mean:.1f}; need {args.min_plddt:g}.')
                provenance = {'catalog': 'Pfam', 'pdb': candidate['id'] if candidate['source'] == 'pdb' else None,
                              'chain': chain, 'segments': segments, 'structure_title': candidate['title'],
                              'interpro_annotation_url': candidate['annotation_url'], 'interpro_segments': pieces,
                              'annotation_sequence_sha256': hashlib.sha256(expected.encode()).hexdigest(),
                              'profile_match': {'bitscore': score, 'evalue': evalue, 'hmm_coverage': coverage},
                              'evidence': 'InterPro Pfam occurrence, sequence-exact coordinates, current-profile gathering thresholds',
                              'source_coordinate_sha256': hashlib.sha256(data).hexdigest()}
                if candidate['source'] == 'alphafold':
                    provenance.update(uniprot=candidate['id'], mean_plddt=mean,
                                      structure_url='https://alphafold.ebi.ac.uk/entry/' + candidate['id'])
                else:
                    provenance['structure_url'] = 'https://www.rcsb.org/structure/' + candidate['id']
                return builder.save_reference(family, selected, work, candidate['source'], args, **provenance)
            except ValueError as error:
                reasons.append(error)
    if not saw_sequence:
        raise CandidateRejected('sequence_mismatch', 'Model/polymer sequence differs from the InterPro-annotated sequence.')
    if not saw_profile:
        raise CandidateRejected('profile_or_coverage', 'No current gathering-threshold match with sufficient profile coverage at this occurrence.')
    if reasons:
        error = reasons[-1]
        if isinstance(error, CandidateRejected):
            raise error
        raise CandidateRejected('modeled_coverage', str(error))
    raise CandidateRejected('no_model', 'No usable coordinate segment.')


def validate_candidate(candidate, family, profiles, cache, work, args, cancel):
    key = candidate_key(candidate)
    folder = work / 'interpro-families' / family['accession'] / key[:16]
    folder.mkdir(parents=True, exist_ok=True)
    profile = builder.family_profile(profiles, family, folder, args, cancel)
    matches = occurrence_hits(candidate, family, profile, folder, args, cancel)
    if candidate['source'] == 'pdb':
        url = 'https://files.rcsb.org/download/' + candidate['id'] + '.cif.gz'
        data = cache.get(url)
        if not data:
            raise CandidateRejected('no_model', 'Linked PDB coordinates unavailable.')
        item = validate_coordinates(data, candidate, family, profile, folder, work, args, cancel, matches=matches)
        item[0]['source_coordinate_url'] = url
        return item
    predictions = cache.json('https://alphafold.ebi.ac.uk/api/prediction/' + candidate['id'])
    if predictions is None or predictions == []:
        raise CandidateRejected('no_model', 'Linked AlphaFold prediction unavailable.')
    if not isinstance(predictions, list) or any(not isinstance(p, dict) for p in predictions):
        raise ValueError('Invalid AlphaFold prediction metadata.')
    errors = []
    for prediction in predictions:
        url = prediction.get('cifUrl', '')
        if prediction.get('uniprotAccession') != candidate['id']:
            continue
        parsed = urlparse(url)
        if parsed.scheme != 'https' or parsed.hostname != 'alphafold.ebi.ac.uk':
            raise ValueError('Unsafe AlphaFold coordinate URL.')
        # Reject obsolete sequences before fetching coordinates, when supplied.
        sequence = prediction.get('sequence')
        if sequence is not None and sequence != candidate['sequence']:
            errors.append(CandidateRejected('sequence_mismatch', 'AlphaFold metadata sequence differs from InterPro.'))
            continue
        data = cache.get(url)
        if not data:
            continue
        try:
            item = validate_coordinates(data, candidate, family, profile, folder, work, args, cancel, matches=matches)
            item[0]['source_coordinate_url'] = url
            return item
        except CandidateRejected as error:
            errors.append(error)
    if errors:
        raise errors[-1]
    raise CandidateRejected('no_model', 'No compatible AlphaFold coordinate model.')


def folded_coverage(output, state, families):
    covered = set(builder.represented(state))
    required = {a for a, f in families.items() if eligible_family(f)}
    excluded = {a for a in families if a not in required}
    report = {'scope': 'Pfam Domain, Family and Repeat entries; experimentally or prediction-validated folded templates',
              'release': '39.0', 'required_families': len(required), 'covered_families': len(required & covered),
              'missing_families': len(required - covered), 'folded_scope_complete': required <= covered,
              'all_pfam_families': len(families), 'all_pfam_covered': len(covered & families.keys()),
              'excluded_by_type': dict(Counter(families[a]['type'] for a in excluded)),
              'short_profiles_still_in_scope': sum(families[a]['length'] < 20 for a in required),
              'notes': ['Family/Repeat types can include non-globular proteins; a failed candidate does not establish disorder.',
                        'Short profiles remain eligible if a natural validated folded segment meets the minimum residue count.',
                        'Unresolved entries remain gaps until evidence supports an explicit exclusion; thresholds are unchanged.']}
    builder.atomic_json(output / 'folded-coverage.json', report)
    with (output / 'folded-families.tsv').open('w') as stream:
        stream.write('accession\tname\ttype\tprofile_length\tstatus\n')
        for a, f in sorted(families.items()):
            status = 'out of folded scope' if a in excluded else 'verified template' if a in covered else 'unresolved'
            stream.write(f'{a}\t{f["name"]}\t{f["type"]}\t{f["length"]}\t{status}\n')
    return report


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--work', required=True, type=Path)
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--profiles', required=True, type=Path)
    parser.add_argument('--models', required=True, type=int)
    parser.add_argument('--sample', type=int, default=40, help='0 visits all missing families in folded scope')
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--page-size', type=int, default=20)
    parser.add_argument('--seed', type=int, default=43)
    parser.add_argument('--report', required=True, type=Path)
    parser.add_argument('--publish-if-complete', action='store_true', help='Publish only when every in-scope family has a verified template')
    parser.add_argument('--integrate-into', type=Path, help='Merge into this library only after complete publication')
    parser.add_argument('--make-default', action='store_true')
    parser.add_argument('--archive', type=Path)
    cfg = parser.parse_args(argv)
    if cfg.sample < 0 or min(cfg.rounds, cfg.workers, cfg.page_size) < 1 or cfg.page_size > 200:
        parser.error('Invalid sample/round/worker/page counts.')
    if cfg.make_default and cfg.integrate_into is None:
        parser.error('--make-default requires --integrate-into.')
    if (cfg.integrate_into or cfg.archive) and not cfg.publish_if_complete:
        parser.error('Integration/packaging requires --publish-if-complete.')
    cfg.out, cfg.work, cfg.profiles, cfg.report = (p.resolve() for p in (cfg.out, cfg.work, cfg.profiles, cfg.report))
    state_path = cfg.work / 'state.sqlite'
    if not state_path.is_file():
        raise ValueError('A prepared Pfam checkpoint state.sqlite is required.')
    families = builder.catalog(cfg.metadata, cfg.models)
    cancel = Event()
    args = builder.parser().parse_args(['--out', str(cfg.out), '--release', '39.0', '--scope', 'folded'])
    cache = builder.PublicCache(cfg.out / 'downloads/structures', cancel)
    started = time.monotonic()
    cfg.report.parent.mkdir(parents=True, exist_ok=True)
    inventory = cfg.out / 'interpro-candidates.sqlite'
    report = {'started': datetime.now(timezone.utc).isoformat(), 'strategy': 'ranked InterPro candidates, breadth-first rounds',
              'completed_visits': 0, 'new_templates': 0, 'outcomes': {}, 'rounds_requested': cfg.rounds,
              'criteria': {'minimum_profile_and_modeled_coverage': .8, 'minimum_mean_plddt': 70, 'minimum_residues': 20}}
    with builder.build_lock(cfg.out), closing(builder.state_database(state_path)) as state:
        with closing(inventory_database(inventory)):
            pass
        scope = folded_coverage(cfg.out, state, families)
        covered = set(builder.represented(state))
        pending = sorted(a for a, f in families.items() if eligible_family(f) and a not in covered)
        random.Random(cfg.seed).shuffle(pending)
        if cfg.sample:
            pending = pending[:cfg.sample]
        report.update(initial_coverage=scope, selected_families=len(pending), stage='running')
        outcomes = Counter()
        last_snapshot = 0
        def snapshot(force=False):
            nonlocal last_snapshot
            report.update(elapsed_seconds=time.monotonic() - started, outcomes=dict(outcomes),
                          request_starts=cache.request_starts)
            builder.atomic_json(cfg.report, report)
            builder.atomic_json(cfg.out / 'progress.json', {'phase': 'InterPro folded-domain fill',
                'round': report.get('round', 0), 'completed_visits': report['completed_visits'],
                'selected_families': len(pending), 'new_templates': report['new_templates'],
                'covered_folded_families': scope['covered_families'] + report['new_templates'],
                'required_folded_families': scope['required_families'], 'outcomes': dict(outcomes),
                'elapsed_seconds': report['elapsed_seconds']})
            if force or time.monotonic() - last_snapshot >= 60:
                folded_coverage(cfg.out, state, families)
                builder.coverage_report(cfg.out, state, families, set(families), args)
                last_snapshot = time.monotonic()
        snapshot(True)
        def visit(accession):
            key = None
            try:
                candidates = load_candidates(families[accession], inventory, cache, cfg.page_size)
                if not candidates:
                    candidates = load_candidates(families[accession], inventory, cache, cfg.page_size, advance=True)
                if not candidates:
                    return accession, None, 'no_unused_candidates', None, None
                # Do not let several incomplete PDB chains consume every round
                # before a readily available full predicted model is tested.
                if report['round'] % 2 == 0:
                    candidates = sorted(candidates, key=lambda row: (row[1]['source'] != 'alphafold', rank(row[1], families[accession])))
                key, candidate = candidates[0]
                record = validate_candidate(candidate, families[accession], cfg.profiles, cache, cfg.work, args, cancel)
                return accession, key, 'verified', record, {'source': candidate['source'], 'id': candidate['id']}
            except InterruptedError:
                raise
            except Exception as error:
                status = 'retryable' if retryable(error) else error.reason if isinstance(error, CandidateRejected) else 'source_or_validation_error'
                return accession, key, status, None, {'error': str(error)}
        try:
            for round_number in range(1, cfg.rounds + 1):
                report['round'] = round_number
                covered = set(builder.represented(state))
                active = [a for a in pending if a not in covered]
                if not active:
                    break
                builder.say(f'InterPro round {round_number}/{cfg.rounds}: {len(active)} folded families')
                for _, (accession, key, status, record, detail) in builder.parallel_jobs(active, visit, cfg.workers, cancel):
                    if record is not None and builder.add_reference(state, record):
                        report['new_templates'] += 1
                    now = datetime.now(timezone.utc).isoformat()
                    if key is not None:
                        with closing(inventory_database(inventory)) as ledger:
                            ledger.execute('INSERT OR REPLACE INTO attempts VALUES(?,?,?,?,?,?)',
                                           (accession, key, VERSION, status, json.dumps(detail), now))
                            ledger.commit()
                    report['completed_visits'] += 1
                    outcomes[status] += 1
                    with (cfg.report.parent / 'attempts.jsonl').open('a') as stream:
                        stream.write(json.dumps({'family': accession, 'candidate_key': key, 'status': status,
                                                 'detail': detail, 'checked': now}) + '\n')
                    snapshot()
                    builder.say(f'InterPro {report["completed_visits"]}: {accession} {status}; '
                                f'{report["new_templates"]} new; {report["elapsed_seconds"]:.0f}s')
        except BaseException:
            cancel.set()
            report['stage'] = 'interrupted'
            raise
        else:
            report.update(stage='pass_complete', finished=datetime.now(timezone.utc).isoformat())
        finally:
            snapshot(True)
        final_coverage = builder.coverage_report(cfg.out, state, families, set(families), args)
        report['final_coverage'] = {k:v for k,v in final_coverage.items() if k not in ('excluded_family_accessions', 'by_type')}
        if cfg.publish_if_complete and final_coverage['target_scope_complete']:
            publication = builder.publish(cfg.out, state, final_coverage, cfg.profiles, cfg.metadata.parent, args, cancel)
            archive_root, archive_manifest = cfg.out, publication
            if cfg.integrate_into:
                from domain_foldseek_ecod import merge
                target = cfg.integrate_into.expanduser().resolve()
                publication = merge([target, cfg.out], target, executable=args.foldseek, replace=True)
                archive_root, archive_manifest = target, publication
                builder.atomic_json(cfg.out / 'integration.json', {'library': str(target),
                    'reference_count': publication['reference_count'], 'pfam_coverage': final_coverage,
                    'manifest_sha256': builder.checksum(target / 'manifest.json', 'sha256')})
                if cfg.make_default:
                    builder.default_library(target, write=True)
            if cfg.archive:
                builder.distribution_archive(archive_root, archive_manifest, cfg.archive)
            report['publication'] = 'complete folded-domain scope published'
        else:
            report['publication'] = 'checkpoint only; missing in-scope families remain' if not final_coverage['target_scope_complete'] else 'complete checkpoint; publication not requested'
        builder.atomic_json(cfg.report, report)
    return report


if __name__ == '__main__':
    main()
