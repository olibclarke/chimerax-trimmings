#!/usr/bin/env python3
"""Add audited Pfam labels to ECOD metadata without changing the search database.

Reuse upstream ECOD family mappings and optional saved gathering-threshold hits
from build_pfam_library. Unresolved entries are explicit, not forced assignments.
"""

import argparse
from collections import Counter
from contextlib import closing
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3

import build_pfam_library as builder
from domain_foldseek_ecod import staged_generation


def natural_pieces(record):
    """Return offsets guaranteed to belong to one natural sequence piece.

    For experimental structures, missing coordinates can shorten the exported
    sequence. Keep only positions inside a piece for every possible placement
    of those missing residues; do not guess where an artificial join falls.
    """
    text = record.get('ecod_seqid_range', '')
    if not isinstance(text, str) or not re.fullmatch(r'(?:[A-Za-z0-9_.-]+:)?\d+-\d+(?:,(?:[A-Za-z0-9_.-]+:)?\d+-\d+)*', text):
        return []
    lengths = []
    for piece in text.split(','):
        first, last = map(int, piece.rsplit(':', 1)[-1].split('-'))
        if not 1 <= first <= last:
            return []
        lengths.append(last - first + 1)
    missing = sum(lengths) - record['length']
    if missing < 0 or (missing and record.get('structure_source') != 'pdb'):
        return []
    offset, pieces = 0, []
    for length in lengths:
        if offset + 1 <= min(offset + length - missing, record['length']):
            pieces.append((offset + 1, min(offset + length - missing, record['length'])))
        offset += length
    return pieces


def scan_reference_sequences(library, profiles, metadata, identifiers, work, executable='hmmsearch', threads=4):
    """Verify selected actual database sequences against the complete Pfam release."""
    import subprocess
    import domain_pfam

    if threads < 1:
        raise ValueError('Reference scan thread count must be positive.')
    work.mkdir(parents=True, exist_ok=True)
    requested = set(identifiers)
    keys = []
    with Path(str(library.database) + '.lookup').open() as stream:
        found = set()
        for line in stream:
            row = line.split()
            if row[1] in requested:
                keys.append(row[0]); found.add(row[1])
    if found != requested or not requested:
        raise ValueError('Requested reference identifiers must exist in the search database.')
    (work / 'keys.txt').write_text(''.join(key + '\n' for key in keys))
    builder.say(f'Verifying {len(requested)} reference sequences against Pfam gathering thresholds')
    with (work / 'scan.log').open('w') as log:
        for command in (["foldseek", 'createsubdb', work / 'keys.txt', library.database, work / 'subset'],
                        ["foldseek", 'convert2fasta', work / 'subset', work / 'references.fasta']):
            subprocess.run(list(map(str, command)), check=True, stdout=log, stderr=subprocess.STDOUT)
    sequences = {}
    for identifier, sequence in builder.fasta_records(work / 'references.fasta'):
        identifier = identifier if identifier in requested else identifier.removesuffix('.pdb')
        if identifier not in requested or identifier in sequences or not sequence or not sequence.isascii() or not sequence.isalpha():
            raise ValueError('Exported reference sequences do not match the requested identifiers.')
        record = library.record(identifier)
        if len(sequence) != record.get('length', len(record.get('sequence', ''))) or (
                record.get('sequence') is not None and sequence != record['sequence']):
            raise ValueError('Exported reference sequence disagrees with its metadata.')
        sequences[identifier] = sequence
    if set(sequences) != requested:
        raise ValueError('Foldseek omitted requested reference sequences.')
    with (work / 'scan.log').open('a') as log:
        subprocess.run(list(map(str, [executable, '--cpu', threads, '--cut_ga', '--noali', '--domtblout', work / 'hits.domtbl',
                                     '-o', work / 'hmmer.txt', profiles, work / 'references.fasta'])),
                       check=True, stdout=log, stderr=subprocess.STDOUT)
    families = domain_pfam.metadata({'metadata': metadata})
    mappings, counts = {}, Counter()
    for family, target, first, last, score, evalue, coverage in builder.parse_domtbl(work / 'hits.domtbl', families, .8):
        if target not in requested:
            raise ValueError('Unexpected reference in sequence scan.')
        f, record = families[family], library.record(target)
        if score < f['ga'][1]:
            raise ValueError('Reference hit is below its Pfam gathering threshold.')
        if f['type'] not in builder.FOLDED_PROFILE_TYPES or record.get('catalog') != 'ECOD':
            continue
        if len(sequences[target]) != record['length']:
            raise ValueError('Reference sequence length changed during scan.')
        if last - first + 1 < 20 or not any(a <= first <= last <= b for a, b in natural_pieces(record)):
            counts['unverifiable_natural_span'] += 1
            continue
        mappings.setdefault(target, []).append({'accession': family, 'name': f['name'], 'description': f['description'],
            'method': 'hmmer', 'segments': [[first, last]], 'bitscore': score, 'evalue': evalue,
            'profile_coverage': coverage, 'profile_version': f['version'], 'clan': f['clan'],
            'evidence': 'Pfam gathering-threshold reference-sequence match'})
        counts['accepted_hits'] += 1
    audit = {'scope': 'Explicit reference IDs, full Pfam profile collection, gathering thresholds',
             'library_manifest_sha256': library.digest, 'metadata_sha256': builder.checksum(metadata, 'sha256'),
             'profiles_sha256': builder.checksum(profiles, 'sha256'),
             'references': {key: hashlib.sha256(sequence.encode()).hexdigest() for key, sequence in sequences.items()},
             'counts': dict(counts), 'mappings': mappings}
    builder.atomic_json(work / 'verified-mappings.json', audit)
    return mappings, audit


def sequence_mappings(checkpoint, library, families, metadata):
    """Import only provenance-matched hits; reject artificial piece-spanning hits."""
    if checkpoint is None:
        return {}, Counter()
    provenance = json.loads((checkpoint.parent / 'provenance.json').read_text())
    accepted = {library.digest, library.manifest.get('pfam_mapping', {}).get('base_manifest_sha256')}
    sources = {i for i, source in enumerate(provenance['reuse'])
               if Path(source['path']).resolve() == library.folder and source['manifest_sha256'] in accepted}
    if not sources or provenance['upstream_md5']['Pfam-A.hmm.dat.gz'] != builder.checksum(metadata, 'md5'):
        raise ValueError('Checkpoint library/profile provenance does not match this mapping request.')
    result, counts = {}, Counter()
    with closing(sqlite3.connect('file:' + str(checkpoint) + '?mode=ro', uri=True)) as db:
        query = '''SELECT h.family,t.original,t.sequence,h.first,h.last,h.score,h.evalue,h.coverage,t.library
                   FROM hits h JOIN targets t ON t.id=h.target'''
        with closing(sqlite3.connect('file:' + str(library.metadata) + '?mode=ro', uri=True)) as refs:
            for family, target, sequence, first, last, score, evalue, coverage, source in db.execute(query):
                if source not in sources or family not in families:
                    continue
                f = families[family]
                if f['type'] not in builder.FOLDED_PROFILE_TYPES:
                    continue
                row = refs.execute('SELECT record FROM refs WHERE id=?', (target,)).fetchone()
                if row is None:
                    raise ValueError('Checkpoint target is missing from the source library: ' + target)
                record = json.loads(row[0])
                if record.get('catalog') != 'ECOD':
                    continue
                if len(sequence) != record['length'] or not 1 <= first <= last <= len(sequence):
                    raise ValueError('Checkpoint reference length or coordinates changed.')
                if not (all(math.isfinite(v) for v in (score, evalue, coverage))
                        and score >= f['ga'][1] and 0 <= evalue and .8 <= coverage <= 1):
                    raise ValueError('Checkpoint hit does not meet the mapping criteria.')
                if last - first + 1 < 20:
                    counts['short_sequence_match'] += 1
                    continue
                if not any(a <= first <= last <= b for a, b in natural_pieces(record)):
                    counts['unverifiable_natural_span'] += 1
                    continue
                result.setdefault(target, []).append({'accession': family, 'name': f['name'],
                    'description': f['description'], 'method': 'hmmer', 'segments': [[first, last]],
                    'bitscore': score, 'evalue': evalue, 'profile_coverage': coverage,
                    'profile_version': f['version'], 'evidence': 'Pfam gathering-threshold reference-sequence match'})
                counts['imported_hits'] += 1
    return result, counts


def annotate_reference(record, families, matches=()):
    record = dict(record)
    if record.get('catalog') != 'ECOD':
        return record
    upstream = record.get('ecod_pfam_accessions', record.get('pfam_accessions', []))
    upstream = list(dict.fromkeys(upstream))
    valid = [a for a in upstream if a in families]
    previous = [m for m in record.get('pfam_mappings', []) if m.get('method') == 'hmmer'
                and m.get('accession') in families and m.get('profile_version') == families[m['accession']]['version']
                and sum(b - a + 1 for a, b in m.get('segments', [])) >= 20]
    matches = list({json.dumps(m, sort_keys=True): m for m in [*previous, *matches]}.values())
    mappings = [{'accession': a, 'name': families[a]['name'], 'description': families[a]['description'],
                 'method': 'ecod', 'evidence': 'Upstream ECOD family-to-Pfam cross-reference'} for a in valid]
    mappings.extend(matches)
    hit_families = {m['accession'] for m in matches}
    if hit_families - set(valid) and valid:
        status = 'conflicting'
    elif len(hit_families) == 1:
        status = 'sequence_supported'
    elif len(hit_families) > 1 or len(valid) > 1:
        status = 'ambiguous'
    elif len(valid) == 1:
        status = 'upstream_single'
    else:
        status = 'unresolved'
    record.update(ecod_pfam_accessions=upstream, pfam_mappings=mappings, pfam_mapping_status=status,
                  pfam_accessions=list(dict.fromkeys([*valid, *(m['accession'] for m in matches)])))
    return record


def publish(library, metadata, models, checkpoint=None, scan_ids=None, profiles=None, scan_work=None, threads=4):
    families = builder.catalog(metadata, models)
    library = builder.Library(library)
    if library.metadata is None:
        raise ValueError('An indexed schema-2 library is required.')
    matches, import_counts = sequence_mappings(checkpoint, library, families, metadata)
    scan_audit = None
    if scan_ids is not None:
        if profiles is None or scan_work is None:
            raise ValueError('Reference scanning requires --profiles and --scan-work.')
        fresh, scan_audit = scan_reference_sequences(library, profiles, metadata, scan_ids, scan_work, threads=threads)
        for target in scan_ids:
            # A full rescan supersedes saved matches for the scanned references.
            matches[target] = fresh.get(target, [])
    with builder.build_lock(library.folder), staged_generation(library.folder) as generation:
        if hashlib.sha256((library.folder / 'manifest.json').read_bytes()).hexdigest() != library.digest:
            raise ValueError('The source library changed during mapping preparation.')
        destination = generation / 'references.sqlite'
        counts, pfams, added, number = Counter(), set(), 0, 0
        with closing(sqlite3.connect(destination)) as db, (generation / 'pfam-mapping.tsv').open('w') as stream:
            db.execute('CREATE TABLE refs(id TEXT PRIMARY KEY,record TEXT)')
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['target', 'ecod_family', 'mapping_status', 'pfam_accessions', 'sequence_supported'])
            for number, source in enumerate(library.records(), 1):
                if scan_ids is not None and source['id'] in scan_ids:
                    source = {**source, 'pfam_mappings': [m for m in source.get('pfam_mappings', []) if m.get('method') != 'hmmer']}
                record = annotate_reference(source, families, matches.get(source['id'], ()))
                db.execute('INSERT INTO refs VALUES(?,?)', (record['id'], json.dumps(record)))
                if record.get('catalog') == 'ECOD':
                    counts[record['pfam_mapping_status']] += 1
                    pfams.update(record['pfam_accessions'])
                    added += bool(record['pfam_accessions']) and not bool(source.get('ecod_pfam_accessions', source.get('pfam_accessions')))
                    writer.writerow([record['id'], record['accession'], record['pfam_mapping_status'],
                                     ','.join(record['pfam_accessions']), bool(matches.get(record['id']))])
                if number % 25000 == 0:
                    db.commit()
                    builder.say(f'Mapped metadata: {number}/{library.manifest["reference_count"]} references')
            db.commit()
            if number != library.manifest['reference_count']:
                raise ValueError('Metadata reference count changed.')
        audit = {'scope': 'All ECOD references; upstream mappings plus reusable validated sequence hits',
                 'ecod_references': sum(counts.values()), 'status_counts': dict(counts),
                 'current_pfam_families': len(pfams), 'newly_mapped_references': added,
                 'checkpoint_import': dict(import_counts), 'complete_sequence_rescan': False,
                 'targeted_sequence_scan': ({k: v for k, v in scan_audit.items() if k != 'mappings'} if scan_audit else None),
                 'notes': ['Unresolved entries have no accepted mapping in the available evidence, not proven absence of a Pfam match.',
                           'Additional sequence mappings require Domain/Family/Repeat type and at least 20 reference residues.',
                           'Multiple or conflicting mappings remain explicit.',
                           'Coordinates, Foldseek structural alphabet, target identifiers and search cutoffs are unchanged.']}
        builder.atomic_json(generation / 'pfam-mapping.json', audit)
        builder.atomic_json(generation / 'source-manifest.json', library.manifest)
        manifest = {**library.manifest, 'reference_metadata': str(destination.relative_to(library.folder)),
                    'reference_metadata_sha256': builder.checksum(destination, 'sha256'),
                    'pfam_mapping': {'schema': 1, 'base_manifest_sha256': library.manifest.get('pfam_mapping', {}).get('base_manifest_sha256', library.digest),
                                     'audit': str((generation / 'pfam-mapping.json').relative_to(library.folder)),
                                     'families': str((generation / 'pfam-mapping.tsv').relative_to(library.folder)),
                                     'summary': audit}}
        builder.atomic_json(library.folder / 'manifest.json', manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', required=True, type=Path)
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--models', required=True, type=int)
    parser.add_argument('--checkpoint', type=Path, help='Saved state.sqlite with matching library/profile provenance')
    parser.add_argument('--scan-ids', type=Path, help='Reference ID list to verify against the full profile collection')
    parser.add_argument('--profiles', type=Path, help='Full Pfam-A.hmm, for optional reference scanning')
    parser.add_argument('--scan-work', type=Path, help='Keep scan evidence and provenance here')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    ids = set(args.scan_ids.read_text().split()) if args.scan_ids else None
    result = publish(args.library, args.metadata, args.models, args.checkpoint, ids, args.profiles, args.scan_work, args.threads)
    print(json.dumps(result['pfam_mapping']['summary'], indent=2))


if __name__ == '__main__':
    main()
