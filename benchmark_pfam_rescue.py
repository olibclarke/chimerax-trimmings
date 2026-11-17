#!/usr/bin/env python3
"""Compare local Pfam scanning with targeted HMMER sensitivity/domain retries.

Input JSON: {"cases": [{"id", "sequence", "candidate_families": ["PF..."],
"structural_regions": [[start, end]], "author_residues": [optional numbers]}]}.
Coordinates are one-based sequence positions. This experiment does not change
the domain command, reference library, or curated acceptance thresholds.
"""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from build_pfam_library import catalog


def command(arguments, folder, timeout):
    started = time.monotonic()
    with (folder / 'command.log').open('a') as log:
        log.write(json.dumps([str(a) for a in arguments]) + '\n')
        log.flush()
        subprocess.run(arguments, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=timeout)
    return time.monotonic() - started


def parse_table(path, families, queries):
    hits = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        row = line.split(maxsplit=22)
        accession = row[1].split('.')[0]  # hmmscan: target profile, query sequence.
        family = families[accession]
        query = queries[row[3]]
        if row[1] != family['version'] or int(row[2]) != family['length'] or int(row[5]) != len(query['sequence']):
            raise ValueError('Profile/query metadata mismatch in HMMER output.')
        start, end = int(row[17]), int(row[18])
        if not 1 <= start <= end <= len(query['sequence']):
            raise ValueError('Invalid sequence positions in HMMER output.')
        first = query['start'] - 1
        sequence_score, domain_score = float(row[7]), float(row[13])
        hits.append({'accession': accession, 'name': row[0], 'query': row[3],
                     'start': first + start, 'end': first + end,
                     'envelope': [first + int(row[19]), first + int(row[20])],
                     'hmm_start': int(row[15]), 'hmm_end': int(row[16]),
                     'profile_coverage': (int(row[16]) - int(row[15]) + 1) / family['length'],
                     'sequence_evalue': float(row[6]), 'sequence_score': sequence_score,
                     'conditional_evalue': float(row[11]), 'independent_evalue': float(row[12]),
                     'domain_score': domain_score, 'bias': float(row[14]),
                     'ga_sequence': family['ga'][0], 'ga_domain': family['ga'][1],
                     'curated': sequence_score >= family['ga'][0] and domain_score >= family['ga'][1],
                     'query_region': [query['start'], query['start'] + len(query['sequence']) - 1]})
    return hits


def placements(hits):
    """Deduplicate repeated searches, retaining distinct occurrences per family."""
    chosen = []
    for hit in sorted(hits, key=lambda h: (-h['domain_score'], h['independent_evalue'])):
        if any(hit['accession'] == old['accession'] and
               max(0, min(hit['end'], old['end']) - max(hit['start'], old['start']) + 1) >=
               .5 * min(hit['end'] - hit['start'] + 1, old['end'] - old['start'] + 1)
               for old in chosen):
            continue
        chosen.append(hit)
    return sorted(chosen, key=lambda h: (h['start'], h['end'], h['accession']))


def new_placements(hits, baseline):
    return [hit for hit in placements(hits) if not any(
        hit['accession'] == old['accession'] and
        max(0, min(hit['end'], old['end']) - max(hit['start'], old['start']) + 1) >=
        .5 * min(hit['end'] - hit['start'] + 1, old['end'] - old['start'] + 1)
        for old in baseline)]


def region_queries(case, baseline, padding=12):
    sequence = case['sequence']
    covered = {p for hit in baseline for p in range(hit['start'], hit['end'] + 1)}
    regions = set()
    start = None
    for p in range(1, len(sequence) + 2):
        if p <= len(sequence) and p not in covered:
            if start is None:
                start = p
        elif start is not None:
            if p - start >= 20:
                regions.add((max(1, start - padding), min(len(sequence), p - 1 + padding)))
            start = None
    for first, last in case.get('structural_regions', []):
        overlap = len(set(range(first, last + 1)) & covered)
        if overlap < .8 * (last - first + 1):
            regions.add((max(1, first - padding), min(len(sequence), last + padding)))
    return {f'region{i}': {'start': a, 'sequence': sequence[a - 1:b]} for i, (a, b) in enumerate(sorted(regions))}


def scan(database, queries, folder, args, models, maximum=False, curated=True):
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / 'queries.fasta'
    fasta.write_text(''.join(f'>{identifier}\n{query["sequence"]}\n' for identifier, query in queries.items()))
    if not queries:
        return {'seconds': 0, 'queries': 0, 'hits': []}
    table = folder / 'domains.tsv'
    options = ['--cut_ga'] if curated else ['-E', '10', '--domE', '10']
    if maximum:
        options.append('--max')
    seconds = command(['hmmscan', '--noali', '--cpu', str(args.threads), '-Z', str(models),
                       '--domtblout', str(table), *options, str(database), str(fasta)], folder, args.timeout)
    hits = parse_table(table, args.families, queries)
    if curated and not all(h['curated'] for h in hits):
        raise ValueError('A --cut_ga result does not meet its model gathering thresholds.')
    return {'seconds': seconds, 'queries': len(queries), 'hits': hits}


def benchmark_case(case, args):
    case_started = time.monotonic()
    folder = args.out / case['id']
    folder.mkdir(parents=True, exist_ok=True)
    full = {'full': {'start': 1, 'sequence': case['sequence']}}
    modes = {'full_database_standard': scan(args.profiles, full, folder / 'full_database_standard',
                                            args, len(args.families))}
    baseline = modes['full_database_standard']['hits']
    candidates = sorted(set(case['candidate_families']) | {h['accession'] for h in baseline})
    unavailable = sorted(set(candidates) - args.families.keys())
    candidates = sorted(set(candidates) & args.families.keys())
    if not candidates:
        raise ValueError('No available candidate profiles for ' + case['id'])
    keyfile = folder / 'candidate_keys.txt'
    keyfile.write_text(''.join(args.families[key]['version'] + '\n' for key in candidates))
    subset = folder / 'candidates.hmm'
    preparation = command(['hmmfetch', '-o', str(subset), '-f', str(args.profiles), str(keyfile)], folder, args.timeout)
    preparation += command(['hmmpress', '-f', str(subset)], folder, args.timeout)
    regions = region_queries(case, baseline)
    (folder / 'regions.json').write_text(json.dumps(regions, indent=2) + '\n')
    for name, queries, maximum in [('target_full_standard', full, False), ('target_full_max', full, True),
                                   ('target_regions_standard', regions, False), ('target_regions_max', regions, True)]:
        print(case['id'], name, len(queries), 'queries,', len(candidates), 'profiles', flush=True)
        modes[name] = scan(subset, queries, folder / name, args, len(args.families), maximum)
    # Exploratory reporting is kept separate from curated rescues, including its time.
    print(case['id'], 'exploratory below-cutoff reporting', flush=True)
    modes['target_full_max_exploratory'] = scan(subset, full, folder / 'target_full_max_exploratory',
                                                args, len(args.families), True, curated=False)
    modes['target_regions_max_exploratory'] = scan(subset, regions, folder / 'target_regions_max_exploratory',
                                                   args, len(args.families), True, curated=False)
    for result in modes.values():
        result['placements'] = placements(result['hits'])
        result['new_curated'] = new_placements([h for h in result['hits'] if h['curated']], baseline)
        result['new_tentative'] = new_placements([h for h in result['hits'] if not h['curated']
                                                 and h['independent_evalue'] <= .01], baseline)
    result = {'id': case['id'], 'length': len(case['sequence']), 'candidate_profiles': candidates,
              'unavailable_profiles': unavailable, 'profile_preparation_seconds': preparation,
              'elapsed_seconds': time.monotonic() - case_started, 'modes': modes}
    for h in [hit for mode in modes.values() for hit in mode['hits']]:
        if 'author_residues' in case:
            h['author_start'] = case['author_residues'][h['start'] - 1]
            h['author_end'] = case['author_residues'][h['end'] - 1]
    (folder / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print('COMPLETE', case['id'], {name: {'seconds': round(m['seconds'], 3),
                                         'curated': len([h for h in m['placements'] if h['curated']]),
                                         'new': len(m['new_curated'])} for name, m in modes.items()}, flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--profiles', required=True, type=Path)
    parser.add_argument('--metadata', required=True, type=Path, help='Matching Pfam-A.hmm.dat.gz')
    parser.add_argument('--models', type=int, required=True, help='Expected number of profiles in the complete release')
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=600, help='Per-command timeout, seconds')
    args = parser.parse_args()
    args.out, args.profiles = args.out.resolve(), args.profiles.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    args.families = catalog(args.metadata, args.models)
    cases = json.loads(args.input.read_text())['cases']
    for tool in ('hmmscan', 'hmmfetch', 'hmmpress'):
        if shutil.which(tool) is None:
            raise ValueError('Missing executable: ' + tool)
    for case in cases:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', case['id']) or not re.fullmatch(r'[A-Z]+', case['sequence']):
            raise ValueError('Invalid benchmark case identifier/sequence.')
    report = {'profiles': str(args.profiles), 'models': args.models, 'threads': args.threads,
              'hmmer': subprocess.check_output(['hmmscan', '-h'], text=True).splitlines()[:3],
              'notes': ['Curated sequence and domain GA thresholds are unchanged.',
                        'Subset E-values use the full profile count (-Z).',
                        'Local hits precede Pfam clan competition; overlapping families are not independent domains.',
                        'Candidate families and structural regions are provided in the input, independently of --max results.',
                        'Region queries use uncovered stretches >=20 residues and unresolved structural matches, with 12 residues of context.'],
              'cases': []}
    for case in cases:
        print('START', case['id'], len(case['sequence']), 'residues', flush=True)
        report['cases'].append(benchmark_case(case, args))
        (args.out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
