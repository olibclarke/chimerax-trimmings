#!/usr/bin/env python3
"""Run a breadth-first Pfam template sweep against an existing checkpoint.

Each visited family gets one seed candidate, chosen for coverage and sequence
completeness. Every coordinate still passes the original builder's current-profile,
coverage, sequence and pLDDT checks. Unsuccessful families are deferred, never
declared structurally impossible. Accepted templates enter the resumable checkpoint;
the published search library is untouched. No exhaustive member crawl is performed.
"""

from collections import Counter
from contextlib import closing
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import sqlite3
from threading import Event
import time

import build_pfam_library as builder


def ranked_seeds(rows, family):
    # Prefer intact, near-profile-length seed fragments without inventing a
    # family assignment from protein names or borrowing another family's hit.
    return sorted(rows, key=lambda row: (row[2].count('X') / max(1, len(row[2])),
                  abs(len(row[2]) - family['length']) / max(1, family['length']), row[0]))


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path, help='Existing Pfam builder output directory')
    parser.add_argument('--work', required=True, type=Path, help='Existing work generation with state.sqlite')
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--profiles', required=True, type=Path)
    parser.add_argument('--models', required=True, type=int)
    parser.add_argument('--sample', type=int, default=40, help='Sample folded-family candidates; 0 visits all missing families')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--report', required=True, type=Path)
    cfg = parser.parse_args(argv)
    if cfg.sample < 0 or cfg.workers < 1:
        parser.error('Sample must be nonnegative and workers must be positive.')
    cfg.out, cfg.work, cfg.report = cfg.out.resolve(), cfg.work.resolve(), cfg.report.resolve()
    families = builder.catalog(cfg.metadata, cfg.models)
    state_path = cfg.work / 'state.sqlite'
    if not state_path.is_file():
        raise ValueError('The work directory must contain an existing prepared state.sqlite checkpoint.')
    cancel = Event()
    args = builder.parser().parse_args(['--out', str(cfg.out)])
    args.no_interpro = True
    args._prediction_gate = None
    cache = builder.PublicCache(cfg.out / 'downloads/structures', cancel)
    started = time.monotonic()
    with builder.build_lock(cfg.out), closing(builder.state_database(state_path)) as state:
        represented = builder.represented(state)
        pending = [a for a, f in families.items() if a not in represented]
        # Motifs, very short profiles and disordered entries remain in the audit;
        # they are not used to estimate throughput for folded-domain templates.
        eligible = [a for a in pending if families[a]['type'] in ('Domain', 'Family', 'Repeat')
                    and families[a]['length'] >= args.min_residues]
        eligible_set = set(eligible)
        random.Random(cfg.seed).shuffle(eligible)
        selected = (eligible + sorted(set(pending) - eligible_set)) if cfg.sample == 0 else eligible[:cfg.sample]
        pilot = cfg.report.parent
        pilot.mkdir(parents=True, exist_ok=True)
        pilot_state_path = pilot / 'candidate-state.sqlite'
        ledger = sqlite3.connect(cfg.out / 'seed-sweeps.sqlite')
        ledger.execute('CREATE TABLE IF NOT EXISTS attempts(family TEXT,uniprot TEXT,seed_sha256 TEXT,outcome TEXT,'
                       'PRIMARY KEY(family,uniprot,seed_sha256))')
        used = set(ledger.execute('SELECT family,uniprot,seed_sha256 FROM attempts'))
        scheduled, seed_keys = [], {}
        with closing(builder.state_database(pilot_state_path)) as candidates:
            candidates.execute('DELETE FROM seeds')
            for number, accession in enumerate(selected, 1):
                rows = state.execute('SELECT uniprot,header,sequence FROM seeds WHERE family=?', (accession,)).fetchall()
                seeds = [row for row in ranked_seeds(rows, families[accession])
                         if (accession, row[0], hashlib.sha256(row[2].encode()).hexdigest()) not in used]
                if seeds:
                    uniprot, header, sequence = seeds[0]
                    candidates.execute('INSERT OR IGNORE INTO seeds VALUES(?,?,?,?)',
                                       (accession, uniprot, header, sequence))
                    scheduled.append(accession)
                    seed_keys[accession] = accession, uniprot, hashlib.sha256(sequence.encode()).hexdigest()
                if number % 1000 == 0:
                    builder.say(f'Preparing sweep: {number}/{len(selected)} families; {len(scheduled)} candidates')
            candidates.commit()
        report = {'started': datetime.now(timezone.utc).isoformat(), 'strategy': 'one ranked seed per family per pass',
                  'seed': cfg.seed, 'pending_families': len(pending), 'eligible_for_sample': len(eligible),
                  'scope': 'all missing families' if cfg.sample == 0 else 'folded-family sample',
                  'scheduled': len(scheduled), 'without_unused_seed': len(selected) - len(scheduled),
                  'low_priority_by_type_or_size': dict(Counter(families[a]['type'] for a in pending if a not in eligible_set)),
                  'criteria': {'minimum_profile_and_modeled_coverage': args.min_coverage,
                               'minimum_mean_plddt': args.min_plddt, 'minimum_residues': args.min_residues},
                  'completed': 0, 'verified': 0, 'checkpoint_families': len(represented), 'request_starts': 0,
                  'families': []}
        builder.atomic_json(cfg.report, report)
        def publish_progress():
            summary = {k:v for k,v in report.items() if k != 'families'}
            builder.atomic_json(pilot / 'progress.json', summary)
            builder.atomic_json(cfg.out / 'progress.json', {
                'phase': 'seed sweep' if cfg.sample == 0 else 'seed sweep sample',
                'completed': report['completed'], 'pending': len(scheduled),
                'covered_families': report['checkpoint_families'], 'total_families': len(families),
                'verified_this_sweep': report['verified'], 'scope': report['scope'],
                'without_unused_seed': report['without_unused_seed']})
        publish_progress()
        last_snapshot = time.monotonic()
        events = pilot / 'families.jsonl'
        def fill(accession):
            try:
                return builder.fill_family(accession, pilot_state_path, families, cfg.profiles,
                                           cache, cfg.work, args, cancel)
            except (OSError, RuntimeError, ValueError) as error:
                return [], [{'source': 'pilot', 'error': str(error)}]
        try:
            for accession, (records, errors) in builder.parallel_jobs(scheduled, fill, cfg.workers, cancel):
                accepted = []
                for record in records:
                    if builder.add_reference(state, record):
                        accepted.append(record[0]['id'])
                state.commit()
                status = 'verified' if accepted else 'deferred; one candidate only'
                if not any(error.get('retryable') for error in errors):
                    ledger.execute('INSERT OR REPLACE INTO attempts VALUES(?,?,?,?)', (*seed_keys[accession], status))
                    ledger.commit()
                detail = {'accession': accession, 'name': families[accession]['name'],
                          'candidate': {'uniprot': seed_keys[accession][1], 'seed_sha256': seed_keys[accession][2]},
                          'status': status, 'references': accepted, 'errors': errors}
                report['families'].append(detail)
                with events.open('a') as stream:
                    stream.write(json.dumps(detail) + '\n')
                report.update(elapsed_seconds=time.monotonic() - started,
                              verified=report['verified'] + bool(accepted),
                              completed=len(report['families']), request_starts=cache.request_starts,
                              checkpoint_families=report['checkpoint_families'] + bool(accepted))
                publish_progress()
                if len(scheduled) <= 100 or report['completed'] % 100 == 0 or time.monotonic() - last_snapshot >= 30:
                    builder.atomic_json(cfg.report, report)
                    last_snapshot = time.monotonic()
                builder.say(f"Sweep {report['completed']}/{len(scheduled)}: {accession} {status}; "
                            f"{report['verified']} verified; {report['elapsed_seconds']:.1f}s")
        except BaseException:
            cancel.set()
            raise
        finally:
            builder.atomic_json(cfg.report, report)
            ledger.close()
            builder.coverage_report(cfg.out, state, families, set(families), args)
        report['finished'] = datetime.now(timezone.utc).isoformat()
        builder.atomic_json(cfg.report, report)
    print(json.dumps({k:v for k,v in report.items() if k != 'families'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
