# Build a full Pfam reference bundle

`build_pfam_library.py` performs the Pfam-first build in one invocation. It
includes **every profile in the pinned Pfam release**, reannotates existing
structural references, and obtains missing representatives from public sources.
It resumes completed work and checks completeness before publication.

Generated archives include the MIT software license, [data attribution](DATA_NOTICE.md)
and [third-party notices/citations](THIRD_PARTY.md), alongside provenance and coverage
reports. Source data retains its upstream terms. PROSITE caches and predictor
environments are not included. To use the existing local library, start with
[the quickstart](README.md#local-domain-library).

Keep it beside `domain_foldseek.py`, `domain_foldseek_catalog.py`, and
`domain_foldseek_ecod.py`. This runs in an ordinary terminal, outside ChimeraX.

## Install dependencies

Use Python 3.11 or later with Gemmi, HMMER 3, Foldseek, and curl on PATH. On macOS:

```bash
brew install hmmer foldseek
python3 -m pip install gemmi
```

Use a virtual environment for Gemmi if your Python installation requires one.

## One-command build

For a new build restricted to folded-domain candidates, add `--scope folded`.
This targets Pfam **Domain, Family and Repeat** types and excludes explicit
**Motif, Disordered and Coiled-coil** types from the structural library.
The complete sequence-profile bundle still includes every Pfam entry. For an
existing checkpoint, use the [InterPro completion pass](#interpro-guided-folded-domain-completion)
below to reuse all previously verified templates and avoid rescanning.

From the Trimmings directory:

```bash
python3 build_pfam_library.py \
  --release 39.0 \
  --reuse .domain_data/ecod-cdd \
  --out .domain_data/pfam-complete \
  --threads 8 --workers 4 \
  --integrate-into .domain_data/ecod-cdd --make-default \
  --archive .domain_data/pfam-39.0.tar.gz
```

`--reuse` accepts an existing Trimmings Foldseek library and can be repeated.
It is optional: without it the builder obtains representatives directly from
PDB/AlphaFold. Choose a separate output directory. `--integrate-into` merges the
published result into an existing library only after verification, retaining its
previous manifest and database generation. `--make-default` selects that merged
library for ChimeraX searches without an explicit `library` argument. Omit these
two options to build a standalone library. The archive contains the merged
library when integration is requested.

The script:

1. Downloads the pinned release's checksums, version, complete family catalog,
   HMM library, and seed alignments. Checksums and inventory/model versions must
   agree. Pfam 39.0 has 53,275 families. The HMM and seed downloads are roughly
   628 MB and 224 MB compressed; uncompressed profiles, indexes, coordinates,
   and checkpoints require additional space. Large downloads use resumable,
   parallel byte ranges when the server supports them; the assembled files must
   still pass the release checksums.
2. Prepares the **full HMMER database**, including `hmmpress` and `hmmfetch`
   indexes. An explicit pilot subset never removes profiles from this bundle.
3. Exports existing Foldseek reference sequences in bulk. HMMER applies the
   pinned Pfam gathering thresholds, starting with a diverse batch of 10,000
   references, then doubling the disjoint reference batches. Families with enough
   verified coordinates leave subsequent scans; unresolved families are searched
   against every remaining reference. By default two scan processes share the
   CPU budget. Completed shards of 250 profiles are checkpoints, so interruption
   does not discard the entire scan.
4. Extracts only needed coordinate targets from existing libraries in batches,
   crops to validated Pfam matches, and writes compact C-alpha references.
   Existing ECOD family names alone do not establish new Pfam memberships.
5. Indexes PDBe's bulk SIFTS Pfam-to-PDB mappings. For remaining families it
   tries PDB chains, all indexed seed representatives in AlphaFold, and then
   paginated InterPro protein matches. Known proteins absent from AlphaFold
   are skipped. All fetched crops are checked against the pinned Pfam profile,
   polymer sequence, modeled coverage, and prediction confidence as applicable.
   There is no default family/candidate sampling limit.
6. Optionally predicts remaining natural seed sequences using a supplied local
   predictor command. It validates those models too.
7. Writes a per-family report, checks that every family has a valid structural
   representative, and atomically publishes the Foldseek library. If requested,
   it also packages the published database, complete profiles, provenance, and
   reports in a `.tar.gz` with a SHA256 file. Downloads and scratch data are
   omitted from this distributable archive.

The default is **one distinct reference sequence per family**, prioritizing
family coverage. `--representatives 3` requests additional diversity; full
family coverage means at least one valid representative for each family.

For an existing checkpoint, the [breadth-first seed sweep](#accelerate-an-existing-checkpoint)
is the recommended next step before an exhaustive public-member crawl.

The terminal reports completed scan shards, coordinate-verification batches,
family counts, elapsed time, and periodic activity while waiting. `progress.json`
records the current phase and counters. Coverage reports refresh after each
reference batch and at completion or interruption.

## Accelerate an existing checkpoint

An exhaustive family-by-family crawl can spend hours trying thousands of
candidates for a short motif or disordered family while ordinary domains wait.
`sweep_pfam_library.py` visits one ranked seed per missing family per pass:

```bash
python3 sweep_pfam_library.py \
  --out .domain_data/pfam-complete \
  --work /path/to/existing/work/generation \
  --metadata .domain_data/pfam-audit/v39.0/Pfam-A.hmm.dat.gz \
  --profiles .domain_data/pfam-complete/profiles/Pfam39.0/Pfam-A.hmm \
  --models 53275 --sample 0 --workers 4 \
  --report .domain_data/pfam-complete/seed-sweep/report.json
```

Use the existing work-generation directory containing `state.sqlite`. Stop the
old builder with Ctrl-C first; verified templates and downloads remain saved.
The sweep shares its output lock and adds accepted references to that exact
checkpoint. It does not rebuild the HMM library, rescan existing references, or
publish an incomplete structural library.

`--sample 0` visits **all missing families**, with ordinary folded-domain
candidates first and motifs, short profiles, coiled coils and disordered families
later. Each visit gets one seed ranked by completeness and similarity in length
to its profile. Current-profile gathering thresholds, 80% profile/modeled coverage,
20-residue minimum and mean predicted pLDDT of 70 remain unchanged. Successes are
checkpointed immediately; failures are **deferred**, not declared absent.

`seed-sweeps.sqlite` remembers attempted seed candidates. Rerunning the sweep
advances to an unused seed for each remaining family. Temporary source failures
remain retryable. Persistent gaps may need PDB candidates, other members or a
local predictor; exhausting seeds is not proof that a family lacks a structure.

Public downloads share a per-host request interval and rate-limit cooldown across
workers. Numeric or HTTP-date `Retry-After` values are respected, including waits
longer than a minute. A throttled family is deferred instead of exhausting all
its remaining candidates against the same unavailable source.

The top-level `progress.json` identifies the active sweep and verified family
count. Detailed results are appended to `seed-sweep/families.jsonl`; compact
progress updates accompany periodic report snapshots. Completing a sweep does
**not** mean structural coverage is complete. The original completeness gate and
coverage audit remain in force before final publication/integration.

In the initial 40-family sample, this strategy added **28 verified templates in
79 seconds**, with unchanged quality criteria. The sample covered ordinary
folded-family candidates and is not an ETA for difficult or unstructured entries.
The first full sweep was started over the remaining **33,609 families**. Its
throughput estimates describe that pass, not guaranteed 100% structural coverage.

## Resume and fill unresolved families
## InterPro-guided folded-domain completion

`fill_pfam_interpro.py` continues an existing checkpoint using InterPro's
experimental-structure links and proteins known to have AlphaFold models.
It discovers candidate sequences and annotated domain boundaries, caches a
persistent inventory, and validates the annotated occurrence with the pinned
Pfam HMM before fetching coordinates. A protein downloaded for several families
is reused through the shared cache; coordinate parsing has a bounded cache.

```bash
python3 fill_pfam_interpro.py \
  --out .domain_data/pfam-complete \
  --work /path/to/existing/work/generation \
  --metadata .domain_data/pfam-audit/v39.0/Pfam-A.hmm.dat.gz \
  --profiles .domain_data/pfam-complete/profiles/Pfam39.0/Pfam-A.hmm \
  --models 53275 --sample 0 --rounds 3 --workers 4 \
  --report .domain_data/pfam-complete/interpro-fill/report.json \
  --publish-if-complete \
  --integrate-into .domain_data/ecod-cdd --make-default \
  --archive .domain_data/pfam-39.0-folded.tar.gz
```

Each round tries one ranked candidate for every still-missing family. Experimental
structures are preferred initially; the next round gives available AlphaFold
models priority, so several incomplete PDB chains cannot consume every attempt.
Candidates favor complete annotated spans and InterPro representatives. New
pages are fetched only when the current inventory is exhausted. Rerunning uses
the persistent attempt ledger to advance to unused candidates; transient source
failures remain retryable. `--sample 40` runs a reproducible pilot instead of all
missing families. The validated initial pilots recovered **23 of 42** previously
uncovered families without lowering acceptance thresholds.

The model's polymer sequence must match the InterPro sequence exactly. Validation
is confined to the specific annotated occurrence, with 12 residues of context.
For discontinuous annotations, the HMM sees the natural intervening sequence,
while exported coordinates preserve the separately annotated pieces. Each
template records the InterPro source URL, boundaries, sequence/coordinate
checksums, current-profile score and coverage, and predicted-domain confidence.

`interpro-candidates.sqlite` stores inventory pages, candidates, attempts and
pagination state. Top-level `progress.json` and `interpro-fill/report.json`
report visits, successes and categorized failures: missing models, sequence
mismatches, profile/coverage failures, low confidence and insufficient modeled
coverage. `interpro-fill/attempts.jsonl` records each completed visit. Earlier
references and the published search library remain usable during the run.

`folded-coverage.json` and `folded-families.tsv` separately audit the structural
target. In Pfam 39.0 it comprises **52,484 Domain/Family/Repeat entries**; **791**
explicit motif/disorder/coiled-coil entries are outside this structural scope.
These types are a practical candidate scope, not proof that every remaining
Family/Repeat entry folds independently. A rejected or low-confidence candidate
stays unresolved; it is never silently removed from the denominator. Short
profiles also remain eligible if a natural validated segment supplies enough
modeled residues. The >=20-residue limit and all confidence/coverage criteria
remain unchanged.

Completeness fields distinguish the full Pfam inventory from the folded target:
`structural_coverage_complete` concerns all Pfam entries, while
`target_scope_complete`/`folded_scope_complete` concerns the explicit folded scope.
A pass finishing does not mean the library is complete. `--publish-if-complete`
publishes only after every in-scope family has verified coordinates; integration
and packaging happen only afterward. Explicit out-of-scope templates already in
the checkpoint are omitted from the published folded structural database, while
the full HMM bundle remains included for sequence annotation. The unclustered
reference collection and accession mapping should be retained when benchmarking
a subsequent representative-first clustered search.

Known conserved motifs and disorder remain sequence features. Existing Pfam
sequence searches retain their corresponding types; this builder does not add
a new low-complexity predictor or imply that low-confidence AlphaFold regions
are experimentally established disorder.

## Exhaustive fallback and local prediction

Run the same command again to resume. Completed scan shards, accepted references,
and successful downloads are reused. Public 404 responses are cached for one day;
temporary errors are retried. Failed families are attempted again on the next run.
An unchanged request that has already integrated successfully returns its existing
result, without scanning the merged library or adding duplicate references. A
missing distribution archive can be regenerated by rerunning the same command.

An optional predictor adapter receives a FASTA path and an output directory:

```bash
python3 build_pfam_library.py \
  --release 39.0 --reuse .domain_data/ecod-cdd \
  --out .domain_data/pfam-complete --threads 8 --workers 4 \
  --predict-command '/path/to/predictor {fasta} {outdir}' \
  --archive .domain_data/pfam-39.0.tar.gz
```

The adapter must write PDB or mmCIF files beneath `{outdir}`, preserve the seed's
amino-acid sequence, and store pLDDT in the C-alpha B-factor field. `{accession}`
is also available. Arguments are substituted without running a shell. For an
installed ColabFold environment, an adapter command can be:

```text
/path/to/colabfold_batch --num-models 1 --num-recycle 3 {fasta} {outdir}
```

Prediction needs the predictor's own models and suitable compute. Local
predictions run **one at a time** by default to avoid competing for GPU memory;
`--prediction-workers` can change this. Adding a predictor or enabling the
InterPro fallback reuses existing accepted structures and scan checkpoints.

Default acceptance criteria are a curated Pfam gathering-threshold match,
at least 80% profile and modeled-span coverage, 20–9,999 modeled C-alpha atoms,
and mean domain pLDDT of at least 70 for fetched or locally generated predictions.
Reused curated library coordinates retain their upstream provenance; the
compressed Foldseek database does not preserve residue pLDDT for remeasurement.
Changing acceptance criteria starts a separate work generation.

## What “complete” means

`build-coverage.json` separately records full **profile** coverage and verified
**structural** family coverage. `build-families.tsv` lists every family, its
reference count, and unresolved/error status. Detailed attempts are retained
under `work/*/families/PF*/attempts.jsonl`.

By default, **any missing structural family causes exit status 2**, and no new
Foldseek manifest or archive is published. A previously published library stays
usable. The verified full profile database remains available under
`profiles/Pfam39.0/`, even when structural coverage is incomplete.

## Use the published Pfam Foldseek library

The verified checkpoint published in this workspace on 16 November 2026 contains
**46,654 structural references**, covering **46,654/52,484** folded-scope Pfam
families (88.9%). It is a usable partial library; the 5,830 remaining gaps stay
listed in its coverage audit. Publishing this snapshot does not resume candidate
discovery or change template validation thresholds. Full sequence profiles remain
available for all 53,275 families.

```chimerax
domains #1 source foldseek library "/path/to/chimerax-trimmings/.domain_data/pfam-complete"
```

The workspace default is configured in `.domain_data/default-library.json`.
When it selects `pfam-complete`, the shorter command uses the same library:

```chimerax
domains #1 source foldseek
```

Each reference carries its validated Pfam accession, so the browser can use its
ordinary background InterPro/Pfam description, exemplar-structure and publication
lookups. Descriptions depend on what InterPro provides for that entry; the local
family description remains available as a fallback. Reused ECOD coordinates retain
their original links and provenance, while the annotation source is **Foldseek /
Pfam**. Structural similarity remains candidate evidence for family membership.

`manifest.json` selects the indexed database and reference metadata;
`coverage.json` and `build-families.tsv` record coverage. The previous ECOD/CDD
library can still be selected explicitly with `library "/path/to/.domain_data/ecod-cdd"`.

This snapshot has one template per family. Smoke tests recover Ras, EF-hand and
SPRY matches, but structural-only searches miss the RCK-like placements in
A0AAE4MDN0 and the EF-hand/channel regions in the large RyR test. Pfam guidance
recovers the latter through sequence assignments. Family coverage therefore does
not establish structural sensitivity: more reference diversity and template-boundary
checks are needed before replacing the ECOD/CDD default.
The primary workspace library instead retains the full ECOD/CDD collection with
[audited Pfam annotations](DOMAIN_LIBRARY.md#add-pfam-annotation-without-rebuilding-the-structures),
reusing the builder's saved profile hits where their provenance and natural spans
can be verified.

## Use the profiles in ChimeraX

Keep `domain_pfam.py` beside `chimerax_trimmings.py`, reload the main helper, and
use the ordinary command:

```chimerax
domains #1 source pfam
```

The default `pfam_backend auto` detects the builder's complete profile bundle
and local `hmmscan`. If local resources are unavailable or a local search fails,
it uses EMBL-EBI HMMER online instead. `pfam_backend local` requires local
resources without submitting sequences online; `pfam_backend online` explicitly
selects the service. Successful local zero-hit searches do not trigger fallback.

For a separately installed database:

```chimerax
domains #1 source pfam pfam_database "/path/to/Pfam-A.hmm"
```

Alternatively set `CHIMERAX_PFAM_DB` before starting ChimeraX. The matching
`Pfam-A.hmm.dat.gz` metadata must be beside `Pfam-A.hmm`, along with all four
current `hmmpress` indexes. A colocated `Pfam.version.gz` identifies the release;
the builder's `Pfam39.0` directory also supplies that label. Published Trimmings
libraries can supply their bundle through the manifest's `profile_database`.
For manually downloaded profiles, prepare them outside ChimeraX:

```bash
hmmpress /path/to/Pfam-A.hmm
hmmfetch --index /path/to/Pfam-A.hmm
```

Local annotation needs no Gemmi, Foldseek, or Python HMMER bindings. HMMER
executables are discovered on PATH and in common Homebrew directories. The
additional `hmmfetch` and `hmmpress` tools enable fast padded retries; full-chain
searching works with `hmmscan` alone. Profile installation is optional and no
downloads or index creation occur in the installed database during a lookup.

Local searches preserve curated sequence/domain gathering scores. A separate
reporting pass retains candidates with independent domain E-value at most 0.01
as tentative. Same-clan competing assignments are flagged, allowing explicit
metadata-defined nesting and separate repeats. Padded retries search unresolved
regions and tentative placements with 12 residues of context against families
already suggested by that sequence. Temporary subset profiles and their indexes
are removed afterward. `--max` is not used.

The Log identifies the selected backend; local details include release, coverage,
and retry coordinates. Cache identities include profile/index metadata and
executable versions by file fingerprint. Auto mode reserves online fallback time
by limiting local work to half the lookup timeout, at most 60 seconds. Explicit
local mode uses the full timeout. Cancellation stops local children and does not
start an online fallback. A fallback refetches the affected sequences together,
so one result set cannot mix local and online releases.

See [the sensitivity experiment](PFAM_RESCUE_BENCHMARK.md) for the measured
benefits and limitations of local profiles and padded retries.

## Structural coverage limitations

All-family structural coverage is a condition to verify, not a promise that
every Pfam entry has an existing or predictable fold. Pfam 39.0 includes 118
profiles shorter than 20 residues, plus motifs and disordered families. Some
may not yield a meaningful reference under the structural search criteria.
The script does not fabricate coordinates, pad a short motif with unrelated
residues, transfer labels from an unrelated fold, or count an HMM as a structure.
Complete sequence-profile coverage is available independently of this limit.

`--families PF00085 PF00622` runs an explicit pilot. It still cannot claim
full-release structural completeness. `--publish-partial` explicitly permits
publishing a **pilot**, with structural gaps retained in its manifest and report;
it is not the default. `--no-interpro` skips the non-seed-member fallback without
waiving the completeness check.

For a quick audit without downloading HMMs, seeds, or coordinates:

```bash
python3 build_pfam_library.py --audit-only --release 39.0 \
  --reuse .domain_data/ecod-cdd --out .domain_data/pfam-audit
```

This audits existing cross-references, not newly validated structural coverage.

## Use the published library

Reload `chimerax_trimmings.py`, then:

```chimerax
domains #1 source foldseek library /absolute/path/to/.domain_data/pfam-complete
```

After `--integrate-into .domain_data/ecod-cdd --make-default`, use:

```chimerax
domains #1 source foldseek
```

The default is stored in `.domain_data/default-library.json` as a relative library
path. An explicit `library` argument overrides it. Without this configuration,
Trimmings uses the original `cdd-pilot` directory.

The builder is not a runtime dependency. A distribution archive can be unpacked
elsewhere and its library directory passed to `library`. To use the complete
profile database independently with HMMER:

```bash
hmmscan --cut_ga --domtblout matches.tsv \
  .domain_data/pfam-complete/profiles/Pfam39.0/Pfam-A.hmm queries.fasta
```

The full profile bundle does not automatically add a local HMMER fallback to
ChimeraX's `source foldseek`; those remain different evidence sources.
