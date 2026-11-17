# Local structural domain matching

For consistent Pfam family names and InterPro descriptions, use the
[published Pfam structural library](PFAM_LIBRARY.md#use-the-published-pfam-foldseek-library).
It includes verified experimental and predicted templates and an explicit family
coverage audit. The database is usable while unresolved families remain gaps.

For a ready-made, complete published structural search database, use the
[bulk ECOD installer and merge workflow](DOMAIN_LIBRARY.md). The CDD-only
builders below remain useful for supplements and coverage experiments.

The Foldseek pilot searches model coordinates against a small library of
experimental domain structures labeled by their CDD family membership.
Results appear in the existing Protein Domains browser as **structural match
candidates**, with alignment coverage, TM-scores, family descriptions, and links
to the matched PDB entries. The query's biochemical function is not established
by a structural match alone.

The current pilot contains **37 domains covering 14 CDD models**. It includes
both NCBI-curated models and Pfam models imported by CDD, with their provenance
preserved. It is a test of the workflow, not a comprehensive CDD database. It
does not yet include AlphaFold reference structures or automatic segmentation.

A separate **CDD-wide builder** now indexes the official NCBI seed alignments
across every source in CDD and can use both PDB structures and AlphaFold
predictions. Its library is built separately in `.domain_data/cdd-full`; the
original pilot remains available for comparisons. All retained reference files
contain only C-alpha atoms, with residue names and sequence order preserved.

## Use in ChimeraX

Keep `domain_foldseek.py` beside `chimerax_trimmings.py`, install the `foldseek`
executable, and reload the main script. No additional Startup command or Python
package inside ChimeraX is needed. The pilot built in this workspace is ready
to use:

```chimerax
domains #1 source foldseek
domains #1:702-855 source foldseek
domains sel source foldseek color false
```

For this source, the requested residues define the structural query. Whole-chain
searches can find local domain matches. To search several proposed domains
separately, supply sequence positions using commas between domains and underscores
between discontinuous pieces of one domain:

```chimerax
domains #1 source foldseek regions "13-345,388-414_604-690,422-558_585-592,702-855"
```

These are full chain **sequence positions**, not author residue numbers. Each
query needs at least 20 modeled C-alpha atoms. An explicit position map retains
missing residues, numbering offsets, insertion codes, and discontinuous pieces.
Foldseek annotations color and select the whole matched sequence range, including
internal alignment gaps. Explicit discontinuous query pieces and the requested
atom selection remain separate limits. Missing residues have no atoms to color.
The paired alignment positions are retained separately for scoring, retry
decisions, and overlap grouping; filling the displayed range does not alter the
search or scores. Details report the number of structurally paired residues.

Where sufficient modeled core support exists, neighboring crossing structural
assignments share one display boundary, so coloring and selection agree.
A local C-alpha contact graph can
move that boundary toward a cut with fewer contacts between domains; cuts through
an assigned helix or strand are avoided when a permitted loop cut exists.
Changes stay near the original boundary, retain domain cores, and respect the
explicit query pieces. Weak contact evidence favors the existing score-based
overlap priority. Inserted domains, discontinuous matches, and sequence-derived
annotations keep their own boundaries. Details identify adjusted display spans;
the original alignments and scores remain unchanged. This is a local display
heuristic, not a guarantee of biologically exact domain boundaries.

The default library is selected by `.domain_data/default-library.json`, with
`.domain_data/cdd-pilot` as the fallback. Override it with
`library /absolute/path/to/library`. Searches run locally in the
background; `wait true`, `refresh true`, `domains cancel`, `gui false`, and
Undo/Redo follow the other domain commands. Cache keys include query coordinates,
sequence, queried regions, the library manifest, and the search algorithm version.
Editing coordinates, rebuilding the library, or updating the search algorithm
causes a fresh search. Model placement alone does not.

By default, chains with identical sequences and the same modeled query positions
share one representative search, including copies in different models. Results
are mapped to each copy's own residues and author numbering. Transferred results
identify the representative in their details. Different missing-residue patterns
or requested regions remain separate queries. Use `deduplicate false` to search
each chain's conformation independently; sequence identity alone does not guarantee
identical structural scores. The Log reports chain and representative counts.

To try sequence-guided structural searching:

```chimerax
domains #1 source foldseek pfam_guided true
```

This uses the same automatic local-Pfam/online fallback as `source pfam`, reusing
the session's compatible Pfam results where available. `pfam_backend online`
selects EMBL-EBI explicitly; `pfam_backend local` requires local profiles.
`pfam_database` selects a different local installation. Only curated,
non-outcompeted matches guide boundaries. The first pass searches those regions with 12 modeled residues
of context at each end, plus all remaining regions; no sequence-negative region
is excluded. Initial screening uses E-value 0.01 and the normal structural
acceptance filters. Unresolved stretches retain the existing loose candidate
search, focused retries, and broader final review, with up to five passes total.
Chains without useful Pfam boundaries retain the ordinary search strategy when
searched alone. Discontinuous domains may still require explicit `regions`.

Curated Pfam assignments take precedence over substantially overlapping
structural rows; their details retain the Foldseek assignments, scores, and links.
Separate repeats and substantially smaller inserted domains remain separate.
Sequence and structural E-values are not compared across methods. The two search
modes and Pfam backend/database identities have separate cache keys. `refresh true`
refreshes both sequence guidance and structural results. If Pfam lookup fails,
the full local Foldseek search runs instead, with a warning; that fallback is
not cached as successful guidance.
Guidance lookup has a maximum of 60 seconds or one-third of the overall timeout,
whichever is smaller. Auto mode reserves part of that budget for online fallback.
Ordinary unguided Foldseek searching remains fully local.

On the large twelve-chain RyR complex (`ATP4CmC_fromfatiguer10_real_space_refined_011.pdb`),
using the same 456,934-reference library and four CPU threads:

| Strategy | Structural search | Pfam lookup |
| --- | ---: | ---: |
| Previous: all twelve chains independently | 284 s | — |
| Default: three sequence representatives | 131 s | — |
| Pfam-guided: three representatives | 66 s | 2.9 s |

These are individual sequential measurements, excluding interface rendering.
Representative-only results matched the previous run's results for those three
chains exactly. They do not represent separate measurements of every conformation.
Guidance recovered all four Pfam RyR repeats per subunit and retained the structural
EF-hand match at author residues 4070–4133. It changed other boundaries and coverage:
the previous broad terminal channel match reached residue 5032, while the guided
browser's curated Ion_trans assignment ends at 4947. Guidance therefore remains
opt-in; this benchmark does not establish equivalent sensitivity or complete
domain coverage. Native mapping, per-chain transfer, cached reruns, and Undo/Redo
were also checked on the complete complex using these real search results.

Pilot candidate filters require E-value <= 0.01, reference coverage >= 50%,
alignment-normalized TM-score >= 0.5, and at least 20 paired query residues. The
best qualifying representative for each distinct family occurrence is retained.
Same-family hits are treated as duplicates when their paired residue footprints
overlap by at least 50% of the smaller footprint. Small boundary overlaps and
separate occurrences of the same family are retained as separate browser rows.
The browser also groups different-family assignments when at least 70% of both
aligned footprints coincide. To accommodate different alignment gaps, it also
groups hits sharing at least 50% of both footprints when their outer spans
overlap by at least 80% of both lengths. Span overlap alone is insufficient.
A weaker partial match also groups when at least 90% of its paired residues
and 90% of its outer span are contained in a stronger primary match, provided
it has at least half as many paired residues or shares an explicit Pfam family
or ECOD homology group with that primary. This handles contained partial-fold
matches and ECOD/Pfam differences in domain granularity. Family relationships
come from reference cross-references and ECOD classification, rather than names.
Related names alone do not establish a match, and a domain aligned within an
unpaired insertion remains separate. Grouping compares each hit directly with
the primary; alternatives cannot expand a group or connect adjacent domains.
The lowest E-value (then highest bit score) defines
the primary footprint and coloring; other family candidates remain linked in its
details. Separate repeats and substantially smaller unrelated inserted domains
remain separate. The search's raw candidates stay cached for reanalysis.
After updating `chimerax_trimmings.py`, reload it with `runscript` before rerunning
`domains` to use updated browser grouping. `refresh true` repeats the structural
search but does not reload the main helper. Cached hits are regrouped when the
command is rerun, so a fresh structural search is unnecessary for grouping updates.
These are search filters, not CDD's calibrated sequence-family confidence thresholds.
E-values depend on the searched library; a small-pilot result cannot be assumed
significant against the full CDD collection.

Foldseek can return only one alignment to a given reference in a search. To find
repeats, the helper searches remaining stretches of at least 20 modeled residues
again, within the overall timeout. This uses the same score
filters on every pass. It also retries up to four distinct plausible placements
per chain per pass when their E-values fail the final cutoff. These must already
have alignment TM-score at least 0.5, reference coverage at least 50%, at least 20
paired residues, and E-value at most 10. The retry spans the aligned query
residues with 12 modeled residues of context on each side, bounded by the original
query region. Placements overlapping accepted annotations by at least half their
paired residues are skipped, and retries with at least 70% mutual window overlap
are combined. Every displayed match still has to pass the original E-value
cutoff of 0.01 and all geometry filters. This prevents an empty first pass from
ending the search before domain-sized regions have been checked.
Normal passes allow 1,000 references per query through Foldseek's prefilter,
matching Foldseek's default candidate budget. A broader final review uses 5,000
for the remaining regions, including regions previously tested at the lower
budget. It also rechecks up to four provisional assignments per original query,
prioritizing the least certain. Their query windows retain local context but
stop at stronger neighboring matches. An early weaker assignment therefore
cannot prevent a better family match from being found in the same region.
The review starts when normal retries finish or as the fourth pass (fifth with
Pfam guidance). After review, newly exposed residual regions continue at the
5,000-candidate depth until no untested regions remain. There is no fixed
four-copy/four-pass limit; tested-region tracking prevents unchanged searches
from looping, and the shared timeout and cancellation still bound the work.
It reports only hits meeting the final E-value cutoff, avoiding
expensive detailed scoring of thousands of rejected alignments. Every displayed
match still requires all geometry filters. This increases search sensitivity
without lowering annotation thresholds. Searches of AFDB50 on the Foldseek webserver use a different reference
collection and may find homologs absent from this ECOD-based library.
Query residues are not filtered by pLDDT/B-factor; changing query confidence
alone does not change the structural cache key. A poorly predicted fold can
nevertheless fail the geometric match filters. The template builders' AlphaFold
quality cutoffs apply to reference-library construction, not the loaded query.
In a controlled real-Foldseek regression check using an artificial chain of
eight copies of a 113-residue globin reference, the old four-pass behavior found
four copies; continued residual searches found all eight in 5.3 seconds.
This used a single-reference test library, not the full ECOD/Pfam database.
Discontinuous domains can still need explicit `regions`;
this is not a general domain-boundary predictor. For a repeat-search hit, the
browser identifies the searched segment: query coverage and query-normalized
TM-score refer to that segment. Residue annotations always use the original
chain positions.

The search calculates exact alignment-normalized TM-scores for every reported
candidate. Query- and reference-normalized TM-scores are display-only and are
calculated with the same exact method for accepted candidates. They do not
control retries, family selection, or overlap handling. Rejected candidates keep
uncalculated display fields unset rather than receiving approximate or invented scores.
This avoids expensive unused calculations while preserving the search depth,
alignments, acceptance thresholds, and full scores of accepted matches. Custom
libraries without a usable target lookup retain the full-score conversion path.

Large exact-score conversions are divided into independent target batches to
balance work between CPU cores. The default remains four CPU threads; workers
share that budget. Candidate discovery retains its established
thread budget. Scoring workers share one overall CPU limit and deadline; their
rows are restored to the original alignment order before retries or match
selection, including equal-score ties. Small results use a single conversion.
This changes scheduling rather than the candidate pool, alignment method, or
score calculations. Cancellation or a failed worker stops the remaining workers.

For one or two long queries against a large, unindexed library, alignment is also
divided across candidate lists. Every worker uses the complete target database,
preserving E-value statistics. The merged alignments restore Foldseek's global
score, target-length, and target-key ordering before exact scoring and adaptive
searches. Many-query batches and indexed custom libraries retain the standard
Foldseek search path. Both paths share the same four-thread default budget.
On the large RyR chain, a back-to-back check took 129.28 seconds with ordinary
alignment and 80.08 seconds with candidate-parallel alignment (38% less time).
All 3,474 candidates, their order, scores, and five search passes were identical.
Timing depends on the query, library, hardware, and filesystem cache.

An opt-in native ChimeraX benchmark, `tests/benchmark_domain_boundaries.py`,
compares display boundaries with downloaded PDB/CATH assignments. Supply local
`ID.pdb` and `ID-cath.json` files with `--inputs`, the structural library with
`--library`, and a results directory with `--out`. It records input hashes and
source URLs, retains raw search results, and reports paired and displayed
footprint overlap separately. Unchanged queries and library metadata reuse the
saved search, allowing boundary rules to be evaluated without rerunning Foldseek.
The eight-structure panel covers Src, adenylate kinase, calmodulin, GroEL, a TIM
barrel, KcsA, a RyR fold, and four tandem fibronectin domains. This is an
architecture/boundary control panel, not a held-out remote-homology benchmark.
Composite references and discontinuous domains remain limitations; CATH and
Pfam also differ in the granularity of some domain definitions.

## Build or repeat the pilot

The builder requires a separate Python with `gemmi`, plus the Foldseek executable:

```bash
python3 domain_foldseek.py build --out .domain_data/cdd-pilot --representatives 3
```

The builder reads CDD's annotation table and retrieves member-database structure
locations through the InterPro API. It downloads experimental mmCIF structures
from RCSB, verifies chain sequence and `label_seq_id` mappings, and extracts
individual domain occurrences. Crops need at least 25 C-alpha atoms and 80%
modeled coverage of the reported location. Exact duplicate reference sequences
within a family are skipped.

To keep this first build small, it examines a limited set of structure candidates
and skips parent chains longer than 2,000 residues. Three requested models lacked
suitable indexed references in that sample; the manifest records missing
families and rejected crops. Expanding coverage will require broader member
selection and AlphaFold representatives. Reference downloads are cached for
rebuilding; library files and downloads are excluded from Git.

## Expand toward the full CDD collection

Keep `domain_foldseek_catalog.py` beside `domain_foldseek.py` for these build
commands. It is not needed for ChimeraX searches. The CDD v3.22 seed archive
contains **71,136 searchable models**. The additional 5,580 `cl` entries are
derived multi-model superfamilies, rather than independent seed alignments.

Download and index the public seed archive once (about 767 MB compressed):

```bash
mkdir -p .domain_data/cdd-catalog
curl --fail --location --continue-at - \
  --output .domain_data/cdd-catalog/fasta.tar.gz \
  https://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz
python3 domain_foldseek.py catalog --out .domain_data/cdd-catalog \
  --archive .domain_data/cdd-catalog/fasta.tar.gz
```

The catalog records the archive's SHA-256 and inventories every model source.
With up to eight candidate sequences per structure source per family, **44,692
models have usable PDB or UniProt seed identifiers**. The other 26,444 currently
have no directly usable candidate in this index. Many require identifier mapping
from RefSeq/GenBank to UniProt; seed fragments outside the supported 25–9,999
residue range are also excluded. This is not evidence that those families lack structures.
NCBI consensus sequences are never treated as natural proteins or AlphaFold
references.

Build a bounded batch, retaining the original pilot references:

```bash
python3 domain_foldseek.py expand --out .domain_data/cdd-full \
  --catalog .domain_data/cdd-catalog/catalog.sqlite \
  --seed-library .domain_data/cdd-pilot --limit 200
```

Repeat the same command to resume with new families. Omit `--limit` to process
all currently resolvable families. Batches interleave source collections;
imported families are completed first and already checkpointed families are
skipped. Ctrl-C checkpoints progress and
publishes completed references. Families with unavailable candidates or transient
errors retain explicit statuses in `build.sqlite`; `--retry-failed` retries
transient failures. A sampled set of candidates can fail even when another
usable reference exists.

The builder verifies an exact, unique placement of each CDD seed fragment in
the structure's polymer sequence, then crops using `label_seq_id`. Identical
copies in different assembly chains are interchangeable, but repeated placements
within one chain are rejected as ambiguous. Crops require at least 25 modeled
C-alpha atoms and 80% modeled coverage. AlphaFold crops additionally require
mean domain pLDDT >= 70; `--min-plddt` changes this build-time filter, and
`--no-alphafold` restricts a batch to experimental structures. Prediction links,
confidence, and seed provenance remain distinct from experimental evidence.

`--representatives` sets the target references per family (default 3); exact
duplicate crop sequences within a family are skipped. `--models pfam00042
smart00020` targets specified models, for example to deepen representation of a
particular fold. Existing imported references are retained. Candidate sampling
is controlled by `catalog --candidates`; increasing it requires rebuilding the
catalog. Changes to catalog data or build settings invalidate affected build
checkpoints.

Each publication creates an immutable Foldseek database and indexed SQLite
metadata, then replaces the small manifest atomically. An interrupted build
leaves the previous publication searchable. Searches load metadata only for
returned hits and keep one database generation across repeat passes. Old
generations are retained so running searches can finish. Full coordinate
downloads remain in the build cache; C-alpha crops and Foldseek database files
are the retained search representations.

Inspect coverage and search the expanded library explicitly:

```bash
python3 domain_foldseek.py status --library .domain_data/cdd-full
```

```chimerax
domains #1 source foldseek library /absolute/path/to/.domain_data/cdd-full
```

The expanded library is still a structural-candidate search, without calibrated
CDD sequence-family confidence. E-values change as the reference collection
grows; the existing search filters are unchanged. Automatic domain segmentation
and complete mapping of RefSeq/GenBank seeds remain separate work.

## Expanded library measurement

The first general expansion contains **540 references across 202 models**:
379 experimental PDB crops and 161 AlphaFold crops. It covers CDD-curated,
Pfam, SMART, COG, TIGRFAM, NCBIFAM, and `sd` profiles. The published Foldseek
files occupy about 0.8 MB and indexed reference metadata 1.0 MB; the C-alpha
crop directory occupies 6.9 MB. Full download caches, the seed catalog/archive,
and previous database generations are additional build storage.

Run the reproducible general-purpose panel with:

```bash
python3 domain_foldseek.py benchmark --library .domain_data/cdd-full \
  --out .domain_data/cdd-full/benchmark/general
```

Queries are public AlphaFold models. References with exactly matching crop
sequences in each query are excluded. This does not remove all close homologs
and is not a remote-homology or calibrated accuracy benchmark.

| Control | Observed structural candidates |
| --- | --- |
| Hemoglobin beta | Globin |
| H-Ras | Ras/small GTPase, plus related Arf and EF-Tu GTPase folds |
| Thioredoxin | TRX-superfamily hit filtered: E-value 0.0638 despite alignment TM-score 0.784 |
| Calmodulin | Two EF-hand occurrences |
| Ferritin heavy chain | Ferritin-like |
| Trypsin | Trypsin-like serine protease |
| Src | SH3 and SH2; kinase coverage remains absent in this batch |
| Lysozyme, an unrepresented fold control | No qualifying hit |

Searches took 0.56–1.29 seconds, median 0.87 seconds, excluding downloads and
constructing the per-query control database. Native ChimeraX checks for Ras,
calmodulin, and trypsin took 0.57–1.13 seconds and verified caching and Undo/Redo.
Names describe the reference families: the related GTPase alternatives illustrate
why a fold match is not sufficient to establish a specific family or function.

For thioredoxin, the best raw match is a 74-residue PDI-family reference crop:
95.9% reference coverage, 64 paired residues, and alignment TM-score 0.784.
Only its E-value fails the 0.01 candidate cutoff. The three current TRX-superfamily
references are short 61–74-residue crops, versus the 109-residue query. CDD seed
fragments can cover the conserved core rather than a complete structural domain;
reference boundaries and representative diversity therefore need validation too.

The default remains the original pilot while broader validation continues.
Increasing the library size changes E-values: D1Z184's best whole-chain RyR hit,
for example, changes from 0.00522 to approximately 0.014 and no longer passes
the fixed 0.01 cutoff. This benchmark has not been tuned to recover that case.
Broader representative diversity, identifier mapping, sensitivity, and competing
family assignments need validation before treating this as a full annotation tool.

Command-line searches retain all raw hits, including weak ones rejected by the
ChimeraX candidate filters:

```bash
python3 domain_foldseek.py search model.cif --chain A \
  --library .domain_data/cdd-pilot --segments '702-855' --out hits.json
```

## Original pilot measurement

With repeat searching enabled and four CPU threads:

| Search | Median total time, three runs |
| --- | ---: |
| A0A662RAG4, four supplied regions | 1.14 seconds |
| A0A662RAG4, whole chain | 1.13 seconds |
| D1Z184, whole chain | 1.72 seconds |

These timings include reading the query mmCIF, exporting queries, starting
Foldseek, searching, and parsing the hits. They do not include retrieving or
predicting domain boundaries. On D1Z184, the native ChimeraX command took
approximately 2.0 seconds initially and 0.02 seconds with cached matches;
Undo/Redo was verified for both repeat annotations.
Those are historical pilot timings. A current check against the 456,934-reference
ECOD/CDD library used the 643-residue A0AAE4MDN0 model, four CPU threads, and all four
search passes. Deferring display-only scores reduced the complete run from 60.98
to 32.04 seconds. All 2,342 candidate alignments and selection scores were
identical, as were the full scores of all 24 accepted candidate hits. Both
TrkA_N occurrences and the RyR fold were retained. The slowest individual pass
fell from 35.44 to 17.75 seconds. These are individual runs on the development
machine, not a general runtime guarantee.

A subsequent back-to-back check measured 32.15 seconds with serial score
conversion and 24.95 seconds with balanced conversion, using the same four-thread
budget. All 2,342 candidates, their order, and their scores were identical. Eight
threads reduced that run to 24.02 seconds, so four remains the default. The
optimizations preserve both TrkA_N occurrences and the RyR fold.

For D1Z184, the original whole-chain pass finds a RyR candidate at 497–580
(E-value 0.00522). Searching the remaining region recovers a second RyR candidate
at 428–488 (E-value 0.00287). Both are retained, without changing score thresholds.
These spans summarize the matched regions. Internal alignment gaps are filled
for coloring and selection, while structural scores retain the paired residues.

For A0A662RAG4 with the supplied regions:

| Region | Best reference family | E-value | Alignment TM-score | Interpretation |
| --- | --- | ---: | ---: | --- |
| 13-345 | Ras | 1.12 | 0.270 | Weak; not annotated |
| 388-414, 604-690 | Globin | 3.27 | 0.251 | Weak; not annotated |
| 422-558, 585-592 | AAA | 0.000488 | 0.604 | Structural candidate |
| 702-855 | RyR | 0.000414 | 0.770 | Partial structural candidate; 44.8% query coverage |

The split queries improve some match scores here, but do not reveal an additional
family beyond the existing AAA and RyR sequence matches. The two other regions
remain unassigned in this small reference collection. Globin controls are
present and do not displace the RyR match.

The pilot has 3,942 reference residues. Its Foldseek search files occupy about
37 KB, cropped PDB files 346 KB, and cached public downloads 59 MB. A full
collection will need substantially more references; the pilot's speed and
specificity should not be extrapolated directly to that collection.
