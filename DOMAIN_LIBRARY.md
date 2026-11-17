# Complete bulk domain library

Use the **published ECOD F40 Foldseek database** as the general-purpose structural
reference library. This avoids collecting individual structures through APIs and
uses upstream structural domain boundaries rather than cropping CDD alignment
cores. The [ECOD publication](https://doi.org/10.1093/nar/gkag1132) describes the
representative sets, bulk search databases, and Pfam mappings. The files are
available from [ECOD distributions](http://prodata.swmed.edu/ecod/distribution).

Use the [ready-to-search Trimmings library and quickstart](README.md#local-domain-library)
for our ECOD/CDD bundle. The commands below install the upstream ECOD library
without the CDD supplement. Install [Foldseek](https://github.com/steineggerlab/foldseek#installation)
first. Archives include [data attribution](DATA_NOTICE.md) and [third-party notices](THIRD_PARTY.md).

The installed v295.2 release contains **456,394 searchable targets**, all joined
to ECOD metadata. It includes experimental and predicted structures, clustered
at 40% sequence identity. Completeness here means **every target in the published
F40 search database**; it does not mean every protein or every CDD model has a
structure. Upstream lists 808 representatives with missing coordinates and 237
with empty coordinates, and excludes other unsuitable entries. These lists are
retained with the download.

## Install

Keep `domain_foldseek.py`, `domain_foldseek_catalog.py`, and
`domain_foldseek_ecod.py` together for installation. Python's standard library
and `curl` are sufficient; no Gemmi or structure prediction is needed to install
the prebuilt database. Foldseek itself is required for searching and merging.

```bash
python3 domain_foldseek.py install-ecod --out .domain_data/ecod
```

The command downloads the pinned release's database archive (about 605 MB),
domain annotations (84 MB), family/Pfam mapping, and checksum manifest. It
verifies the archive and database files, indexes all targets, and publishes an
atomic library manifest. Downloads resume after interruption; a failed import
preserves the previous publication. Unpublished failed generations are removed.

The installer also downloads CDD's small model-description table for a model-by-
model cross-reference audit. An existing catalog can be reused instead:

```bash
python3 domain_foldseek.py install-ecod --out .domain_data/ecod \
  --release v295.2 --cdd-catalog .domain_data/cdd-catalog/catalog.sqlite
```

Active search files and indexed annotations occupy about **1.3 GB**. Cached
downloads add about 0.7 GB. Previous published generations remain available to
running searches and consume additional space. Coordinates are held in the
prebuilt Foldseek database; individual full-atom models are not downloaded.

## Use in ChimeraX

Reload `chimerax_trimmings.py`, then specify the installed library:

```chimerax
domains #1 source foldseek library /absolute/path/to/.domain_data/ecod
```

Only `domain_foldseek.py` needs to accompany the main script during searches;
the build helpers are not runtime dependencies. The Protein Domains browser
shows ECOD family/group names, structural scores, and links to reference domains.
Pfam cross-references are retained in indexed metadata. These are structural
match candidates, with family and topology assignments kept distinct. A
topology-only ECOD entry does not acquire a specific family name by inference.

Exact TM-scores are used. An upstream alignment-score normalization error can
still make near-perfect matches exceed 1; the helper corrects only cases inside
the identified off-by-one bound and retains the reported score in the details.
Other invalid scores remain errors. Query sequences and coordinates remain local. Session caching, cancellation,
model switching, and Undo/Redo work as for the original pilot.

During Foldseek searches, the status bar shows elapsed seconds, model count,
and the current search pass/query-region count, updating about once per second.
A separate completion message in the status bar and Log gives elapsed time and
match/model counts. Cancellation, timeout, and failure replace the running
message. Cached results get a clearly marked completion message.

All uncached models/chains and query regions are batched into one Foldseek
invocation per pass, using four CPU threads. Repeat passes run sequentially,
because each uses the previous matches to choose remaining regions; the regions
within each pass are batched together again. Pfam searches distinct sequences
sequentially; InterPro retrieves precomputed matches in batches of up to 100 hashes.

Domain descriptions are prepared in the background and reused when navigating.
Where an ECOD family explicitly maps to Pfam, the browser prefers InterPro/Pfam
function descriptions and example PDBs; CDD/Pfam structural hits use their member
database accession. Local descriptions remain the fallback when richer details
are unavailable, and the actual matched reference remains linked. Multiple
family cross-references are retained; names alone are never used to infer a mapping.
This enriches reference-family information without changing structural match
assignments or boundaries. Pfam/InterPro details prefetch across analyzed models,
with the selected entry taking priority. Successful remote metadata is cached
for 14 days in
`.domain_data/domain-descriptions` beside the script; **Reload details** bypasses
that cache. If this directory is not writable, session caching still works.

## Add or modify references

Keep the downloaded base library intact and merge it with a custom library:

```bash
python3 domain_foldseek.py merge \
  --libraries .domain_data/ecod .domain_data/cdd-full \
  --out .domain_data/ecod-cdd
```

The existing CDD builder produces compatible custom libraries; see
[CDD build commands](DOMAIN_FOLDSEEK.md). The merger concatenates amino-acid,
3Di, C-alpha, and header databases in the same order, combines indexed metadata,
and publishes a separate manifest. It retains original catalog provenance for
each reference and rejects duplicate target identifiers. The installed ECOD plus
540-reference CDD supplement contains **456,934 targets**. Additional merges can
be chained. To update ECOD, install the new release rather than merging two ECOD
releases with overlapping IDs.

```chimerax
domains #1 source foldseek library /absolute/path/to/.domain_data/ecod-cdd
```

## CDD completeness is audited separately

The installer writes `coverage.json` and a `cdd_coverage.tsv` inside the published
database generation. For CDD v3.22, **16,791 of its 71,136 searchable models have
direct Pfam-family cross-references**; the other 54,345 have no direct mapping in
this join. That is not a count of missing folds. CDD contains overlapping
subfamilies, superfamilies, motifs, and whole-protein profiles with a different
classification from ECOD. A Pfam cross-reference is also not an experimental
demonstration of the query's CDD membership.

For exact coverage of all CDD profiles, the complete NCBI CDD sequence-profile
database is a separate annotation layer. Adding CDD labels to the structural
library should use bulk, validated sequence-profile matches and retain an
explicit per-model coverage report. A full ECOD download alone must not be called
a complete CDD-labeled library.

The CDD audit is also **not a current Pfam coverage audit**. A direct comparison
with the official **Pfam 39.0** catalog (53,275 families) finds explicit ECOD
cross-references for **19,987 families (37.5%)**, leaving **33,288** without a
direct family mapping. The CDD supplement adds no extra current Pfam families.
The previous 16,791/24,073 comparison used only the older Pfam entries in CDD.
Missing mappings do not establish that the corresponding folds are absent.
The local audit, including one row per family and a missing-family list, is in
`.domain_data/pfam-audit/v39.0/`. The official source is
[Pfam's release catalog](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/).

To target complete Pfam representation, use a pinned official Pfam catalog as
the inventory. First annotate existing reference sequences with its profiles;
then add sequence-verified PDB/AlphaFold domain crops for families still missing.
The [one-command Pfam builder](PFAM_LIBRARY.md) implements this pipeline,
including parallel retrieval, resumable profile scans, optional local prediction,
and a strict all-family structural-completeness gate before publication.
The full Pfam HMM library provides every sequence profile; it does not supply
coordinates for every family. Disordered, motif, and coiled-coil entries also
need to remain explicit in a structural-coverage report.

Foldseek also supports [ProstT5 databases built from protein FASTA sequences](https://github.com/steineggerlab/foldseek#create-custom-database-from-protein-sequence-fasta).
That can bypass protein-ID mapping and coordinate retrieval for CDD seed
representatives, but yields **predicted 3Di sequences**, not C-alpha structures.
TM-score and LDDT validation are unavailable for that representation; it is not
a drop-in replacement for the coordinate-backed library used by this browser.

## Add Pfam annotation without rebuilding the structures

`map_ecod_pfam.py` enriches an installed schema-2 ECOD or ECOD/CDD library using
a pinned Pfam catalog. It retains every target, original domain boundary, coordinate
and Foldseek structural-alphabet entry. Only reference metadata is published in a
new generation; the previous metadata and manifest remain available for inspection.

Using metadata already downloaded by the Pfam builder:

```bash
python3 map_ecod_pfam.py \
  --library .domain_data/ecod-cdd \
  --metadata .domain_data/pfam-audit/v39.0/Pfam-A.hmm.dat.gz \
  --models 53275
```

Optionally add `--checkpoint /path/to/work/generation/state.sqlite` to reuse
sequence-profile hits from an existing Pfam build. The checkpoint must refer to
this exact source-library manifest and the same Pfam metadata checksum. Imports
require curated gathering thresholds and at least 80% profile coverage and 20 matched reference residues. Hits that
cross artificial joins between discontinuous ECOD pieces, or whose natural spans
cannot be verified, are excluded. Missing experimental coordinates are allowed
only where a signature is guaranteed to remain inside one natural piece for
every possible placement of the missing residues. Upstream family mappings
remain available as cross-references for discontinuous domains. This command needs Python and the accompanying builder
modules; it does not download structures, require Gemmi or rerun Foldseek indexing.

To verify selected reference IDs against the entire installed Pfam release,
add `--scan-ids ids.txt --profiles /path/to/Pfam-A.hmm --scan-work /path/to/audit`.
This optional preparation step requires Foldseek and HMMER `hmmsearch`, exports
the actual selected database sequences, and keeps profile matches, sequence
hashes, and provenance. It publishes metadata atomically without rebuilding the
structural database. Ordinary structural searches acquire no new dependencies.

The browser promotes a Pfam name/accession only when a gathering-threshold
signature covers at least 80% of its profile, accounts for a majority of the
paired reference residues, and has at least 50% of its signature positions
represented in the alignment. This permits genuine structural insertions around
a complete signature without transferring small partial motifs to a whole fold.
Repeated copies of one signature are considered together. Overlapping signatures
within one curated Pfam clan can resolve to the strongest sequence evidence;
unrelated competing families remain ambiguous. A family-level ECOD cross-reference alone is insufficient:
it can describe the source protein rather than the matched structural domain.
Localized sequence evidence can override an inconsistent upstream cross-reference,
with the disagreement disclosed in the details. Unsupported or ambiguous mappings
use the ECOD topology or homology-group label. Their candidate Pfam links remain explicitly unverified and
do not supply the region's primary function description. Details retain the
mapping evidence and original ECOD link.
Cached labels, accessions and descriptions are recomputed together from the actual
matched reference. Unverified links never select a functional description, and
InterPro responses must match the requested member/integrated accession.
The query is still a structural-match candidate, not a sequence-confirmed member.

The mapping audit and per-target TSV are linked from `manifest.json` under
`pfam_mapping`. An unresolved entry has no accepted mapping in the available
evidence; it has not necessarily been exhaustively searched against Pfam. Cached
scan reuse is recorded separately from a complete sequence rescan. No nearest-family
labels are inferred from fold similarity alone.

The initial workspace mapping pass added evidence for **3,090 previously unmapped ECOD
references**, representing **22,550 current Pfam families** across the full
collection. Subsequent targeted reference scans are recorded separately in the
manifest audit; they are not a complete rescan of every template.
The structural database remains **456,934 targets**, including the CDD
supplement. RCK-like placements in A0AAE4MDN0 retain their structural matches;
individual Pfam labels require the localized support above. The smaller Pfam-only library remains an
explicit alternative rather than replacing ECOD's reference diversity.

## Verification

In the full merged library, thioredoxin is recovered at E-value approximately
**4 × 10⁻¹³**, without relaxing the 0.01 cutoff. Native ChimeraX checks also recover
Src SH3, SH2, and kinase domains, and the lysozyme fold that was absent from the
small pilot. These three examples took roughly **1.6–9 seconds** on this machine,
including repeat searches and exact TM-scores; Undo/Redo was verified.

Residual searches reserve only the best nonredundant placements, so alternative
family hits cannot consume adjacent unassigned regions and prevent a neighboring
domain from being found. In the browser, different-family candidates sharing at
least 70% of both aligned footprints appear under one primary row; alternatives
remain listed with scores and links in its details. Small boundary overlaps and
distinct repeats remain separate. Nearby/Selected tabs identify annotations at
the center of rotation or overlapping a residue selection without moving the view.
This remains a structural search, without calibrated function-specific confidence.
