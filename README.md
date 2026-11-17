# ChimeraX Trimmings

A compact UCSF ChimeraX setup for interactive model and map work.

[MIT licensed](LICENSE). Database attribution and optional-service terms: [THIRD_PARTY.md](THIRD_PARTY.md).

This repository contains:

- [chimerax_trimmings.txt](chimerax_trimmings.txt): startup commands, aliases, keybindings, and button-panel setup
- [chimerax_trimmings.py](chimerax_trimmings.py): companion ChimeraX Python commands, including fast shape alignment and experimental whole-structure gradient scoring
- [map_salami.py](map_salami.py): map/model-fit panel and movie export

## Installation

1. Download or clone this repository to a permanent folder. Keep its Python files together.
2. Edit the first two lines of [chimerax_trimmings.txt](chimerax_trimmings.txt) so the `runscript` paths point to your local copies of the Python scripts. Remove any older separate `runscript` line for `shape_align.py`; its commands are now included in the main script.
3. In ChimeraX on macOS, open `UCSF ChimeraX -> Preferences... -> Startup`.
4. Paste the contents of [chimerax_trimmings.txt](chimerax_trimmings.txt) into the startup commands box, or source it from a startup file you already maintain.
5. Restart ChimeraX.

Windows note: the Startup preferences are typically reached from the `Favorites` menu rather than the application menu.

## Local domain library

1. Install [Foldseek](https://github.com/steineggerlab/foldseek#installation)
   (macOS: `brew install foldseek`; Linux: follow the linked instructions).
2. Download the [library](https://github.com/olibclarke/chimerax-trimmings/releases/download/domains-2026-11-17/trimmings-ecod-cdd-v295.2.tar.gz)
   ([SHA256](https://github.com/olibclarke/chimerax-trimmings/releases/download/domains-2026-11-17/trimmings-ecod-cdd-v295.2.tar.gz.sha256)).
   From this repository's folder, extract it:

```bash
mkdir -p .domain_data
tar -xzf ~/Downloads/trimmings-ecod-cdd-v295.2.tar.gz -C .domain_data
```

3. In ChimeraX, after loading the helper:

```chimerax
domains #1 source foldseek library "/absolute/path/to/chimerax-trimmings/.domain_data/ecod-cdd"
```

The bundle has 456,934 ECOD/CDD references with Pfam annotation support;
it does not cover every Pfam family. Searching needs neither HMMER nor Gemmi.
Browse matches, fetch/align exemplar domains, or show their full-structure context.
Queries stay local; descriptions/exemplars use public services. Stop with `domains cancel`.
[Library details and upstream alternative](DOMAIN_LIBRARY.md) · [Search options](DOMAIN_FOLDSEEK.md) · [Build Pfam templates](PFAM_LIBRARY.md).

## Companion Script

The Python script registers:

- `nextmodel`
- `prevmodel`
- `nextitem`
- `previtem`
- `lowermapthreshold`
- `raisemapthreshold`
- `togglemaps`
- `togglemodels`
- `cyclelighting`
- `cyclemodeldisplay`
- `exposedhydrophobics`
- `hydrophobicpatches`
- `patchbrowser`
- `shapealign`
- `domains`
- `domains cancel`
- `domainbrowser`
- `tmsegments`
- `tmsegments cancel`

`togglemaps` and `togglemodels` are stateful hide/restore commands. `cyclemodeldisplay` cycles the currently shown atomic models through the saved original display state plus several preset representations.

## Fast Shape Alignment

Use `shapealign` to superimpose atomic models and density maps without sequence
alignment. Model/model, model/map, map/model, and map/map fitting use the same
syntax:

```chimerax
shapealign #2 to #1
shapealign #2 to #1 objective gradient
shapealign #2-5 to #1
shapealign #2/A to #1/B
shapealign #2 to #1 radius_search
shapealign #2 to #1 radius_search 60
shapealign #2 to #1 radius_refine 30
shapealign #2 to #1 color_fit true
```

The first input moves; `to` identifies one fixed reference. Multiple moving
inputs align independently. Atomic selections restrict fitting atoms but the
whole moving model is repositioned.
The reference is skipped automatically if it is included in the moving range.

Known AlphaFold database models are automatically weighted by pLDDT to reduce
the influence of uncertain tails. For local prediction files, use
`confidence_moving plddt`, `confidence_reference plddt`, or `confidence plddt`
for both inputs. These settings expect pLDDT in the B-factor field, not
experimental B-factors. `confidence off` disables weighting. The whole model
still moves rigidly; maps are unchanged.

Global atomic batches with at least four moving models automatically use the
available logical CPU count minus one, capped by the number of moving models.
Use `workers 1` to run serially or `workers N` to request a process count.
Map fits remain serial. Results retain one Undo action.

Unless a CoFR sphere is requested, maps use full density, independent of contours
and display regions. Map values,
grid metadata, styling, and visibility are unchanged; only the moving owner's
rigid transform changes.

The default `objective density` compares density overlap. Try the experimental
`objective gradient` when bulk overlap gives an implausible orientation.
`radius_search` globally searches a reference-map sphere around the current
center of rotation, keeping the moving selection intact. Omit its value to use
1.5 times the moving input's enclosing radius. `radius_refine` locally
refines an existing alignment using both inputs' density inside that sphere.
Both radii are in Angstroms and cannot be combined. `color_fit true` colors moving atomic
models by local density agreement; it is not available for moving maps.

Type `shapealign` or `shapealign help` to open the Log with a short, formatted
summary. Click **All options** or run `usage shapealign` for complete syntax.
The command is included in `chimerax_trimmings.py` and supports Undo/Redo.
See [Shape Alignment](SHAPE_ALIGN.md) for installation, all options, examples,
and practical limitations.

## Domain Coloring and Annotation

`domains` colors protein residues by domain/family and opens a dockable
**Protein Domains** browser. It is included in `chimerax_trimmings.py`.

```chimerax
domains #1
domains #1 source interpro
domains #1 source foldseek
domains #1-3
domains #1/A color false
domains #1 source pfam include_tentative false
domains #1 source pfam pfam_backend online
domains #1 features false
domains #1 source features
domains #1 source features prosite_profiles true
domains #1 source tmsegments
domains #1 source tmsegments results /path/predicted_topologies.3line
domainbrowser #2
domains cancel
```

**Sources.** Pfam automatically uses an installed local Pfam/HMMER bundle when
available, otherwise it sends each distinct full chain sequence to the
[EMBL-EBI HMMER API](https://hmmer-web-docs.readthedocs.io/en/latest/api.html).
Both use curated family cutoffs. Local searches also retry padded regions with
ordinary HMMER settings and retain reported tentative candidates. InterPro sends sequence
checksums to the [InterPro Matches API](https://github.com/ProteinsWebTeam/interpro-matches-api)
and retrieves existing matches from databases including SUPERFAMILY, CATH-Gene3D,
Pfam, SMART, and CDD. It can provide broader coverage, but requires the exact
sequence to be indexed; use Pfam for novel or modified sequences.

Local Pfam is optional: keep `domain_pfam.py` beside the main helper and install
HMMER plus a pressed Pfam database with matching family metadata. The
[Pfam builder](PFAM_LIBRARY.md#use-the-profiles-in-chimerax) prepares this bundle
before structural-library construction finishes. Automatic discovery checks
`CHIMERAX_PFAM_DB`, profiles bundled with the configured structural library, and
`.domain_data/*/profiles/Pfam*/Pfam-A.hmm`. Use `pfam_database "/path/Pfam-A.hmm"`
for another installation. No profiles are downloaded automatically.

`pfam_backend auto` is the default. An unavailable or failed local backend falls
back to the online service, with the backend recorded in the Log. A valid empty
local result stays local. `pfam_backend local` requires local resources and never
submits sequences online; `pfam_backend online` always uses the service. These
settings also apply to `pfam_guided true` Foldseek searches. Local details show
the profile release, profile coverage, and any padded retry range. Changing the
backend, database, or local executables invalidates the search cache. Online and
local results can differ because their Pfam releases differ.

An optional **local Foldseek search** compares the loaded coordinates with
domain references. The workspace default retains ECOD/CDD's complete structural
reference collection with [audited Pfam mappings](DOMAIN_LIBRARY.md#add-pfam-annotation-without-rebuilding-the-structures).
Matches with localized sequence-profile support use Pfam names, accessions and InterPro descriptions; original
ECOD classifications and reference links remain in the details. Localized,
complete profile signatures take precedence over inconsistent family-wide
cross-references. Unsupported mappings retain fold-level ECOD labels and do not
supply Pfam functional descriptions. Names, accessions, descriptions and example
structures use the same validated family, including when reusing cached hits.
A [Pfam-only library](PFAM_LIBRARY.md#use-the-published-pfam-foldseek-library)
uses consistent Pfam accessions for InterPro descriptions and example structures,
including templates whose coordinates originally came from ECOD. The verified
workspace snapshot covers **46,654 of 52,484** folded-scope families; its coverage
audit records the remaining gaps. A [bulk installer](DOMAIN_LIBRARY.md) provides the complete
published ECOD F40 database (456,394 targets) and supports merging custom
references. A small experimental CDD pilot is available, plus a
resumable CDD-wide builder using sequence-verified PDB and AlphaFold crops.
A [Pfam-first builder](PFAM_LIBRARY.md) bundles every profile in a pinned Pfam
release and fills structural-reference gaps with a resumable, audited pipeline.
The [InterPro-guided completion pass](PFAM_LIBRARY.md#interpro-guided-folded-domain-completion)
uses annotated PDB/AlphaFold occurrences to fill gaps in a folded-domain target,
keeping motifs and disorder separate from structural templates. Full Pfam
sequence profiles remain available regardless of structural scope.
Both use C-alpha-only references and display structural candidates, scores,
and reference links in the same browser. Foldseek retries unassigned stretches
and plausible candidate regions at domain scale, using the same acceptance
cutoffs on every pass. A broader final review revisits provisional assignments
so weaker early matches cannot block better alternatives. Exact display-only
scores are deferred until a candidate passes the annotation filters. Independent
exact-score calculations are balanced across the existing four-thread CPU budget.
Sequence-identical copies with the same modeled query positions share one
representative search by default; results map back to each copy's own residues.
Use `deduplicate false` to compare their individual conformations.
An optional `pfam_guided true` first searches distinct sequences with the selected
Pfam backend, then searches padded Pfam regions and every remaining stretch with Foldseek.
Curated Pfam assignments retain overlapping Foldseek evidence in their details.
Unresolved regions still receive the extensive structural fallback search.
See
[Local structural domain matching](DOMAIN_FOLDSEEK.md) for setup, selected-region
queries, the pilot benchmark, and its limited coverage.

Matches describe domains, repeats, families, and superfamilies; they may overlap
or leave regions unannotated. Boundaries start from database alignments, retaining
missing residues and insertion codes. Foldseek colors and selects the complete matched range, filling
internal alignment gaps while preserving explicitly separate query pieces.
Adjacent structural matches can share a display boundary, refined locally using
C-alpha contacts and secondary structure when supported. Inserted domains and
sequence annotations remain separate. Scoring retains the actual paired residues;
details identify adjusted display boundaries. There is no adjustable confidence threshold. Matching annotations
share colors across models; unassigned target residues are gray. Only residues
in the requested selection receive new annotations and attributes.
Family colors are allocated together using perceptual color and hue separation,
with extra contrast for neighbors in the sequence or in space.
Different signatures keep separate colors even when they share an integrated
InterPro entry; copies of the same signature retain the same color across models.
If a new neighboring assignment requires a color change, existing analyzed copies
are updated together, with Undo/Redo support. `color false` preserves colors.
Browser swatches reflect the actual ribbon colors. Fully covered alternative
assignments have a crossed, empty swatch and an explanation instead of an unused
family color; higher-priority assignments control their overlapping residues.

**Browser.** Results persist across searches. Choose a model from the **Model**
dropdown or use **< / >** to cycle through analyzed models. The model arrow
buttons automatically show only the selected atomic model; Undo restores prior
visibility. Each model remembers its current domain. **Show only selected model**
also hides other atomic structures. Click a domain or use **Previous / Next** to select
its residues. **Auto-center on domain**, **Focus domain**, and **Restore view**
control the camera. Switching models alone keeps the current view and selection.
The list is organized by domain/feature signature, including for a single model.
Each expandable parent shows its copy count and collates occurrences across chains
and models. Selecting a parent selects all listed copies; expand it to select an
individual occurrence and inspect its scores. **Previous / Next** moves between
groups when a parent is selected, or between occurrences when a child is selected.
The model dropdown and **Nearby / Selected** filters apply before grouping, so
copy counts describe the current view. A group spanning multiple models has no
single model for **Show only selected model**; select a child or choose a model.

Use **Nearby** to identify annotated domains around the current center of rotation
(default radius **6 Å**, measured to the nearest domain atom in visible models).
Use **Selected** to list every annotation overlapping the selected residues.
Both tabs update automatically and respect the Model filter. Automatic updates
show details without changing the selection or camera; clicking a row uses the
usual domain selection and focus controls. **All domains** restores the full list.

Details prefetch in the background across analyzed models: available function
descriptions, InterPro/GO annotations, database links, and example PDB structures.
The selected entry takes priority. Successful remote lookups are cached on disk
for 14 days; descriptions and formatted details are reused within the session.
For Foldseek hits with localized Pfam support, available InterPro/Pfam
descriptions and example PDBs take priority. Unverified ECOD family-level links
remain candidate links; the displayed region keeps its structural classification.
Local library descriptions provide a
fallback, and the matched structural reference remains linked.
Keep `domain_exemplars.py` beside the main helper to enable the **Structures**
controls where source coordinates and domain boundaries are known.
This includes PDB-backed ECOD matches and AlphaFold-derived ECOD references
identified as `<UniProt>_F1_nD<number>`, including their discontinuous boundaries.
ECOD domains use the exact **Download Domain PDB** endpoint from their ECOD page;
the full PDB or AlphaFold model is fitted to that domain for hidden context.
AlphaFold source identity, model coverage and chain are checked against its API;
unknown prediction identifiers are not guessed. The current AlphaFold model is
retrieved, which may differ in version from ECOD's original template.
Other ECOD identifiers still offer their exact domain download; when the full
source cannot be identified, only the domain is loaded and context is marked unavailable.
**Fetch & align domain** downloads the structure in the background,
fits its mapped domain to the selected occurrence with ChimeraX Matchmaker,
and shows a separate domain model. The full structure is also loaded, hidden
by default, with exactly the same alignment. **Show full context / Hide context**
switches between them; the domain is gold and its surrounding structure gray.
Fetching or realigning frames the visible exemplar; switching context refits
the view and clipping to the domain or full structure and centers rotation there.
Each has its own model entry. **Realign domain** updates the fit after moving
the query. Alignment uses mapped domain residues, preserves discontinuous
regions, and reports the fitted residue count and RMSD. Coordinate downloads
are cached for 14 days. The query's coordinates and colors are unchanged.
PDB exemplars include DOI/PubMed links to their deposited primary structure paper
when available. Publications load separately in the background and are cached per
PDB entry across families and models.
These describe the family; functions and specificity can differ between members.
**Reload details** bypasses the description cache without repeating the search.
Example PDBs require a matching domain mapped onto a crystallized chain;
membership of the parent protein or a shared paper is insufficient. Details show
the matched family, chain, and residue range. Older unverified structure lists
are refreshed in the background. Citation cleanup also repairs separators left
in cached descriptions. After reloading the helper, `domainbrowser` replaces an
open panel from the previous script instance without rerunning domain searches.

**Sequence features.** Domain searches also run SEG low-complexity detection and
PROSITE pattern matching by default, using Python's standard library. Local scans
run alongside the domain lookup, share cancellation, and apply together in one Undo
action. Use `features false` to skip these additional scans, or
`domains #1 source features` to run feature analysis alone. The current browser
view is preserved when domains and features are analyzed together.
The feature helper downloads the public PROSITE definitions once (about 26 MB), then works offline;
no sequences are sent for local scans. Use `prosite_database "/path/prosite.dat"`
to supply definitions offline. If the download fails, SEG results are still returned.
The original SIB definition file and its license notice are retained unchanged.

**Transmembrane features.** `domains #1 source tmsegments` predicts membrane
segments with the existing DeepTMHMM/BioLib backend and preserves domain colors.
Use `results "/path/predicted_topologies.3line"` to import predictions offline, or
`python "/path/to/biolib/venv/bin/python"` for another configured client. Both
`tmsegments` and this source add separate **TM helix** and **TM beta strand** groups
to **Sequence features**, with per-chain segment numbers, selectable copies, and
**Nearby / Selected** filtering. Later domain/local-feature searches retain them.
DeepTMHMM runs only when explicitly requested; ordinary domain searches do not
submit remote topology predictions. Both routes reuse the same session cache and
complete-chain mapping, and support Undo/Redo. Either `domains cancel` or
`tmsegments cancel` stops a prediction started through `domains`.
See [DeepTMHMM setup and imports](#transmembrane-helices-and-beta-strands).

Choose **Sequence features** or **All annotations** in the browser to see motifs,
sites, low-complexity regions, and available Pfam disorder/coiled-coil annotations.
Pfam motifs and InterPro sites appear automatically after their respective searches,
without an additional search. An older InterPro session cache is refreshed once to
recover previously discarded sites. **Nearby** and **Selected** work in every view.
Features can overlap domains and each other; selecting one preserves domain colors.
Separate `sequence_feature_names` and `sequence_feature_accessions` residue attributes
support selection, and results/attributes support Undo/Redo. Running feature analysis
retains domain results; subsequent domain searches retain independently scanned features.

Local PROSITE matches show the pattern, matched sequence, description and documentation
link. They are candidate sequence signatures, not confirmed functional sites. Common,
nonspecific signatures are hidden unless `prosite_common true` is set. SEG uses the
standard 12/2.2/2.5 settings with probability-based boundary trimming; low complexity
does not imply disorder. Repeated and overlapping motif occurrences are retained.

Optional `prosite_profiles true` also sends distinct sequences to the
[ScanProsite API](https://prosite.expasy.org/scanprosite/scanprosite_doc.html) for scored
profile matches at its standard cutoff. These appear as **Profile-supported region**;
individual catalytic/binding residues are shown only when supplied as site annotations
by InterPro. Online failure retains the local results. Local results and successful
online results, including empty results, are cached on disk for 14 days by sequence,
definition checksum, algorithm and options. Automatic features reuse this cache
even when `refresh true` repeats a domain search; use
`domains #1 source features refresh true` to refresh feature results and public
definitions. Each identical sequence is analyzed once across models. A feature
backend failure leaves valid domain results available and is reported in the Log.
See [PROSITE](https://prosite.expasy.org/) and its
[license](https://prosite.expasy.org/prosite_license.html).
PROSITE data is CC BY-NC-ND 4.0; commercial use requires a separate license.
It is downloaded unchanged on demand and is not part of our structural library.
Without the necessary permission, use `features false` to avoid the automatic
PROSITE feature pass. Do not include its definitions or caches in an unrestricted
database release. See [redistribution notes](THIRD_PARTY.md).

| Option | Purpose | Default |
| --- | --- | --- |
| `source` | `pfam`, `interpro`, local `foldseek`, local `features`, or DeepTMHMM `tmsegments` | `pfam` |
| `results` | Import a DeepTMHMM `.3line` file; applies to `source tmsegments` | Predict or reuse cached topology |
| `python` | BioLib client Python; applies to `source tmsegments` | Existing DeepTMHMM environment |
| `features` | Also run local SEG and PROSITE scans with any domain source | `true` |
| `prosite_profiles` | Add online ScanProsite profile matches to enabled feature scans | `false` |
| `prosite_common` | Include common/nonspecific PROSITE signatures | `false` |
| `prosite_database` | Existing `prosite.dat` file for feature scans | Download and cache public definitions |
| `library` | Reference-library directory; applies to `foldseek` | Configured local library, otherwise `cdd-pilot` |
| `regions` | Quoted sequence ranges to search separately; applies to `foldseek` | Requested residues as one query per chain |
| `pfam_guided` | Use Pfam boundaries to guide local Foldseek searches | `false` |
| `pfam_backend` | `auto`, `local`, or `online`; applies to Pfam and Pfam-guided Foldseek | `auto` |
| `pfam_database` | Local `Pfam-A.hmm` file or containing directory | Automatic discovery |
| `deduplicate` | Share Foldseek searches across identical sequences with the same modeled query positions | `true` |
| `color` | Color domain atoms and cartoons; sequence features preserve domain colors | `true` |
| `include_overlaps` | Include competing Pfam assignments | `false` |
| `include_tentative` | Show reported below-cutoff Pfam candidates with E-value ≤ 0.01 | `true` for Pfam |
| `refresh` | Repeat the search instead of using session-cached results | `false` |
| `gui` | Open or refresh the browser | `true` |
| `wait` | Wait for completion before subsequent script commands | `false` |
| `timeout` | Maximum lookup time in seconds | `600`; `1200` for `source tmsegments` |

Pfam's clan competition can hide a weaker assignment even for a one-residue
alignment overlap within the same clan. Use `include_overlaps true` to display
significant competing hits. Pfam includes reported below-cutoff candidates with
independent domain E-value ≤ 0.01 by default, labeled **tentative** in the browser
and log. Use `include_tentative false` for curated matches only.
Details show their score and the domain
gathering cutoff. A tentative hit is hidden if at least 50% of its sequence
footprint overlaps the combined footprints of curated matches. Small boundary
overlaps are retained. This rule also applies with `include_overlaps true`.
Curated assignments take priority for coloring remaining overlaps. Both options
can be combined, and toggling either reuses cached hits.
The service may omit weaker candidates entirely. InterPro retains all returned
locations. Foldseek retains separate occurrences of a family and searches
remaining stretches for repeat matches. Structural assignments sharing at least
**70% of both aligned residue footprints** collapse into one browser row, choosing
the lowest E-value, then highest bit score. Shorter, strongly contained matches
also group when they have comparable size or share an explicit Pfam family or
ECOD homology group with the primary. The row's details retain the alternative
names, accessions, scores, and links. Small overlaps, separate repeats, and
unrelated small inserted domains remain separate. Primary footprints
are not expanded to cover the alternatives. Rerunning the same command applies
this grouping to session-cached hits without repeating the structural search.

Searches run in the background in graphical sessions; use `wait true` in scripts
that depend on the result. Foldseek shows an elapsed-seconds counter in the
status bar, updated about once per second, with the current search pass and
query-region count. Completion posts a separate status and Log message with
elapsed time, match count, and model count. Cancellation, timeout, and failure
replace the running status; cached results are identified explicitly.
Colors, attributes, and browser results support
Undo/Redo. Reanalyzing a chain replaces its results while retaining other models.
Results are cached per source for the current session; rerun `domains` after
loading a saved session to rebuild the browser.

Residue attributes `domain_name` and `domain_accession` hold the primary match;
`domain_names` and `domain_accessions` retain overlapping matches separated by
` | `. For example: `select #1 & ::domain_accession=PF00071`.
`domain_tentative` records whether the primary assignment is tentative; select
those residues with `select #1 & ::domain_tentative=true`.


## Transmembrane Helices and Beta Strands

`tmsegments` runs [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM) on BioLib and
annotates both predicted **TM helices** (`M`, orange) and **TM beta strands**
(`B`, cyan). Other residues are gray. Both atom and cartoon colors are updated;
use `color false` to assign attributes while retaining existing colors. Predictions
do not change the model's secondary-structure assignments or coordinates.

```chimerax
tmsegments #1
tmsegments #1/A
tmsegments #1-3 wait true
tmsegments #1 color false
tmsegments #1 gui false
tmsegments #1 refresh true timeout 1800
tmsegments cancel
```

The atom specification identifies **complete chains** to analyze: touching one
residue requests its whole protein chain. Full chain sequences, where available,
are sent to BioLib; coordinates stay local. Identical sequences are predicted
once and results are cached for the current ChimeraX session. `refresh true`
requests fresh predictions. GUI submissions run in the background; `wait true`
waits for completion and is useful in scripts. Headless sessions always wait.
While a remote job is pending, the status bar shows a running message with an
elapsed-seconds counter, updated about once per second. This includes time queued
on BioLib. Completion replaces it with the elapsed time and counts of TM helices
and TM beta strands. Cancellation, timeout, and failure also update the status bar.
The default timeout is 1200 seconds. Cancellation/timeout stops the local client;
an already submitted remote job may continue on BioLib.

### BioLib setup

DeepTMHMM has separate academic/commercial terms; obtain the applicable
[DTU/BioLib license](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/)
for commercial use. The MIT license of the BioLib Python client does not grant
rights to the predictor. No predictor weights or containers are distributed here.

The command is included in `chimerax_trimmings.py`; no extra startup script is
needed. The BioLib client runs in a separate Python environment. On macOS, create
the default environment in Terminal (using an installed Python 3):

```bash
python3 -m venv ~/.local/share/chimerax-trimmings/deeptmhmm
~/.local/share/chimerax-trimmings/deeptmhmm/bin/python -m pip install pybiolib
```

To use another environment, pass its Python executable explicitly:

```chimerax
tmsegments #1 python /absolute/path/to/venv/bin/python
```

Alternatively set `CHIMERAX_DEEPTMHMM_PYTHON` in the environment from which
ChimeraX is launched. If BioLib requires authentication, log in once from Terminal
with that environment's `biolib login`, or supply `BIOLIB_TOKEN` in ChimeraX's
launch environment. Tokens are not command arguments. See the
[BioLib Python documentation](https://biolib.com/docs/using-applications/python/).
The client interface was checked against pybiolib 1.4.647.

### Import and inspect predictions

An existing DeepTMHMM `predicted_topologies.3line` file can be imported without
installing the BioLib client or submitting a remote job:

```chimerax
tmsegments #1 results /absolute/path/predicted_topologies.3line
select #1 & ::tm_type=helix
select #1 & ::tm_type=strand
select #1 & ::tm_topology=S
select #1/A & ::tm_segment=2
```

The file must contain an **exact sequence match** for each requested chain;
identifiers need not match model or chain names. Extra sequences in a batch file
are allowed. Imported predictions take precedence over cached remote results
for that run and do not replace the remote cache.

| Residue attribute | Meaning |
| --- | --- |
| `tm_topology` | Original DeepTMHMM label: `M`, `B`, `S`, `I`, `O`, or `P` |
| `tm_type` | `helix`, `strand`, or an empty string for other residues |
| `tm_segment` | One-based TM segment number within the chain; `0` outside TM segments |

Signal peptides (`S`) retain their own topology label and are excluded from TM
segment counts. `I`, `O`, and `P` retain DeepTMHMM's inside/outside/periplasm labels;
these are sequence predictions, not an estimated membrane plane in 3D.

The Log shows a clickable list of both segment types, including sequence ranges
and their corresponding structure residue specifications. Missing coordinates,
residue-number offsets, and insertion codes are respected. A segment with no
modeled residues is still reported. Unknown/nonstandard sequence letters are
rejected rather than silently deleted or replaced. If the full polymer sequence
is unavailable, the helper warns before using ChimeraX's inferred chain sequence.

Results also open the **Sequence features** browser view, grouping helices and
beta strands separately across models and chains. These browser records and
`sequence_feature_names` / `sequence_feature_accessions` update together with
the topology attributes in one Undo action. `gui false` suppresses opening the
browser. Segments without modeled residues remain in the Log. `domains #1 source
tmsegments` offers the same prediction/import workflow while retaining domain colors;
the original `tmsegments` command keeps its optional membrane coloring.

Annotations and coloring support Undo/Redo. Reanalyzing a chain replaces its TM
annotations, including clearing earlier assignments when no TM segments are
predicted; other chains' records remain available. Predictions for a chain that
was deleted or whose sequence/residue mapping changed during submission are
skipped. The combined records are available as `session._trimmings_tm_records`.
Attributes/records are session annotations; this helper does not add a custom
mmCIF export or persist the prediction cache across ChimeraX sessions.

## Exposed Hydrophobic Residues

`exposedhydrophobics` replaces the selection with whole protein residues whose
relative sidechain solvent accessibility exceeds a fractional cutoff.

```chimerax
exposedhydrophobics #1
exposedhydrophobics #1/A minimum_exposure 0.4
exposedhydrophobics #1 residue_names LEU,ILE,VAL,PHE probe_radius 1.4
exposedhydrophobics #1 minimum_exposure 0.25 minimum_area 10
```

| Argument | Meaning | Default |
| --- | --- | --- |
| Atom specification | Candidate residues, e.g. `#1`, `#1/A`, or `sel` | Required |
| `minimum_exposure` | Relative sidechain exposure must exceed this fraction, from 0 to 1 | `0.25` (25%) |
| `minimum_area` | Optional additional sidechain SASA cutoff, in square Angstroms | `0` |
| `residue_names` | Comma-separated residue names to consider | `ALA,VAL,ILE,LEU,MET,PHE,TRP,PRO` |
| `probe_radius` | Solvent probe radius in Angstroms | `1.4` |

The calculation uses all non-water atoms in each candidate's model, including
hidden chains, nucleic acids, ligands, and ions. Restricting candidates to one
chain does not remove burial by other chains in that model. Separate models are
calculated independently, so overlapping alternative structures do not bury each
other. Existing atom radii and any hydrogens present are used. Sidechain area is
summed using ChimeraX's `sideonly` definition, excluding backbone atoms and CA.
Glycine therefore has zero sidechain area.

Relative exposure is the sidechain SASA in the full model divided by its sidechain
SASA with only that residue present. The reference retains the residue's own
backbone and conformation and uses the same atom radii, hydrogens, and probe radius
as the full-model calculation. This compares exposed fractions rather than imposing
the same absolute area on large and small sidechains. It is an isolated-residue
reference, not a published maximum-accessibility table or tripeptide reference.
Residues with zero reference area have zero relative exposure. Incomplete
sidechains are normalized using only the atoms present, so a highly exposed
fragment does not imply that the missing sidechain would also be exposed.

The command assigns three attributes to every protein residue in the models it
calculates: `sidechain_sasa` (square Angstroms), `sidechain_sasa_reference` (square
Angstroms), and `relative_sidechain_sasa` (0 to 1). They support related selections
without recalculating, for example:

```chimerax
select #1 & :ARG,LYS,ASP,GLU & ::relative_sidechain_sasa>0.25
```

Rerun `exposedhydrophobics` after changing coordinates, radii, model composition,
or the probe radius to refresh those values. It does not create a surface or
change model visibility, coloring, or coordinates. Selection changes are undoable.

## Binder-Accessible Hydrophobic Patches

`hydrophobicpatches` screens for exposed hydrophobic clusters that a globular
binder, approximated as a rigid sphere, could contact. It selects the reported
patches and, by default, gives each patch a distinct atom/cartoon color.

```chimerax
hydrophobicpatches #1 binder_radius 15
hydrophobicpatches #1/A binder_mass 15 context #2
hydrophobicpatches #1 binder_radius 15 samples 128 highlight false
select #1 & ::hydrophobic_patch_rank=1
```

Supply exactly one of `binder_radius` (Angstroms) and `binder_mass` (kDa).
Mass is converted to an equivalent dry sphere at an assumed protein density of
1.35 g/cm^3; 15 kDa corresponds to about 16.4 Angstroms radius. Use a radius when
you have a better estimate of the binder's effective steric size.

| Option | Meaning | Default |
| --- | --- | --- |
| Atom specification | Candidate residues; the full containing model remains the steric environment | Required |
| `binder_radius` | Binder sphere radius in Angstroms | Supply this or `binder_mass` |
| `binder_mass` | Approximate globular binder mass in kDa | Supply this or `binder_radius` |
| `context` | Additional blocker atoms from separate models, e.g. glycans or assembly partners | None |
| `minimum_exposure` | Minimum relative sidechain exposure with the solvent probe | `0.25` |
| `probe_radius` | Solvent probe radius in Angstroms, distinct from the binder radius | `1.4` |
| `minimum_residues` | Minimum hydrophobic residues in a patch | `3` |
| `contact_distance` | Maximum atom-surface to binder-surface gap in Angstroms | `3` |
| `link_distance` | Maximum distance between exposed sidechain atom centers linking residues | `6` |
| `samples` | Deterministic sphere-placement directions per exposed heavy sidechain atom | `64` |
| `max_patches` | Maximum number of distinct, non-overlapping patches reported | `10` |
| `approach` | Require a straight, collision-free approach from outside | `true` |
| `highlight` | Color patch atoms and cartoons | `true` |
| `gui` | Automatically open or refresh the patch browser in GUI sessions | `true` |
| `residue_names` | Hydrophobic residue names | `ALA,VAL,ILE,LEU,MET,PHE,TRP,PRO` |

The screen first applies the same relative sidechain-exposure definition as
`exposedhydrophobics`, using the solvent probe (1.4 Angstroms by default). It then samples sphere
placements tangent to exposed sidechain atoms. Every reported patch must have
one shared clash-free placement near all its reported residues, and its residues
must form a connected cluster at `link_distance`. With `approach true`, the sphere
must also have an unobstructed straight outward path to infinity. This is stricter
than testing whether individual residues are independently accessible.

All non-water atoms in the containing model contribute to exposure and steric
occlusion, including hidden chains, ligands, and modeled glycans. Each candidate
model is screened independently; `context` adds blockers from separate models.
Patches may span chains within a model, but are not combined across models.

Patches are ranked by the summed water-probe SASA of hydrophobic sidechain atoms
near the sampled binder sphere. This is nearby exposed area, not predicted buried
interface area or binding energy. The highest-ranked patch is kept first, then
lower-ranked overlapping alternatives are omitted. Residue attribute
`hydrophobic_patch_rank` is the rank, or zero for unassigned protein residues in
the analyzed models. Placement centers and approach directions in scene coordinates
are available in `session._trimmings_hydrophobic_patches`.

Use these as candidate epitopes for inspection, not predictions of affinity,
specificity, or binder design success. A sphere omits binder shape, sidechain
flexibility, missing glycans, and unmodeled membrane/assembly context. Finite
direction sampling can miss accessible placements; increase `samples` to refine
the screen. A failed straight-approach test does not rule out a curved approach.
Incomplete target sidechains retain the exposure limitations described above.

Color and selection changes are undoable together; undo a highlight before
comparing another set of parameters. `highlight false` leaves colors untouched.
The command uses NumPy and SciPy supplied with ChimeraX, creates no model copies
or probe models, and leaves coordinates, map contours, and display masks unchanged.

### Patch Browser

In a graphical ChimeraX session, a search automatically opens a dockable
**Hydrophobic Patches** panel. Its list-and-navigation layout follows ResNotes,
but ResNotes does not need to be installed.

- Click a patch in the ranked list, or use **Previous** / **Next**, to select it.
  Navigation wraps at either end.
- **Auto-center on patch** fits the patch and adjusts local clipping. Uncheck it
  to change selection while keeping the current camera.
- **View from binder approach** orients focused views along the sampled approach
  direction. Uncheck it to retain the current orientation when centering.
- **Focus patch** centers the current patch manually.
- **Restore view** restores the camera, clipping, and rotation center captured
  when the panel opened. It does not restore or change model positions.
- **Refresh list** loads the latest search results. New searches refresh the
  existing panel, and closed-model results are removed automatically.

The details show the residue specification, nearby hydrophobic area, and binder
sphere radius. Focusing follows a moved target's current transform; rescreen after
structural edits or changes in the target's position relative to blocker models.
Hidden models are left hidden, and the panel does not change display masks.

Reopen the latest results with `patchbrowser`. Suppress automatic opening with:

```chimerax
hydrophobicpatches #1 binder_mass 15 gui false
```

Without a GUI, searches remain command-only. Browser state is not saved in
session files; rerun the search after reopening a session.

## Active and Frozen Models

The opt-in `modelmotion` commands provide Chimera-style activation independently
of atom/model selection. Structures, maps, and independent surfaces are supported.
Groups expand to their movable members; surfaces owned by a structure or map move
with that owner. Specify the owner rather than its surface submodel.

```chimerax
modelmotion active #1
modelmotion unfreeze #3
modelmotion freeze #1
modelmotion status
modelmotion all
modelmotion off
```

`active` replaces the movement set, `unfreeze` adds models, and `freeze` removes
models. If `freeze` is the first command, all other current models remain active.
Newly opened models are active by default; existing frozen models stay frozen.
An empty active set prevents rotation/translation; clearing selection has no effect.

Mouse and trackpad input work simultaneously and use the same active/frozen model
set. With the Trimmings startup bindings and ChimeraX's standard middle-button and
wheel bindings, mouse controls are:

| Mouse input | Behavior |
| --- | --- |
| Left drag | Rotate active models together |
| Middle-button drag | Translate active models in the screen plane |
| Control + middle-button drag | Translate active models along screen Z |
| Scroll wheel | Keep the existing whole-view zoom binding |
| Right drag | Keep the Trimmings clipping control |
| Alt + right drag | Keep the Trimmings contour control |

Custom mouse bindings are respected: existing model rotation/translation bindings
are redirected to the active set when `modelmotion` is enabled. Selected-model
movement bindings then use the active set independently of selection.

| Trackpad gesture | Movement |
| --- | --- |
| Two-finger swipe | Rotate the active models together |
| Two-finger twist | Rotate active models around the screen-normal axis |
| Three-finger swipe | Translate active models in the screen plane |
| Control + three-finger swipe | Translate active models along screen Z |
| Pinch | Keep the existing zoom binding |

Trackpad gestures must be enabled in ChimeraX and available to it from macOS.
Mouse and trackpad rotation use the current ChimeraX center of rotation, including
the pivot set by `cofron`, `cofroff`, `cofr sel`, or an explicit `cofr` position.
Mouse bindings that currently
rotate/translate models are redirected to the same active set; selection, clipping,
and contour bindings retain their existing behavior. Trackpad movement is undoable
after a short pause. `modelmotion off` restores the previous mouse and trackpad
bindings changed on activation.
Freezing applies to these movement modes; commands such as `fitmap` can still move
models. Activation is temporary session state and is not saved in session files.

## Keybindings

- `F1`: `prevmodel`
- `F2`: `nextmodel`
- `F3`: `centersel`
- `F4`: `view all`
- `F5`: `togglemaps`
- `F6`: `togglemodels`
- `F7`: `cyclelighting`
- `F8`: `cyclemodeldisplay`
- `F9`: `prevmap`
- `F10`: `nextmap`
- `F11`: `lowermapthreshold`
- `F12`: `raisemapthreshold`

## Button Panel

The startup file creates a `Shortcuts` button panel with quick access to commonly used tools and commands.

It currently includes buttons for:

- `Vol_Viewer`: open the Volume Viewer tool
- `Model_Panel`: open the Model Panel
- `Log`: open the Log tool
- `default_disp`: restore the default ribbon/cartoon molecular display
- `map_sphere`: zone maps to a 15 A sphere around the center of rotation
- `map_unsphere`: remove map zoning
- `cofron`: show the center-of-rotation pivot marker
- `cofroff`: hide the center-of-rotation pivot marker while keeping the CoFR centered on the view
- `cootmode`: apply the dark-background Coot-style display preset
- `mark_cofr`: place a marker at the current center of rotation
- `togglemaps`: hide shown maps, then restore the same set on the next use
- `togglemodels`: hide shown atomic models, then restore the same set on the next use
- `reset_mouse_and_help`: restore the default mouse bindings and open the mouse-mode help page
- `previous_model_map`: show the previous unskipped atomic model or volume map
- `next_model_map`: show the next unskipped atomic model or volume map

## Alias Reference

### `cofron`

Centers the rotation point on the current view and shows the pivot marker.

Usage: `cofron`

```chimerax
alias cofron cofr centerofview showpivot 7,0.25
alias usage cofron synopsis "Center the rotation point on the view and show the pivot"
```

### `cofroff`

Recenters the rotation point on the current view and hides the pivot marker.

Usage: `cofroff`

```chimerax
alias cofroff cofr centerofview showpivot false
alias usage cofroff synopsis "Recenter the rotation point on the view and hide the pivot marker"
```

### `symclip`

Sets symmetric clip planes around the center of rotation.

Usage: `symclip <half-distance>`

```chimerax
alias symclip cofr centerofview; clip near -$1 far $1 position cofr
alias usage symclip synopsis "Set symmetric clip planes around the center of rotation" $1 "half-distance"
```

### `modelsel`

Selects atomic models using the `##num_residues` selector.

Usage: `modelsel`

```chimerax
alias modelsel select ##num_residues
alias usage modelsel synopsis "Select atomic models"
```

### `mapsel`

Selects non-atomic models such as maps using the inverse selector.

Usage: `mapsel`

```chimerax
alias mapsel select ~##num_residues
alias usage mapsel synopsis "Select non-atomic models such as maps"
```

### `sidesel`

Selects protein sidechains and nucleic-acid bases for the given atom specification.

Usage: `sidesel <atom-spec>`

```chimerax
alias sidesel select sidechain & $1 | ~backbone & nucleic & $1
alias usage sidesel synopsis "Select protein sidechains and nucleic-acid bases" $1 "atom-spec"
```

### `sidedisp`

Displays protein sidechains and nucleic-acid bases for the given atom specification.

Usage: `sidedisp <atom-spec>`

```chimerax
alias sidedisp select sidechain & $1 | ~backbone & nucleic & $1; display sel
alias usage sidedisp synopsis "Display protein sidechains and nucleic-acid bases" $1 "atom-spec"
```

### `cootmode`

Applies a dark-background Coot-style display preset for model building.

Usage: `cootmode`

```chimerax
alias cootmode set bgColor black; surface cap false; surface style solid; nucleotides #* atoms; lighting flat; graphics silhouettes false; style stick; ~rib; color ##num_residues gold; color byhet ; disp;  ~disp @H*; style ions sphere; style solvent ball; size ballscale 0.2;  size stickradius 0.07; transparency 70; cofr centerofview; clip near -10 far 10 position cofr; color ~##num_residues cornflower blue
alias usage cootmode synopsis "Apply the dark-background Coot-style display preset"
```

### `cootmode_white`

Applies a white-background variant of `cootmode`.

Usage: `cootmode_white`

```chimerax
alias cootmode_white set bgColor white; surface cap false; surface style solid; nucleotides #* atoms; lighting flat; graphics silhouettes false; style stick; ~rib; color ##num_residues orange; color byhet ; disp;  ~disp @H*; style ions sphere; style solvent ball; size ballscale 0.2;  size stickradius 0.07; transparency 85; cofr centerofview; clip near -10 far 10 position cofr; color ~##num_residues medium blue
alias usage cootmode_white synopsis "Apply the white-background Coot-style display preset"
```

### `cootmode_mesh`

Applies a mesh-map variant of `cootmode`.

Usage: `cootmode_mesh`

```chimerax
alias cootmode_mesh surface cap false; surface style mesh; lighting flat; nucleotides #* atoms; graphics silhouettes false; style stick; ~rib; color ##num_residues gold; color byhet ; disp;  ~disp @H*; style solvent ball; style ions sphere; size ballscale 0.2;  size stickradius 0.07; cofr centerofview; clip near -10 far 10 position cofr; color ~##num_residues #3d60ffff; transparency 50
alias usage cootmode_mesh synopsis "Apply the mesh Coot-style display preset"
```

### `ca_and_sidechains`

Shows a CA trace for proteins or phosphate trace for nucleic acids, plus protein sidechains and nucleic-acid bases.

Usage: `ca_and_sidechains <model-spec>`

```chimerax
alias ca_and_sidechains ~rib $1; ~surf $1; ~disp $1; disp @CA&protein&$1; disp @P&nucleic&$1; style $1 stick; disp sidechain&$1; disp ~backbone&nucleic&$1; size stickradius 0.1; size pseudobondradius 0.1
alias usage ca_and_sidechains synopsis "Show CA or phosphate trace plus sidechains or bases" $1 "model-spec"
```

### `ca_trace`

Shows a CA trace for proteins and a phosphate trace for nucleic acids.

Usage: `ca_trace <model-spec>`

```chimerax
alias ca_trace ~rib $1; ~surf $1; ~disp $1; disp @CA&protein&$1; disp @P&nucleic&$1; style $1 stick; size stickradius 0.1; size pseudobondradius 0.1
alias usage ca_trace synopsis "Show a CA or phosphate trace" $1 "model-spec"
```

### `expand_sym`

Expands an asymmetric unit to biological assembly `1` as a new model and shows the result as ribbon copies. This is most useful for mmCIF structures where assembly `1` corresponds to the full biological unit.

Usage: `expand_sym <model-spec>`

Note: this alias assumes the desired biological assembly is `1`.

```chimerax
alias expand_sym ~disp $1; rib $1; sym $1 assembly 1 copies false newModel true
alias usage expand_sym synopsis "Expand an asymmetric unit to biological assembly 1 as ribbon copies" $1 "model-spec"
```

### `map_sphere_15`

Zones maps to a 15 A sphere around the center of rotation.

Usage: `map_sphere_15`

```chimerax
alias map_sphere_15 volume unzone ~##num_residues; sel; close #10000; marker #10000 position cofr; sel ~sel; volume zone ~##num_residues nearAtoms sel minimalBounds true range 15; close #10000
alias usage map_sphere_15 synopsis "Zone maps to a 15 A sphere around the center of rotation"
```

### `cubic_map`

Creates a spherical mask map on a cubic grid.

Usage: `cubic_map <sphere-radius>`

```chimerax
alias cubic_map shape sphere radius $1 modelid #10000; volume onesmask #10000 spacing 1 border -0.5
alias usage cubic_map synopsis "Create a spherical mask map on a cubic grid" $1 "sphere radius"
```

### `map_unsphere`

Removes zoning from all maps.

Usage: `map_unsphere`

```chimerax
alias map_unsphere volume unzone ~##num_residues
alias usage map_unsphere synopsis "Remove map zoning from all maps"
```

### `fit_by_chain`

Splits a model by chain, fits each chain independently into a map, then recombines the result.

Usage: `fit_by_chain <model-spec> <map-spec>`

```chimerax
alias fit_by_chain split $1; fitmap $1.* inmap $2 eachmodel true; combine $1.* close true
alias usage fit_by_chain synopsis "Split a model by chain, fit each chain into a map, and recombine" $1 "model-spec" $2 "map-spec"
```

### `rainbow_new`

Colors models with the `RdYlBu-5` rainbow palette.

Usage: `rainbow_new`

```chimerax
alias rainbow_new rainbow palette RdYlBu-5
alias usage rainbow_new synopsis "Color models by rainbow using the RdYlBu-5 palette"
```

### `default_mol_display`

Restores a default ribbon/cartoon molecular display with chain rainbow coloring.

Usage: `default_mol_display`

```chimerax
alias default_mol_display ~disp; rib; rainbow chain palette RdYlBu-5; lighting soft
alias usage default_mol_display synopsis "Restore the default cartoon molecular display"
```

### `local_fitmap`

Fits a model or map into a local zone cut from a target map around the center of rotation.

Usage: `local_fitmap <model-or-map> <target-map> <zone-radius>`

```chimerax
alias local_fitmap ~sel; close #10000-10001; marker #10000 position cofr; sel #10000; volume zone $2 nearAtoms sel range $3 newMap true modelid 10001 minimalbounds true; fitmap $1 inmap #10001 eachmodel true; close #10001; close sel; show $2
alias usage local_fitmap synopsis "Fit a model or map into a local zone around the center of rotation" $1 "model or map to fit" $2 "target map" $3 "zone radius"
```

### `hidemodels`

Hides all atomic models.

Usage: `hidemodels`

```chimerax
alias hidemodels hide ##num_residues models
alias usage hidemodels synopsis "Hide all atomic models"
```

### `showmodels`

Shows all atomic models.

Usage: `showmodels`

```chimerax
alias showmodels show ##num_residues models
alias usage showmodels synopsis "Show all atomic models"
```

### `caps_off`

Turns clipped surface caps off.

Usage: `caps_off`

```chimerax
alias caps_off surface cap false
alias usage caps_off synopsis "Turn clipped surface caps off"
```

### `caps_on`

Turns clipped surface caps on.

Usage: `caps_on`

```chimerax
alias caps_on surface cap true
alias usage caps_on synopsis "Turn clipped surface caps on"
```

### `open_vseries`

Opens a file browser with volume-series mode enabled.

Usage: `open_vseries`

```chimerax
alias open_vseries "open browse vseries true"
alias usage open_vseries synopsis "Open a file browser with volume series enabled"
```

### `selbetween`

Selects residues between the current selection endpoints.

Usage: `selbetween`

```chimerax
alias selbetween ks ri
alias usage selbetween synopsis "Select residues between the current selection endpoints"
```

### `helix`

Marks residues as helix.

Usage: `helix <residue-spec>`

```chimerax
alias helix setattr $1 res is_helix true
alias usage helix synopsis "Mark residues as helix" $1 "residue-spec"
```

### `strand`

Marks residues as strand.

Usage: `strand <residue-spec>`

```chimerax
alias strand setattr $1 res is_strand true
alias usage strand synopsis "Mark residues as strand" $1 "residue-spec"
```

### `coil`

Clears helix and strand assignment for residues.

Usage: `coil <residue-spec>`

```chimerax
alias coil setattr $1 res is_strand false; setattr $1 res is_helix false
alias usage coil synopsis "Clear helix and strand assignment for residues" $1 "residue-spec"
```

### `rock_movie`

Records a rocking movie to `~/Desktop/rock_movie.mp4`.

Usage: `rock_movie`

```chimerax
alias rock_movie cofr showpivot false; movie record; rock y 50; wait 600; movie encode ~/Desktop/rock_movie.mp4; stop
alias usage rock_movie synopsis "Record a rocking movie and save it to the default Desktop file"
```

### `centersel`

Centers on the current selection and applies local clipping around it.

Usage: `centersel`

```chimerax
alias centersel cofr sel; clip near 10 position cofr; clip far -10 position cofr; cofr centerOfView showpivot true; view sel
alias usage centersel synopsis "Center the view on the selection and apply local clipping"
```

### `volume_project`

Displays a map as a projected image volume.

Usage: `volume_project <map-spec>`

```chimerax
alias volume_project set bgColor black; volume $1 step 1 sd_level 0,0 sd_level 50,1 color white color white style image projection_mode auto maximum_intensity_projection true bt_correction true linear_interpolation true; lighting depth_cue false
alias usage volume_project synopsis "Display a map as a maximum-intensity projection" $1 "map-spec"
```

### `split_diff_map`

Replaces a map with a scaled positive/negative mesh difference map.

Usage: `split_diff_map <input-map> <output-model-id>`

```chimerax
alias split_diff_map volume scale $1 rms 1 modelId $2; close $1; volume $2 capFaces false meshLighting false squareMesh false level -3 color #da1200000000 level 3 color #0000bda00000; lighting depthCue true; volume $2 style mesh
alias usage split_diff_map synopsis "Replace a map with a scaled positive and negative mesh difference map" $1 "input map" $2 "output model-id"
```

### `local_diff_map`

Builds a local difference map from two maps around a model.

Usage: `local_diff_map <model-spec> <first-map> <second-map>`

Depends on: `split_diff_map`

```chimerax
alias local_diff_map view name tmp; close #1001,1002,1003,1004,1005; fitmap $1 inMap $2; vol zone $2 near $1 range 5 newMap true modelId 1000 minimalBounds true; fitmap $1 inMap $3; vol zone $3 near $1 range 5 newMap true modelId 1001 minimalBounds true; volume scale #1000 sd 0.1 modelId 1002; volume scale #1001 sd 0.1 modelId 1003; fitmap #1003 inMap #1002; volume resample #1003 onGrid #1002 modelId 1004; volume subtract #1002 #1004 modelId 1005 minRms true; volume #1005 step 1; split_diff_map #1005 #1006; close #1000,1001,1002,1003,1004,1005; view tmp
alias usage local_diff_map synopsis "Build a local difference map between two fitted maps around a model" $1 "model-spec" $2 "first map" $3 "second map"
```

### `local_diff_map_sphere`

Builds a local difference map from two maps around the center of rotation.

Usage: `local_diff_map_sphere <first-map> <second-map> <radius>`

Depends on: `split_diff_map`

```chimerax
alias local_diff_map_sphere ~sel; close #10000,#1001,#1002,#1003,#1004,#1005; marker #10000 position cofr; sel #10000; vol zone $1 near sel range $3 newMap true modelId 1000 minimalBounds true; vol zone $2 near sel range $3 newMap true modelId 1001 minimalBounds true; volume scale #1000 sd 0.1 modelId 1002; volume scale #1001 sd 0.1 modelId 1003; fitmap #1003 inMap #1002; volume resample #1003 onGrid #1002 modelId 1004; volume subtract #1002 #1004 modelId 1005 minRms true; volume #1005 step 1; split_diff_map #1005 #1006; close #10000,#1000,#1001,#1002,#1003,#1004,#1005
alias usage local_diff_map_sphere synopsis "Build a local difference map between two maps around the center of rotation" $1 "first map" $2 "second map" $3 "radius"
```

### `binary_mask`

Creates a binary mask map at a chosen contour level.

Usage: `binary_mask <map-spec> <contour-level>`

```chimerax
alias binary_mask volume threshold $1 minimum $2 set 0 maximum $2 setMaximum 1 modelId 2000; volume #2000 level 0.5
alias usage binary_mask synopsis "Create a binary mask map at a chosen contour level" $1 "map-spec" $2 "contour level"
```

### `soft_mask`

Creates a soft-edged mask map from a chosen contour level.

Usage: `soft_mask <map-spec> <contour-level>`

```chimerax
alias soft_mask volume threshold $1 minimum $2 set 0 maximum $2 setMaximum 1 modelId 2000; volume falloff #2000 iterations 20 modelId 2001; volume #2001 level 0.5; close #2000
alias usage soft_mask synopsis "Create a soft-edged mask map at a chosen contour level" $1 "map-spec" $2 "contour level"
```

### `delh`

Deletes all hydrogen atoms.

Usage: `delh`

```chimerax
alias delh delete @H*
alias usage delh synopsis "Delete all hydrogen atoms"
```

### `delh_sel`

Deletes hydrogen atoms from the current selection.

Usage: `delh_sel`

```chimerax
alias delh_sel delete sel & @H*
alias usage delh_sel synopsis "Delete hydrogen atoms from the current selection"
```

### `carve`

Carves displayed map surfaces to the current selection by the given distance.

Usage: `carve <distance>`

```chimerax
alias carve surface zone ~##num_residues nearAtoms sel distance $1
alias usage carve synopsis "Carve displayed map surfaces to the current selection" $1 "distance"
```

### `uncarve`

Removes surface carving from displayed map surfaces.

Usage: `uncarve`

```chimerax
alias uncarve surface unzone ~##num_residues
alias usage uncarve synopsis "Remove surface carving from displayed map surfaces"
```

## Notes

- The Python file is meant to be sourced inside ChimeraX with `runscript`; it is not a standalone Python program.
- `nextmodel` and `prevmodel` use the Model Panel tool and therefore depend on ChimeraX's current Model Panel implementation.
- If you move the repository, update the `runscript` path in [chimerax_trimmings.txt](chimerax_trimmings.txt).
