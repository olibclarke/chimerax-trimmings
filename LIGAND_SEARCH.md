# Ligand Search

Screen CCD components and ChEBI metabolites against an isolated positive-density blob and browse a
running shortlist. This experimental tool suggests ligand hypotheses; density
fit alone does not establish chemical identity.

Keep `ligand_search.py` and `ligand_chebi.py` beside `chimerax_trimmings.py`, then load the main script
as usual. No additional Python packages are required beyond ChimeraX.

## Quick Start

```chimerax
ligandsearch
```

1. Choose the density map and enter the approximate local resolution in Angstroms.
2. Display the desired positive map contour.
3. Click **Pick blob**, then click its map surface. Picking highlights the connected
   surface in green, using ChimeraX's native blob-coloring API. Search is enabled
   only after a pick; CoFR is not used. Contour and working region are handled
   automatically. If needed to separate the blob, the contour is raised with a
   message and the displayed contour is updated too.
4. Choose **All compounds** or select chemical groups, then click **Start search**.
   The first library search downloads one combined CCD + ChEBI archive;
   progress is shown in the tool. Subsequent searches use the cache.
5. Browse hits with **Previous / Next** or the table. **Focus** centers the preview;
   **Keep hit** creates an atomic model with Undo support. **Stop** retains the
   completed shortlist. **Export CSV** saves the current top hits.

To choose a map and activate blob picking:

```chimerax
ligandsearch #2 resolution 2.6 top 20
ligandsearch #2 resolution 2.6 fast true
ligandsearch #2 resolution 2.6 fast true torsions true conformers 8
ligandsearch #2 resolution 2.6 fast true refine true
ligandsearch stop
```

Map data and existing atomic models are not changed. Only the picked contour is
adjusted if needed; other contours, transparency and display style are preserved.
The temporary green highlight is cleared on a new pick, map change
or closing the tool, restoring the previous surface colors. An automatically
raised contour is retained. One temporary ligand preview is reused and follows
the highest-ranked hit as results arrive.
Selecting a table row or using **Previous / Next** pauses live following so you
can inspect that candidate, even if it leaves the shortlist. Check **Follow top
hit** to resume. Each new search enables live following again.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| map specifier | GUI selection | Density map; supplying it activates blob picking, not an automatic search. |
| `resolution` | `2.6` | Approximate local map resolution, in Angstroms. |
| `top` | `20` | Maximum retained CCD hits, from 1 to 100; one best pose per component. |
| `workers` | `0` | CPU workers; zero chooses automatically, leaving one available core free where possible. |
| `fast` | `false` | Opt-in fast search: 8 extra orientation starts and up to 1,500 density samples, instead of 24 and 3,000. Same candidate set and final scorer, but poses and rankings may change. |
| `frequency` | `true` | Combine volume match with PDB occurrence frequency for search ordering. `false` uses volume alone. |
| `observed` | `true` | Try the CCD example observed-coordinate conformer when available and distinct, even without torsion sampling. Use `observed false` to prefer ideal coordinates only. |
| `torsions` | `false` | Generate a sparse, deterministic set of additional conformers around flexible single bonds. Slower; not exhaustive or energy-minimized. |
| `conformers` | `8` | Maximum total conformers per ligand when `torsions true`, including its starting conformers. Range 2-32 for library searches or 2-256 for a single ligand. |
| `refine` | `false` | After screening, run a broader sparse torsion search and local rigid fitting on the retained top hits. |
| `refine_conformers` | `64` | Maximum conformers per shortlisted hit in the final pass, including the starting conformation. Range 2-256 for library searches or 2-512 for a single ligand. |
| `codes` | All indexed entries | Optional comma-separated CCD codes, for example `"ATP,ADP,AMP"`. |
| `ccd_code` | none | Fit one CCD entry, using the local index or downloading only that entry. Omit `codes`, `groups`, `sources` and `smiles`. |
| `smiles` | none | Fit one custom ligand using ChimeraX's built-in online SMILES converter, without a library search. Quote the string; omit `codes`, `groups` and `sources`. |
| `groups` | `all` | Comma-separated chemical-family or role keys. Intersects with `codes` when both are supplied. |
| `group_match` | `any` | `any` accepts any selected group (OR); `all` requires every selected group (AND). |
| `sources` | `ccd_chebi` | Search both libraries, or use `ccd` to search only CCD. |
| `library` | ChimeraX download cache | Override the local ligand-library folder. Quote paths containing spaces. |

```chimerax
ligandsearch #2 resolution 2.8 top 10 workers 4 fast true frequency true torsions true conformers 8 codes "ATP,ADP,AMP" library "~/Desktop/ligand library"
```

**Fast screening** is also available as a checkbox in the browser. It does not
reduce the CCD library or change the map resolution, conformer choices or final
density-scoring formula. Standard search remains the default; omit `fast` or
use `fast false` for the original settings. CSV exports record the search mode
in `fast_search` for each hit.

## Fit One Ligand

Choose **Single ligand (CCD code)** in **Compounds**, enter a code such as
`ATP`, pick the blob and start the search. The tool uses the local CCD index
if the entry is present, otherwise it downloads just that entry from RCSB.
No full library download is needed. Both torsion passes are available.

```chimerax
ligandsearch #2 ccd_code ATP resolution 2.6 torsions true conformers 64 refine true
```

Choose **Single ligand (SMILES)** in **Compounds**, enter a stereospecific
SMILES, pick the blob and start the search. Or pass the SMILES on the command line:

```chimerax
ligandsearch #2 smiles "CCOc1ccccc1" resolution 2.6 torsions true refine true
```

The built-in ChimeraX converter generates the initial 3D coordinates online.
Only the SMILES string goes to the conversion service, never the density map or
model; don't use confidential structures. Enter one molecule rather than a salt
or mixture. No CCD or ChEBI download is needed. Both torsion passes work as
usual, and the CSV records the input string.

Both single-ligand modes allow up to 256 screening conformers and 512 in the
final pass. The defaults are still 8 and 64; higher caps can take longer.

## Chemical Groups

Choose groups from the browser's **Chemical groups** menu. The adjacent selector
chooses **Any (OR)** or **All (AND)**. For example,
select **Organic nitrogen compounds** and **Metabolites**, then choose AND to
require both. The same choice is available in commands:

```chimerax
ligandsearch #2 groups nucleotides
ligandsearch #2 groups "steroids,fatty_acids"
ligandsearch #2 groups metabolites
ligandsearch #2 groups "organic_nitrogen,metabolites" group_match all
ligandsearch #2 groups all
```

Common families: `nucleotides`, `nucleosides`, `nucleobases`, `amino_acids`, `peptides`,
`lipids`, `fatty_acids`, `phospholipids`, `glycerolipids`, `sphingolipids`,
`steroids`, `sterols`, `isoprenoids`, `carbohydrates`, `glycosides`, `porphyrins`.
Nucleotides, amino acids and carbohydrates include annotated derivatives.
Additional presets include `benzenoids`, `heterocycles`, `organic_acids`,
`organic_oxygen`, `organic_nitrogen`, `organohalogens`, `organosulfur`,
`organophosphorus`, `hydrocarbons`, `alkaloids`, `phenylpropanoids`,
`organometallics` and `inorganic`.

**More structural classes...** opens a searchable list of detailed ChemOnt
families with counts, such as pyrimidines or carbonic acid monoesters. Select
any combination. Each row also shows its command key, for example
`chemont:0001655`; these keys are accepted by `groups` too.

Roles: `metabolites`, `drugs`, `inhibitors`, `agonists`, `antagonists`, `toxins`,
`cofactors`, `vitamins`, `detergents`, `buffers`. A role describes a reported use
or biological activity, not physiological presence in a particular sample.

Groups overlap: cholesterol can be a lipid, steroid and sterol. OR searches their
union; AND searches their intersection, without duplicating entries. Counts show membership
in the indexed library, not counts of distinct hits across several groups.
Group selection changes eligibility, not density scores or search ordering.
Changing groups also filters the retained shortlist, density preview and CSV
export immediately, without discarding the other retained hits. The browser
shows which groups the underlying search covered. Start a new search to screen
the newly selected groups fully: filtering an existing top-hit buffer cannot
recover compounds that were never tested or retained. Groups still combine by
the chosen OR/AND rule, including any selected detailed structural classes.

Map picking, resolution, compound source and chemical filters stay above
**Top Hits**. **Search options...** opens a separate dialog for torsion sampling,
screening speed, CPU workers, compound IDs and search notes. **Library** contains
cache, install and update actions. Hit details scroll independently, so long
chemical names do not squeeze the results table down to its header.
Compound names appear beside their IDs, ahead of the numeric metrics. Hover
over truncated names or column headings for full text and explanations. The
progress bar tracks screening and the final torsion pass. Hit actions are enabled
only when applicable; CSV export still permits stale results, marked as such.
Entry fields ignore wheel/trackpad scrolling, including while focused. Click the
arrows, type a value, or explicitly open a dropdown to change it.

## Library Download

On the first library search, Ligand Search downloads and verifies the
[combined library](https://github.com/olibclarke/chimerax-trimmings/releases/tag/ligands-chebi257)
(61 MB download, 275 MB installed). It includes 51,155 CCD components and
14,888 additional ChEBI metabolites, plus chemical groups and PDB occurrence
counts. Searches work offline after installation. **CCD + ChEBI metabolites**
is the default. Choose **CCD** to restrict the search to the CCD.

The library is stored in `LigandSearch` inside ChimeraX's configured download
cache (typically `~/Downloads/ChimeraX/LigandSearch` on macOS). Existing caches
are reused and are never replaced by the prebuilt archive. To install the
combined library alongside an older cache, set a new **Cache folder** under
**Library → Library status and cache...**, or pass `library` when opening the tool.

**Library → Install library archive...** installs a completed combined-library
archive into an empty cache folder. Use this if automatic download is blocked
or if you need an offline install. If a download is interrupted, start the
search again to retry. See [data sources and attribution](DATA_NOTICE.md#ligand-search-archive).

Chemical families are computed directly from the CCD structure using the
[GNPS2 ClassyFire2 service](https://apps.gnps2.org/classyfire), independently of
whether the compound has a ChEBI record. The primary ChemOnt lineage supplies
family tags; speculative alternative matches are not used for filtering.
These are automated structural classifications, not expert-validated assignments.
The broad **Organic nitrogen compounds** preset instead requires an explicit
carbon–nitrogen bond, so it also includes nitrogen-containing amino acids,
nucleotides and other families outside the primary organic-nitrogen taxonomy
branch. Existing caches are corrected locally, without downloading coordinates
again. Detailed ChemOnt classes remain available separately.

[ChEBI](https://www.ebi.ac.uk/chebi/downloads) supplements the structural families
and provides curated role annotations via [RCSB cross-references](https://data.rcsb.org/)
and structure matching. Connectivity-level matches contribute only tags shared
by all matching ChEBI structures. CCD monomer types also identify amino-acid,
nucleotide and saccharide components. Compound names are not used to guess tags.
Hit tooltips and CSV exports retain the classification source and identifiers.

The browser reports actual class-level classification coverage, separately from
biological-role coverage. Unsupported or invalid structure strings are not assigned
invented categories; the classifier retries an existing alternative CCD SMILES
where available. `unclassified` selects entries with no tags from either source.
**All compounds** still includes every indexed component.

**Update chemical groups** refreshes annotations without rebuilding coordinates.
The annotations are also refreshed by **Expand/Update CCD library** and fetched
automatically if a filtered search needs them. Cached searches work offline;
an interrupted update preserves the previous complete group cache and resumes
completed structural classifications on retry. Results are cached by structure
and classifier rule version. Only public CCD structure strings are sent to the
classification service, never the experimental map or fitted model. No additional
Python packages are needed in ChimeraX, and fitting/searching uses the local cache.

## ChEBI Metabolites

With **CCD + ChEBI metabolites** selected (the default), try:

```chimerax
ligandsearch #2 resolution 2.6 groups metabolites
ligandsearch #2 groups "organic_nitrogen,metabolites" group_match all
```

The expansion uses defined molecules with a ChEBI metabolite-role
annotation, including inherited annotations. Exact identity matches reuse CCD
coordinates and appear only once. Additional compounds use PubChem3D conformers;
2D records, changed stereochemistry/charge states and disconnected salts are not
silently substituted. Molecules without usable 3D coordinates are reported as
unavailable, not included as unvalidated structures. No RDKit or Open Babel is
needed.

No separate ChEBI download is needed with the combined archive. With an older
CCD-only cache, choose **CCD** for an immediate search, install the combined
archive in a new cache, or use **Library → Expand/Update ChEBI metabolites** to
build the expansion in place. **Stop** preserves completed entries so the build
can resume. Searches can use a partial cache, with a warning.

Both sources share the volume/PDB-frequency queue and density scorer; candidates
are not prioritized merely because their coordinates came from CCD. Group
filters apply to both sources. Counts in the group menu describe the selected
library, without counting reused entries twice. Expanded entries carry ChEBI ontology tags and the
bond-based organic-nitrogen tag, not automatically inferred ChemOnt classes.
CSV exports include `library_source`, `chebi_ids`, `pubchem_cid`, `inchikey` and
`formal_charge` where available. The compound-ID filter also accepts ChEBI IDs,
for example `codes "CHEBI:15422"`. Public identity queries send no experimental
map or fitted coordinates to ChEBI, RCSB or PubChem.

## Sparse Conformers

Enable **Sample torsions during screening** in the browser, or use `torsions true`.
The browser enables **Search conformers** only while sparse torsion sampling is
selected and no job is running. The label shows the cap for the selected mode.

CCD screening tries the ideal conformer and a distinct observed conformer when
available, with or without torsion sampling. The default sampling cap is eight
total conformers per ligand, including these originals. The library is not
modified and needs no rebuild.

Additional conformers apply a few reproducible single-bond rotations, including
combinations of up to three torsions. Ring bonds, terminal heavy-atom bonds,
metal coordination and common restricted/conjugated linkages are excluded
conservatively. Ring shapes, bond lengths, bond angles and tetrahedral handedness
are preserved. Duplicate shapes and newly introduced gross internal heavy-atom
overlaps are rejected. Ligands without eligible bonds simply use their existing
conformers. This is not an exhaustive torsion grid or a force-field calculation.

Each accepted conformer gets the same rigid search and final density scorer.
The browser retains one best result per CCD entry and reports how many conformers
were tested. CSV exports include the sampling flag, tested/skipped counts and
the winning conformer's bond rotations in degrees relative to its CCD seed.
No receptor clash tests are added. Sampling can multiply fitting time for flexible
ligands, and more conformers can also increase the opportunity to fit noise.

### Final Top-Hit Search

Enable **Final torsion search on top hits** and set **Final conformers**,
or use:

```chimerax
ligandsearch #2 fast true refine true refine_conformers 64 top 20
```

After screening finishes, the retained hits (normally 20) get a larger sparse
search starting from their fitted conformations. This allows large torsion
changes, including combinations of up to three rotations, not just small-angle
adjustments. Each variant is locally fitted from the existing pose; the global
search is not repeated. Chemical exclusions and internal-overlap checks are the
same as for initial sparse sampling. Rigid ligands get one local rigid fit.

The original pose remains eligible and is retained unless a variant improves
the same density score. Hits are reordered as the final pass completes. Details
and CSV report before/after scores, tested/skipped conformers and winning final
torsion changes. **Stop** keeps completed improvements and the remaining original
hits. This is independent of `torsions`: you can screen rigid conformers quickly,
then sample torsions only on the shortlist. It does not revisit compounds outside
that shortlist or perform force-field energy minimization.

## Ordering And Scores

The queue favors components whose approximate density volume matches the blob.
With `frequency true`, a logarithmic PDB-frequency boost, capped at threefold,
also brings commonly observed components forward. Frequency counts distinct PDB
entries, not ligand copies within an entry. Branched components and nonstandard
polymer monomers are included; standard residues within protein/RNA chains do
not inflate their counts. A component with no recorded occurrences remains
eligible.

**Frequency and volume determine search order only.** Hits are ranked by density
correlation, with a soft completeness factor for atom support and blob coverage.
The browser shows the score, raw density CC, coverage, supported atom fraction,
volume ratio and PDB occurrence count separately. These are comparative measures,
not identification probabilities or binding affinities. Scores remain provisional
while candidates are still being tested.

The CCD index is derived from the [public CCD archive](https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz).
Occurrence counts follow RCSB's
[bulk ligand-composition approach](https://www.rcsb.org/docs/general-help/ligandexpo-shutdown).
**Expand/Update CCD library** refreshes coordinates, occurrence counts and chemical
groups. A failed update keeps the previous
complete cache; missing occurrence data falls back to volume-only ordering with
a warning. Searching an already downloaded library works offline.

## Limitations

- Rigid CCD conformers can miss flexible ligands. Optional sparse torsion sampling
  covers only a small subset of shapes; ring conformations and restricted bonds
  are not sampled, and generated conformers are not energy-minimized.
- An observed CCD conformer can originate from the structure being investigated;
  it is not necessarily an independent validation control.
- Components need complete, noncollinear heavy-atom coordinates with at least four
  atoms. Ions and tiny fragments are not screened.
- Map regions must have orthogonal grid axes. Resample skew crystallographic grids
  in ChimeraX first. The supplied resolution cannot be finer than map sampling.
- No receptor clash filtering is performed. Inspect chemistry, map support,
  experimental context and possible receptor rebuilding before accepting a hit.
- A complete CCD search may take substantial time. Results populate progressively;
  stop when you have enough candidates to inspect.
