# Ligand Search

Screen CCD components against an isolated positive-density blob and browse a
running shortlist. This experimental tool suggests ligand hypotheses; density
fit alone does not establish chemical identity.

Keep `ligand_search.py` beside `chimerax_trimmings.py`, then load the main script
as usual. No additional Python packages are required beyond ChimeraX.

## Quick Start

```chimerax
ligandsearch
```

1. Choose the density map and enter the approximate local resolution in Angstroms.
2. Set **Blob contour** high enough to separate the ligand blob from nearby density.
3. Click **Pick blob**, then click its map surface. Alternatively center the view
   on the blob and use **Use CoFR**.
4. Click **Start search**. The CCD library and PDB occurrence counts download on
   first use, with progress shown in the tool. Subsequent searches use the cache.
5. Browse hits with **Previous / Next** or the table. **Focus** centers the preview;
   **Keep hit** creates an atomic model with Undo support. **Stop** retains the
   completed shortlist. **Export CSV** saves the current top hits.

To start directly around the current center of rotation:

```chimerax
ligandsearch #2 resolution 2.6 level 0.03 top 20
ligandsearch #2 resolution 2.6 level 0.03 fast true
ligandsearch #2 resolution 2.6 fast true torsions true conformers 8
ligandsearch stop
```

The map, contour, colors and existing atomic models are not changed. One temporary
ligand preview is reused and follows the highest-ranked hit as results arrive.
Selecting a table row or using **Previous / Next** pauses live following so you
can inspect that candidate, even if it leaves the shortlist. Check **Follow top
hit** to resume. Each new search enables live following again.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| map specifier | GUI selection | Density map; supplying it starts a search around CoFR. |
| `resolution` | `2.6` | Approximate local map resolution, in Angstroms. |
| `level` | Positive displayed contour | Contour used to identify the connected blob; does not change the displayed contour. |
| `radius` | `14` | Map region radius around the picked point or CoFR, in Angstroms. Increase it if the blob reaches the region boundary. |
| `top` | `20` | Maximum retained CCD hits, from 1 to 100; one best pose per component. |
| `workers` | `0` | CPU workers; zero chooses automatically, leaving one available core free where possible. |
| `fast` | `false` | Opt-in fast search: 8 extra orientation starts and up to 1,500 density samples, instead of 24 and 3,000. Same candidate set and final scorer, but poses and rankings may change. |
| `frequency` | `true` | Combine volume match with PDB occurrence frequency for search ordering. `false` uses volume alone. |
| `observed` | `false` | Also try the CCD example observed-coordinate conformer. Slower; this is not a comprehensive conformer search. |
| `torsions` | `false` | Generate a sparse, deterministic set of additional conformers around flexible single bonds. Slower; not exhaustive or energy-minimized. |
| `conformers` | `8` | Maximum total conformers per CCD entry when `torsions true`, including selected CCD originals. Range 2-32; rigid ligands may yield fewer. |
| `codes` | All indexed entries | Optional comma-separated CCD codes, for example `"ATP,ADP,AMP"`. |
| `library` | ChimeraX cache | Folder containing the local CCD index and occurrence counts. Quote paths containing spaces. |

```chimerax
ligandsearch #2 resolution 2.8 level 0.03 radius 14 top 10 workers 4 fast true frequency true observed false torsions true conformers 8 codes "ATP,ADP,AMP" library "~/Desktop/ligand library"
```

**Fast search** is also available as a checkbox in the browser. It does not
reduce the CCD library or change the map resolution, conformer choices or final
density-scoring formula. Standard search remains the default; omit `fast` or
use `fast false` for the original settings. CSV exports record the search mode
in `fast_search` for each hit.

## Sparse Conformers

Enable **Sparse torsion sampling** in the browser, or use `torsions true`.
The browser enables **Max conformers per ligand** and the observed-coordinate
conformer checkbox only while sparse torsion sampling is selected and no job
is running.

The default cap is eight total conformers per ligand, including the original
ideal conformer and, with `observed true`, the selected CCD observed conformer.
Originals are retained; the library is not modified and needs no rebuild.

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

The [public CCD archive](https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz)
is indexed locally. Occurrence counts follow RCSB's
[bulk ligand-composition approach](https://www.rcsb.org/docs/general-help/ligandexpo-shutdown).
**Expand/Update CCD library** refreshes both. A failed update keeps the previous
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
