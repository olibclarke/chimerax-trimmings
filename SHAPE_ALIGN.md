# Shape Alignment

`shapealign` performs rigid shape-based alignment of atomic models and density
maps. Use it for overall superposition without sequence alignment, for aligning
large complexes, or for refining an existing alignment around a local region.
It supports model/model, model/map, map/model, and map/map fitting.

It is not a sequence-alignment tool or a flexible fit. Use Matchmaker when you
need sequence-based residue correspondence.

## Load the Command

The command is included in the main Trimmings script:

```chimerax
runscript /path/to/chimerax_trimmings.py
```

It is available automatically if that script is already in your ChimeraX Startup
preferences. No separate alignment script is required.

## Basic Usage

```chimerax
shapealign #2 to #1
```

The first specification is the moving input; `to` identifies the fixed reference.
Either ID can refer to an atomic model or a map. The reference stays in place.

To align several inputs independently:

```chimerax
shapealign #2-5 to #1
```

Moving inputs can mix atomic models and maps. They do not move together as a
single rigid assembly.

For atomic models, restrict the fitting atoms to chains, residue ranges, or a
selection when needed:

```chimerax
shapealign #2/A to #1/B
shapealign #2/A:100-220 to #1/A:100-220
shapealign sel to #1
```

The selected atoms determine the fit, but the entire moving model is repositioned.
Hidden atoms in the specification are included. Hydrogens and water are excluded.
The reference is skipped automatically if it is included in the moving selection.

## Prediction Confidence

`confidence auto` (the default) weights models identified by ChimeraX's
AlphaFold database metadata. Local files without that metadata are not guessed
from their name, colors, or B-factor range. For local predictions with pLDDT
stored in the B-factor field, enable it explicitly:

```chimerax
shapealign #2 to #1 confidence plddt
shapealign #2 to #1 confidence_moving plddt
shapealign #2 to #1 confidence_reference plddt
shapealign #2 to #1 confidence off
```

`confidence` sets the policy for both atomic inputs; `confidence_moving` and
`confidence_reference` override it independently. Each accepts `auto`, `off`, or
`plddt`. Use the side-specific options when only one input is a local prediction.
Do not force `plddt` on experimental B-factors. Maps are never confidence-weighted.
The Log reports which atomic inputs were weighted.

Weights rise smoothly from 0 at pLDDT 50 to 0.5 at 70 and 1 at 90; scores below
50 are ignored and scores above 90 have full weight. This affects centering,
orientation seeds, automatic fitting scales, and simulated densities throughout
search and refinement. It works with either objective, local refinement,
model/map fits, and parallel atomic batches. The entire moving model still moves
rigidly, including ignored atoms; no atoms are deleted or hidden.

At least four non-collinear atoms with positive weight must remain in each
atomic fitting region. Optional local correlation scores use the weighted
densities and are assigned only to retained fitting atoms. pLDDT does not measure
relative domain placement uncertainty, so inspect multi-domain fits carefully.

## Choose the Objective

The default `objective density` compares density overlap. For atom-only fits, it
can automatically use partial fitting when the structures have different sizes
or the initial fit is poor.

```chimerax
shapealign #2 to #1 objective density
```

The experimental `objective gradient` emphasizes matching shape features rather
than bulk density overlap. Try it when density fitting gives an implausible
orientation, particularly with different compositions. It always uses a
whole-density objective, without automatic partial or per-chain fitting.

```chimerax
shapealign #2 to #1 objective gradient
```

Gradient scoring can take longer and does not guarantee a better fit. Always
inspect the result. Scores from the two objectives are not directly comparable.

## Fitting Maps

Unless a CoFR sphere is requested, maps use their full density, independently of the displayed contour, surface
zone, displayed subregion, or display step. Positive and negative values are
retained. There is no automatic background subtraction, thresholding, or sign
inversion.

When a map is involved, both objectives use whole-density fitting. A broader
search may be tried for a poor fit, but it keeps the same objective.
`partial true` is available only for atom-only fits.

Only the moving input's rigid transform changes. Map voxel values, grid origin,
voxel spacing, contours, styling, and visibility remain unchanged. The command
does not create replacement volumes or use `volume zone`.

Nonzero origins, anisotropic voxel sizes, and rotated or skewed grids are
supported. Input units and voxel spacing must be correct. A moving map's surfaces
move with the map; a separately opened map is not moved automatically when its
associated atomic model moves.

Large nonzero backgrounds or noise can affect a full-density fit. Use an
appropriately prepared map or restrict fitting with a local sphere.

## Regional Global Search

For a small model in a large map, position the ChimeraX center of rotation at
the expected binding site and search the surrounding reference density:

```chimerax
shapealign #2 to #1 radius_search
shapealign #2 to #1 radius_search 60
shapealign #2 to #1 radius_search 60 objective gradient
```

With no value, `radius_search` uses 1.5 times the maximum distance from the moving
fitting atoms' centroid to any fitting atom, giving 50% padding around an enclosing
sphere. This is independent of the model's position or orientation and respects
atomic selections. For moving maps, it uses the farthest corner of the full
nonzero density bounds from the density center. In a batch, the automatic radius
is calculated separately for each moving input. The chosen radius is logged.

`to` must specify a map. `radius_search` is in Angstroms and selects reference
voxel centers inside the starting CoFR sphere. The entire moving selection is
used, even if it starts outside that sphere. A global orientation and translation
search is always performed; `fallback_correlation 0` does not disable this
explicitly requested search. `extra_starts` can extend the orientation search.

The radius selects the reference density used for scoring, not a maximum movement
distance or a requirement that every moving atom stay inside the sphere. Choose
a sphere large enough to contain the expected fit. Only temporary numerical
density is cropped: the input map, its surfaces, and display remain unchanged.

## Local Refinement

Use `radius_refine` after the inputs are already approximately aligned. Position
the ChimeraX center of rotation at the region of interest, then run:

```chimerax
shapealign #2 to #1 radius_refine 30
shapealign #2 to #1 radius_refine 30 objective gradient
```

The radius is in Angstroms. Only atom or voxel centers inside the starting sphere
contribute. The sphere is fixed for the command; each whole moving input still
moves rigidly. There is no global orientation search in this mode.

Atomic specifications can further restrict the sphere. For maps, the input map
is not cropped or zoned. `radius_refine` cannot be combined with `partial true`
or nonzero `extra_starts`. `radius_search` and `radius_refine` are mutually
exclusive. The old `radius` spelling remains an alias for `radius_refine`;
specify only one spelling.

## Resolution and Search

Resolution is chosen automatically from the physical size of the fitting
regions. To set the final fitting scale explicitly:

```chimerax
shapealign #2 to #1 resolution 6
```

Values are in Angstroms and must be at least 2. Larger values give smoother,
coarser fits and can reduce memory use. For maps, this scale adds Gaussian blur
to the existing density; it does not estimate or replace the map's intrinsic
resolution. The original map is unchanged.

For a more extensive global orientation search, increase `extra_starts`:

```chimerax
shapealign #2 to #1 objective gradient extra_starts 64
```

This costs more time. `fallback_correlation` controls when an additional search
is attempted; its default depends on the objective. Set it to 0 to disable
quality-triggered extra searches:

```chimerax
shapealign #2 to #1 objective gradient fallback_correlation 0
```

Batch fits reuse reference densities, gradient fields, and FFTs when their
settings match. The temporary cache is limited to 256 MiB per fitting process
and cleared after the command; fitting scratch space is additional. The benefit
depends on how much reference preparation the models share.

For atom-only density fitting, `partial true` explicitly requests a partial fit;
`partial false` disables partial fitting. Omitting it allows automatic selection.

## Parallel Batches

The reference is automatically skipped if it falls inside the moving range,
so broad specifications such as `shapealign #2-1000 to #13` are supported.
Selecting only the reference does nothing and adds no Undo action.

Global atomic-model batches with at least four moving models automatically use
the process-available logical CPU count minus one, capped by the number of
moving models (minimum one worker). CPU affinity limits are respected where
supported. To select a process count or force serial fitting:

```chimerax
shapealign #2-100 to #1 workers 4
shapealign #2-100 to #1 workers 1
```

`workers 0` selects automatic mode. Explicit counts are limited by the
available CPUs and number of moving models, without an additional fixed cap.
Local refinement stays serial in
automatic mode, but can use explicit workers. Map and mixed-input fits remain
serial; use `workers 0` or `workers 1` for those.

The command waits for the batch and shows progress. Worker fits use the same
settings as serial fits. Results are applied in model order after all fits
succeed, with one Undo action; a failed worker leaves the models unchanged.
If a fitting model changes during computation, the batch is not applied.
More workers consume more memory, including a separate preprocessing cache
and fitting scratch space for each process.

## Optional Fit Coloring

For a moving atomic model, optionally color the fitted atoms and cartoons by
local density correlation:

```chimerax
shapealign #2 to #1 color_fit true
shapealign #2 to #1 objective gradient color_fit true
```

The reference can be a model or a map. Colors run from red at 0 or lower, through
white at 0.5, to blue at 1. Gray indicates an undefined correlation.

Use `correlation true` to save values without recoloring. Values are stored as
the atom attribute `shape_correlation`; residues receive the mean over their
defined fitted atoms.

This is a local density-agreement measure, not RMSD or a confidence estimate.
It is separate from the gradient fitting score. Both options default to false,
require `move true`, and are unavailable for moving maps. They add calculation
time.

## Options

Type `shapealign` or `shapealign help` to open the Log with a short, formatted
quick reference and a status-bar prompt. Click **All options** or run
`usage shapealign` for the complete syntax.

| Option | Meaning | Default |
| --- | --- | --- |
| Moving specification | Atomic models, atom subsets, maps, or a mixture | Required |
| `to` | One fixed atomic model/atom subset or map | Required |
| `objective` | `density` or experimental `gradient` | `density` |
| `resolution` | Final fitting scale in Angstroms, at least 2 | Automatic |
| `max_points` | Maximum source density samples per stage, at least 1000 | `12000` |
| `extra_starts` | Additional global orientation starts | `0` |
| `workers` | `0` uses available CPUs minus one, `1` serial, or an explicit process count for atomic-model batches | `0` |
| `confidence` | `auto`, `off`, or `plddt` weighting for atomic inputs | `auto` |
| `confidence_moving` | Override confidence weighting for moving models | Inherit `confidence` |
| `confidence_reference` | Override confidence weighting for the reference | Inherit `confidence` |
| `move` | Apply the fit; false only calculates and reports it | `true` |
| `radius_search` | Global search against a reference-map CoFR sphere; optional radius in Angstroms, or omit the value for 1.5 times the moving input's enclosing radius | Off |
| `radius_refine` | Local refinement using both inputs' CoFR spheres, in Angstroms | None |
| `radius` | Compatibility alias for `radius_refine` | None |
| `partial` | `true` requests atom-only partial density fitting; `false` disables it | Automatic for atom-only density fitting |
| `fallback_correlation` | Cutoff for an additional search; 0 disables quality-triggered searches | `0.7` density; `0.5` gradient |
| `correlation` | Save local density CC on moving atoms/residues | `false` |
| `color_fit` | Save local density CC and color moving atoms/cartoons | `false` |
| `correlation_resolution` | Density scale for optional local CC, in Angstroms | `8` |
| `correlation_window` | Neighborhood width for optional local CC, in Angstroms | `12` |

## Undo and Practical Tips

Applied fits support ChimeraX Undo/Redo. A multi-input command is one undo action;
if a fit fails, earlier movements in that batch are rolled back. With
`move false`, positions remain unchanged and no undo action is added.

Symmetric complexes can have several equivalent placements. Different
conformations or compositions may not superimpose exactly with one rigid
transform. Restrict atomic selections or use local refinement to focus on the
region you need.

If memory limits are reached, use a coarser `resolution`, a smaller atomic
selection, a regional `radius_search`, or a local `radius_refine`. The command
cannot fit instanced model copies or treat both a parent and its descendant as
independent moving inputs.

No finite search guarantees the global optimum. Inspect the aligned structures
and maps rather than treating the reported score as proof of a correct fit.
