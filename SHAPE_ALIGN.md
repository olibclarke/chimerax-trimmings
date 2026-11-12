# Shape Alignment

`chimerax_trimmings.py` includes `shapealign`, a sequence-independent rigid
alignment helper for atomic structures. It is an alternative to Matchmaker when
the task is overall superposition, not sequence correspondence. Loading the main
Trimmings script registers this command along with the other helpers; there is no
separate alignment script or additional startup line to install.

If upgrading from the separate helper, remove any old `runscript` line for
`shape_align.py` from Startup preferences. `map_salami.py` remains separate.

## Usage

```chimerax
runscript /path/to/chimerax_trimmings.py
open 6qzp
open 9s3d
shapealign #1 to #2
```

The moving selection can contain one or more atomic models; `to` must identify
exactly one fixed reference model, excluded from the moving selection. Models
move independently, not as a rigid assembly. A chain, residue range, or selection
can restrict the fitting data, for example:

```chimerax
shapealign #1/A to #2/B
shapealign #1 & nucleic to #2 & nucleic
shapealign sel to #2 resolution 6 max_points 20000 extra_starts 64 move false
shapealign #2-10 to #1
shapealign #2-10 to #1 radius 30 color_fit true
shapealign #1 to #2 objective gradient
```

Each whole source model moves, including any child models. Separate maps are not
moved along with it. Hydrogens and water do not contribute to fitting. Display
masks do not restrict fitting; use explicit atom specifications instead.

| Option | Meaning | Default |
| --- | --- | --- |
| `objective` | Existing density scoring or experimental whole-density gradient scoring | `density` |
| `resolution` | Override the final Gaussian resolution in Angstroms, at least 2 | Automatic when omitted |
| `max_points` | Maximum source density samples per stage, at least 1000 | `12000` |
| `extra_starts` | Extra broad-search orientations, nonnegative | `0` |
| `move` | Apply the calculated rigid model transform | `true` |
| `radius` | Refine only the existing pose using atoms inside this radius of the current CoFR, in Angstroms | None (global search) |
| `correlation` | Save local density CC attributes without recoloring | `false` |
| `color_fit` | Save local density CC attributes and color fitted atoms/cartoons | `false` |
| `correlation_resolution` | Gaussian density resolution for scoring, in Angstroms; independent of fitting resolution | `8` |
| `correlation_window` | Gaussian correlation-neighbourhood FWHM, in Angstroms | `12` |
| `partial` | Force partial-template search (`true`) or disable it (`false`) | Automatic when omitted |
| `fallback_correlation` | Quality cutoff for additional search; 0 disables it | `0.7` density; `0.5` gradient |

`move false` leaves all model transforms unchanged and does not add an undo
action. Applied fits are undoable and redoable. Atomic coordinates, reference
position, model visibility, camera, and existing maps are unchanged. Coloring
changes only with `color_fit true`.
No temporary models are opened, and no model IDs are consumed.

A batch processes models in ID order and forms one undo/redo action, including
any correlation attributes and colors. If a fit fails, previous changes in that
batch are undone automatically. Moving models must not contain one another.
Python callers receive `(transform, stats)` for a single model or a list of those
pairs for a batch; stats identify the moving and reference model IDs.


## Experimental Gradient Objective

```chimerax
shapealign #1 to #2 objective gradient
shapealign #1 to #2 objective gradient fallback_correlation 0
shapealign #1 to #2 objective gradient resolution 4
```

This opt-in mode uses a whole-structure density-gradient objective, not a partial
or per-chain fit. All fitted atoms contribute to simulated density, and both
complete densities contribute to normalization. One rigid transform moves each
requested source model. The existing `objective density` is unchanged.

The continuous gradient inner product can be evaluated by integration by parts:

```text
score(T) = <density_m,T, -Laplacian(density_r)> / sqrt(E_m * E_r)
E = <density, -Laplacian(density)>
```

The implementation uses an isotropic finite-difference Laplacian with zero
boundary values on the padded density grids. Negative values around density
peaks are important penalties and are not clipped away. A scalar transformed
target field allows ChimeraX's compiled sum-product fitter to optimize the same
objective that ranks candidates. This emphasizes matching density features over
broad envelope overlap; it does not encode residue identity or TM annotations.
Finite-grid interpolation and sampling make normalization approximate. The score
is not Pearson correlation, a probability, or comparable to density-mode scores.

The initial size-derived pyramid tries principal-axis poses with up to 48 distinct
poses retained between levels. If its final score is below `fallback_correlation`
(default 0.5 for this objective), a broader gradient FFT search generates up to
48 starting poses. With at least three pyramid levels, that search starts at
the second level, bypassing the broadest blur that can lose helix detail. The
remaining levels refine those poses and retain the original solution unless a
higher score is found. FFT normalization uses complete density energies, not
template-mask Pearson statistics. The smaller input may be used as the numerical
FFT template for efficiency; neither input is trimmed from scoring.

`fallback_correlation 0` disables the FFT backup. `partial true` is rejected
rather than silently changing the objective. No automatic partial or chain-seeded
fallback is used. `radius` still enables current-pose-only local refinement with
no global search. Optional correlation attributes/coloring remain the independent
local density Pearson metric described below, not the gradient objective.

The log identifies gradient scoring and any FFT backup. Stats include `objective`,
`gradient_search`, and, if searched, `gradient_search_resolution`,
`gradient_search_starts`, `gradient_search_orientations`,
`gradient_search_fft_voxels`, `gradient_initial_score`, and
`gradient_search_improved`. Memory limits still apply. This finite experimental
search is not guaranteed to find a global optimum.

## Partial Search and Backup

This section describes the existing `objective density` mode only.

```chimerax
shapealign #2 to #1
shapealign #2 to #1 partial true
shapealign #2 to #1 partial false
shapealign #2 to #1 fallback_correlation 0.65
```

With `partial` omitted, a smaller/larger fitted atom-count ratio below 0.75
selects partial search directly. Otherwise the fast whole-shape method runs
first; if its normalized density-overlap score is below `fallback_correlation`,
partial search is used as a backup, retaining the initial fit as another seed.
The cutoff is a tunable heuristic, not calibrated confidence. It is not the
atom-level local correlation coloring score. A partial score is normalized
differently; its additional chain-seed quality gate is capped at 0.4, as
described below, rather than re-running the whole search recursively.

`partial true` forces the partial search; `partial false` disables the partial
and chain-seeded searches. `fallback_correlation 0` disables quality-triggered
searches but retains size-triggered partial search. Local CoFR
refinement never starts a global search, and `partial true` with `radius` is
rejected. Logs identify the chosen search and any fallback trigger.

The numerical search always places the smaller template inside the larger
target. When the user-selected moving model is larger, the result is inverted
before applying it: the requested moving model still moves, not the reference.
No sequence alignment or atom correspondence is used in fitting.

Linear FFT correlations scan all coarse-grid translations for principal-axis
and reproducible broad rotations. Target neighbourhoods containing roughly a
template-sized atom count provide extra principal-axis orientation hints. These
maintain capture at finer small-object resolutions without deliberately blurring
the object into a featureless envelope. Every orientation still scans every
translation; the hints do not restrict fitting to those neighbourhoods.

Masked Pearson normalization uses the template density above 0.1% of peak as
support. Unmatched target density outside that support is not included in its
denominator. Multiple translation peaks and distinct orientations are retained,
then refined through the size-derived pyramid with ChimeraX's native
correlation-about-mean optimizer. The finite search is not an exhaustive global
guarantee. Its FFT workspace is capped at 4 million voxels; restrict the regions
or use a coarser explicit resolution if necessary.

## Local Refinement

After globally aligning a pair, place the center of rotation at the region of
interest using your normal CoFR workflow, then run:

```chimerax
shapealign #1 to #2 radius 30
shapealign #1 to #2 radius 30 resolution 4
shapealign #1/A to #2/B radius 20 resolution 4 move false
```

`radius` is a positive, finite distance in Angstroms. It switches the command to
local refinement rather than merely cropping the input to a new global search:

- Capture the existing CoFR in scene coordinates; do not set `cofr centerofview`
  or replace a manually chosen pivot.
- Independently select atom centers inside that same sphere in each input model,
  intersecting the supplied atom specifications. Hidden atoms still contribute;
  water and hydrogens do not. At least four non-collinear heavy atoms are needed
  in each cropped region.
- Freeze both sets of fitting atoms before optimization. The ROI is not moved,
  recentered, or resampled to add new atoms during the fit. Gaussian density
  tails extend beyond those atom centers, but outside-sphere atoms contribute
  nothing to the objective.
- Start from the existing relative model pose only. There is no principal-axis
  search, extra orientation search, or preliminary centroid alignment. Refine
  through the size-derived local range, no broader than twice its fine endpoint.
  Nonzero `extra_starts` is rejected in this mode.
- Apply the best local fit as one rigid transform of the entire moving model,
  including atoms outside the sphere. This does not move just the selected atoms
  or introduce local structural deformation. The reference remains fixed.

Ordinary parent-model transforms are included when locating the sphere. Local
fits support the same undo/redo and dry-run behavior as global fits. The log also
reports the radius and number of fitting atoms from each model.
Batch fits capture the sphere once, before any model moves, even when CoFR is
normally computed automatically from the scene.
Neither size mismatch nor a low score triggers global reinitialization here.

This mode assumes a good global placement. A small, featureless or symmetric
patch may not define a unique local pose, and hard sphere boundaries can include
different atoms in the two structures. Enlarge or reposition the sphere if
needed. The radius restricts fitting data, not the maximum permitted translation
or rotation; this is local optimization, not a hard motion restraint.

## Optional Fit Coloring

```chimerax
shapealign #2 to #1 color_fit true
shapealign #2 to #1 correlation true
shapealign #2-10 to #1 radius 30 resolution 4 color_fit true
shapealign #2 to #1 color_fit true correlation_resolution 6 correlation_window 10
```

`color_fit true` implies correlation scoring. `correlation true` computes the
same attributes but leaves colors unchanged. Both default to false; no density
correlation field, attribute assignment, or coloring is done on the normal path.
They require `move true`: hypothetical dry-run scores must not be assigned to a
model still in a different pose.

After fitting, two equal-atom-weight Gaussian densities are generated on the
same grid in the reference frame, using only the fitting atoms. In local mode,
these are the atom sets captured in the initial CoFR sphere. Full-model atoms
outside that sphere do not enter the score, even near its boundary.

For a Gaussian window W, the field is the weighted Pearson correlation:

```text
meanA = W(A); meanB = W(B)
CC = (W(A*B) - meanA*meanB)
     / sqrt((W(A*A) - meanA*meanA) * (W(B*B) - meanB*meanB))
```

Separable compiled SciPy filters compute the five moments once per grid, instead
of extracting a neighbourhood for every atom. The window is Gaussian with
`correlation_window` FWHM and a three-sigma cutoff. The shared grid uses spacing
`correlation_resolution/3` and enough padding for the density and window kernels.
Input density arrays are reused as scratch space to limit peak memory.

The field is interpolated at each fitted moving atom center. Defined values lie
between -1 and 1 and are stored as atom attribute `shape_correlation`. The same
attribute on each fitted residue is the arithmetic mean over its defined fitted
atoms, not over unfitted atoms or undefined samples. This is spatial density
agreement, not matched-residue RMSD or a confidence probability. Pearson
correlation ignores positive density scaling and uniform offsets; it does not
measure absolute mass agreement. Repeated motifs or symmetry-related subunits
can agree in shape without representing the same sequence position.

Near-zero variance in either density makes correlation undefined. The calculation
uses both an absolute variance floor (1e-10 after normalizing each map to unit
peak) and a relative floor (variance greater than 1e-5 times squared local mean)
to avoid unstable divisions in nearly flat regions. Interpolation requires valid
support from all contributing voxels (at least 0.999 valid weight). Undefined
samples have no attribute, rather than a fabricated zero score.

Colors use a fixed scale: CC <= 0 is red, 0.5 white, and 1 blue; undefined is gray.
Atoms use their own values and cartoons use residue means. Existing transparency
is retained. Only fitted atom colors and the corresponding residue cartoons are
changed; reference colors, surfaces, hydrogens, and water remain untouched.
Native `color byattribute` can reuse the saved atom (`a:shape_correlation`) or
residue (`r:shape_correlation`) values with other palettes.

Each scored fit replaces previous correlation attributes on that moving model,
clearing stale values outside a new selection/ROI while leaving their colors
unchanged. Scores are snapshots: rerun scoring after changing model positions,
coordinates, selections, or the reference. An unscored fit leaves older values
and colors untouched. Undo/redo restores positions, attributes (including their
absence), and any changed colors together.

Scoring is optional because it has a real cost. On 6QZP/9S3D, default scoring and
coloring added about 1.1 seconds to a 2.1-second numerical fit (about 3.4 seconds
total including preparation). 217,455 of 219,419 fitted heavy atoms had defined
scores in that trial. The log reports fitting and scoring/coloring times
separately; timings vary with grid size and model size. Increasing
`correlation_resolution` makes scoring coarser and usually faster.

The correlation grid is capped at 8 million voxels because several moment arrays
are needed simultaneously; `max_points` does not limit this full-grid work. If
the cap is exceeded, use coarser `correlation_resolution` or a smaller fitting
region. Correlation is computed before moving the model, so a scoring-grid error
does not leave a partially applied fit. No displayed map or model is created.

## Method

1. Read non-water heavy-atom coordinates in each model's local frame. This avoids
   scene-coordinate cache issues and prevents scene rotations from changing
   density voxelization. Each atom contributes equally, rather than using
   molmap's element-number weighting.
2. For global alignment, compute both principal-axis frames from the full selected
   atom clouds. Try
   all 24 proper signed axis permutations, never reflections. Also include
   identity/centering starts and the existing relative scene placement. Optional
   extra orientations use unit quaternions with a fixed random seed.
3. At each scale, bin nearby atoms into cells of width resolution/6 and retain
   each cell's centroid and atom count. This reduces Gaussian-deposition work
   without discarding the relative density of RNA and protein regions.
4. Generate two in-memory grids using ChimeraX's compiled Gaussian deposition.
   Gaussian sigma is resolution/(pi*sqrt(2)), matching molmap's convention; grid
   spacing is resolution/3, with a four-sigma border and Gaussian cutoff.
5. Fit source density voxels against the reference using the same native
   `locate_maximum` routine underlying fitmap. Maximize the density sum product,
   not a nearest-atom objective. Samples include the low-density outskirts down
   to 0.1% of peak density. If necessary, reproducibly subsample to `max_points`.
6. Refine through the generated resolution pyramid. Retain up to 12 distinct
   coarse placements and six intermediate placements. Keeping multiple
   placements avoids deciding between nearly symmetric alternatives too early.
7. Apply only the best final rigid transform, using the correct parent-model
   coordinate system. Store the old/new local positions for native undo/redo.

The log reports sampled density overlap (or partial density correlation),
elapsed fitting time, fitted size, and the actual resolution range. The whole
score divides estimated grid overlap by the full grid norms.
It is an approximate density-overlap score, not an atom-matched RMSD, probability,
or validated confidence threshold. Sampling and interpolation introduce error;
even self-fit scores need not be exactly one. Do not compare scores across very
different selections or resolutions as though they were calibrated confidence.

## Limitations

- Shapes need substantial shared structure. Partial search handles unequal
  complexes, but factors, missing domains, and different compositions can still
  create ambiguous optima. Fit shared regions explicitly when needed.
- Every fit is rigid. Domain motions, ribosome subunit rotations, and RyR gating
  changes cannot all be superimposed simultaneously with one transform.
- Symmetry-equivalent placements may permute chains. Shape cannot identify which
  chemically identical subunit should be called A. Atom or sequence correspondence
  requires a separate alignment step if needed downstream.
- Multi-start local fitting is not a guaranteed exhaustive global search.
  `extra_starts` can help when principal axes are unreliable, but inspect the fit.
- Instanced structures, fitting a model onto itself, and ancestor/descendant
  source/reference models are rejected. Ordinary transformed parent groups work.
- Each temporary grid is capped at 16 million voxels (64 MB of float32 data).
  Larger requests fail before that grid is allocated; use a coarser resolution
  or smaller fitting region. `max_points` limits fitting samples, not grid size.
- The helper uses ChimeraX Python fitting APIs rather than only command syntax.
  Their compatibility should be rechecked when upgrading ChimeraX.




