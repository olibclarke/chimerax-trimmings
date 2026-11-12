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

## Automatic Resolution Range

The whole search pyramid, not just its first level, is generated from the size
of the fitted atom clouds. It uses no GPCR/ribosome presets or atom-count-based
resolution choices. Chain, residue, and CoFR-sphere restrictions are applied
before measuring size.

For each cloud, diagonalize its coordinate covariance. A uniform solid ellipsoid
has variance a^2/5 along semiaxis a, giving equivalent diameters
`2*sqrt(5*variance)`. The geometric mean of the three diameters is an
equivalent-volume diameter D. Axes thinner than 2 Angstroms use a 2-Angstrom
floor so planar clouds still have finite thickness. The smaller of the two
fitted diameters determines the range; using all axes prevents a long thin
object being treated as a large globular particle.

- Automatic fine endpoint: D/12, bounded to 2-8 Angstroms. These bounds avoid
  unnecessary atom-level detail at the fine end and preserve helix-scale
  features for large objects.
- Coarse endpoint: D/4, or the fine endpoint if that is coarser.
- Intermediate levels: geometric spacing, with no greater than a twofold jump.
- `resolution 4` overrides only the fine endpoint. The initial resolution still
  follows object size; it does not revert to a fixed 4x multiplier.
- Local refinement uses the cropped ROI's size and caps its coarse endpoint at
  twice its fine endpoint, since its initial placement is already aligned.

Both endpoints vary with size inside the stated bounds. For example, 1F88/B
versus 2RH1/A selects 11.7 -> 6.74 -> 3.89 Angstroms; 8XP3 versus an 80S model
selects 44.3 -> 25 -> 14.1 -> 8 Angstroms. The log prints the effective fitted
size and actual range. Stats include `resolutions`, `resolution_mode`,
`object_size`, and both input size estimates. Optional correlation-coloring
resolution and neighbourhood remain independent, explicitly set physical scales.

This is a geometric heuristic, not a guarantee that one range is optimal for
every fold. Outlying domains and conformational differences can affect the
ellipsoid. The larger reference still determines the search volume and memory
cost, even when a smaller template determines the resolution range.

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

### Non-Nested Complexes

Neither full input need be a subset of the other: for example, 1F88 contains two
rhodopsins, whereas 2RH1 contains one receptor with an inserted lysozyme domain.
Whole-object coarse blurring can lose the receptor placement even when its
fine-scale density score would be better.

When a partial score is below `min(0.4, fallback_correlation)`, the command can
generate extra poses from selected polymer chains on either side. Each eligible
chain contributes at least 20% but less than 90% of its input's fitted heavy
atoms. At most four chains, largest first, are searched against the full other
input. These bounds avoid tiny-ligand matches and unbounded all-chain pairing.
This is a heuristic, not a guarantee for arbitrary common domains within chains.

The proposed poses enter native refinement directly at the original final
resolution, bypassing the coarse stages that lost the correct basin. Final
ranking uses the same complete input densities, template direction, samples,
and resolution as the existing partial fit, not incomparable chain scores.
The existing pose is retained unless the common score improves. The requested
moving model receives one rigid transform; it is never split or copied. Atom
selections still restrict all seed atoms, and local CoFR mode never uses seeds.
An optional seed that exceeds an allocation guard is skipped with a warning.

The log reports accepted chain-seeded improvements; stats include
`chain_backup_attempts`, `chain_backup_used`, and, when accepted,
`chain_backup_seed` and `chain_backup_initial_score`. Good partial fits skip the
extra search. `partial false` or `fallback_correlation 0` disables it.

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

## Gradient Benchmarks

These opt-in `objective gradient` runs use complete structures in fitting and
scoring, with no automatic partial or chain-seeded fallback. Native Matchmaker
is used only to establish fixed comparison pairs, not to initialize the fit.

| Whole moving -> whole reference | Matchmaker reference | Shape RMSD on those retained pairs | Placement RMS difference from Matchmaker | Runtime |
| --- | --- | --- | --- | --- |
| 1F88 -> 2RH1 | 86 pairs, B/A; MM RMSD 1.245 A | 2.177 A | 1.782 A | 5.2-5.8 s |
| 2RH1 -> 1F88 | 82 pairs, A/A; MM RMSD 1.212 A | 2.017 A | 1.614 A | 5.7 s |
| 8XP3 -> 9S3D | 1,142 pairs, S2/18S; MM RMSD 0.778 A | 0.791 A | 0.144 A | 14.0 s |

The forward GPCR test repeats the same placement across three random starting
poses. Reversing the inputs selects the other rhodopsin copy: both are legitimate
alternatives for this dimer/monomer comparison. These retained-core comparisons
measure placement accuracy, not the gradient objective, which uses whole inputs.
The entire dimer and fusion construct cannot superimpose exactly.

For whole 6QZP -> 9S3D, gradient mode takes about 2.5 seconds with a score of
0.685 and no FFT backup. The equal-size pair passes the gradient cutoff of 0.5.
Using the density mode's 0.7 cutoff instead spent about 17.6 seconds on a broad
search with negligible score improvement; the two metrics require different
heuristics. The subcomplex case needs the broader search and is slower than the
existing partial-density method. No universal speedup or accuracy claim is made.

## Density Benchmarks

Tested in the installed ARM64 ChimeraX Daily (Python 3.14) on this Mac:

- 6QZP: human 80S ribosome, 219,463 loaded atoms.
- 9S3D: human ribosome-nascent-chain complex with associated factors, 218,816
  loaded atoms. The compositions are not identical.
- The equal-size 6QZP/9S3D pair uses the fast whole-shape path. The adaptive
  revision retains the same approximate final overlap score (0.916).
- Nearest-neighbour distances are not used to establish a correct placement:
  a misplaced subcomplex can be close to unrelated ribosomal atoms. Instead,
  native Matchmaker supplies a reference transform and fixed atom pairs.
- A separate default Matchmaker benchmark on the same loaded pair did not finish
  within a 30-second process budget and was stopped. No exact Matchmaker runtime
  or completed-fit speedup ratio is claimed.
- Unit tests include a synthetic fourfold-symmetric shape, arbitrary rigid
  transforms, transformed parent groups, deterministic starts, native command
  parsing, display preservation, dry runs, undo/redo, and memory/input guards.
  Local tests verify spherical scene-coordinate selection, omission of outside
  distractors, current-pose-only refinement, transformed groups, atom-specification
  intersection, invalid/empty spheres, state preservation, and undo/redo.
  Additional tests compare the correlation field with direct weighted Pearson
  calculations, verify undefined regions and optional-path isolation, and check
  multi-model fitting, grouped undo/redo, and rollback on batch failure. Resolution
  tests cover size scaling, rotation invariance, duplicate-atom independence,
  thin and planar clouds, explicit overrides, and ROI-derived ranges. A real RyR
  pair has not yet been benchmarked.

Gradient regressions check discrete integration by parts, self-normalization,
negative target values, unmatched-density normalization, rigid fitting, real FFT
translation search, local-mode isolation, unchanged input arrays and reference,
native command parsing, undo/redo, and rejection of partial fitting. A worse FFT
result cannot replace the existing gradient solution.

Matchmaker comparisons use the same retained pairs for both methods, without
re-fitting after the shape placement:

| Shape-fitting inputs | Automatic range (Angstroms) | Matchmaker reference pairs | Matchmaker RMSD | Shape RMSD on those pairs | Placement RMS difference from Matchmaker |
| --- | --- | --- | --- | --- | --- |
| 1F88/B -> 2RH1/A | 11.7, 6.74, 3.89 | 86 | 1.245 A | 1.978 A | 1.540 A |
| Whole 1F88 -> whole 2RH1 | 15.2, 8.76, 5.06 | 86, from B/A | 1.245 A | 1.991 A | 1.557 A |
| Whole 8XP3 -> whole 6QZP | 44.3, 25, 14.1, 8 | 1,235, from S2/18S | 0.796 A | 0.936 A | 0.492 A |
| Whole 8XP3 -> whole 9S3D | 44.3, 25, 14.1, 8 | 1,142, from S2/18S | 0.778 A | 0.982 A | 0.600 A |

The GPCR result repeated across three arbitrary starting rotations/translations,
taking about 2-3 seconds per fit. Shape fitting uses only the requested B/A chains
in that chain-restricted benchmark. The whole-structure GPCR test uses both
rhodopsin chains and the whole fusion construct, taking about 7 seconds with
chain-seeded backup. It reproduces the placement across three starting poses.
The two rhodopsin copies give alternative valid placements; benchmark baselines
can explicitly include both with `--equivalent-moving-chain A`. The 40S/80S fits use the full
complexes; their ground truth comes only from the shared 18S RNA, not from a
nearest-neighbour match. Partial ribosome searches take roughly 4-5 seconds.

These are not chemically identical structures. In the GPCR comparison,
Matchmaker retained 86 of 258 initially paired positions; the full pairing has
large outliers. The listed 1.978-Angstrom shape RMSD applies only to the retained
core, not every receptor residue. Likewise the full 18S pairings include mobile
regions, and their RMSDs are much larger than the retained-core RMSDs. Listed
placement differences measure the two transforms on the same retained source
points; they are not nearest-neighbour distances or guarantees for other inputs.

## GPU Acceleration

This version is CPU-only: ChimeraX's compiled Gaussian/interpolation routines do
the numerical work, with NumPy handling preprocessing and bookkeeping and the
bundled SciPy handling partial FFT searches and optional correlation filtering. No new
packages are required. The installed ChimeraX Python has no PyTorch, CuPy, or MLX.

GPU acceleration is possible, but not by simply enabling GPU rendering or giving
fitmap a GPU flag. A GPU backend would need GPU-side Gaussian deposition, batched
orientation scoring, interpolation/gradients, and pose refinement. On Apple
Silicon that likely means a Metal-capable backend; CUDA is not available on this
Mac. It would add dependencies and packaging work. At the measured runtime,
improving the CPU search and preprocessing is useful before adding that burden.
Batch alignments or substantially larger searches would be better candidates
for benchmarking a GPU backend.

## Running Tests

Use the Python supplied with ChimeraX, not a generic system Python:

```bash
/Applications/ChimeraX_Daily.app/Contents/bin/python3.14 -B -m unittest discover -s tests -v
```

The tests import the alignment helpers directly from `chimerax_trimmings.py`.
Integration checks also load that main script with a session and verify that
both new and existing commands register without a sibling alignment script.
The tests do not access the network or need a GUI. For the real pair, provide
downloaded mmCIF files to the benchmark:

```bash
/Applications/ChimeraX_Daily.app/Contents/bin/python3.14 -B tests/benchmark_shape_align.py /path/to/6qzp.cif /path/to/9s3d.cif
```

Add `--color-fit` or `--correlation` to measure optional scoring overhead.
`--correlation-resolution` and `--correlation-window` set its scales.
Add `--objective gradient` to either benchmark to test the experimental objective.

For fixed-correspondence validation against Matchmaker:

```bash
/Applications/ChimeraX_Daily.app/Contents/bin/python3.14 -B tests/benchmark_shape_partial.py /path/to/1f88.cif /path/to/2rh1.cif --moving-chain B --reference-chain A --fit-chains
/Applications/ChimeraX_Daily.app/Contents/bin/python3.14 -B tests/benchmark_shape_partial.py /path/to/8xp3.cif /path/to/6qzp.cif --moving-chain S2 --reference-chain S2
```

`--fit-chains` restricts shape fitting to those chains as well. Without it, the
whole complexes are shape-fitted and the selected chains supply only Matchmaker
ground truth. `--resolution` overrides the fine endpoint; `--runs` controls the
number of randomized starting poses.
