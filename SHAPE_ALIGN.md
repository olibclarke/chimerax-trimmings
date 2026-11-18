# Shape alignment

Align atomic models or maps by shape:

```chimerax
shapealign #2 to #1
shapealign #2-5 to #1
shapealign #2/A to #1/B
shapealign sel to #1
```

The first input moves; `to` stays fixed. Atomic selections choose fitting atoms,
but the whole model moves. Fits support Undo/Redo.

## Fit a copy

```chimerax
shapealign #2 to #1 copy true
shapealign #2 to #1 radius_search 60 copy true
```

Fits a new copy of the entire moving atomic model, leaving the original's
position, colors, and display unchanged. Atomic selections still restrict the
fitting region. The Log reports the new model ID. Undo removes the fitted copy;
Redo restores it without recalculating the fit. With multiple moving models,
each gets its own copy. Maps are never copied.

`copy true` requires `move true`. Repeating an identical search can find the
same site again: existing fitted copies are not masked or excluded. To fit
another site, reposition the center of rotation and use `radius_search`.

## Multiple non-overlapping copies

```chimerax
shapealign #2 to #1 tiled true copies 6
shapealign #2 to #1 radius_search 80 copies 3
```

Without symmetry, `copies N` searches for up to N independent placements per moving atomic model in a reference
map. With N greater than one, fitted model copies are created and the original
is left unchanged; `copy true` is unnecessary. Tiled searches pool sites across
all tiles before checking overlap, including sites crossing tile boundaries.
All copies from one fit share one Undo/Redo action. `move false` reports the
placements without creating models.

Placements are selected greedily by score. `score_fraction` defaults to 0.8:
remaining placements must score at least 80% of the best candidate. Fewer copies
are returned if the retained candidates do not pass the score/overlap checks.
This is a relative cutoff, not a confidence estimate or a guarantee that every
copy in the map will be found.

Overlap uses selected non-water heavy atoms within `overlap_distance` (default
4 A). With confidence weighting enabled, atom-pair contributions use both atoms'
relative confidence weights; uncertain loops and termini count less, and
zero-weight atoms do not count. Without confidence weighting, all fitting atoms
count equally. `max_overlap` (default 0.1) is the maximum weighted overlapping
fraction, checked in both directions. Only copies generated in this command
are checked against one another, not other models already open.

## Known symmetry

```chimerax
shapealign #2 to #1 tiled true copies 6 symmetry C6
shapealign #2 to #1 tiled true copies 5 symmetry h,8,30
shapealign #2 to #1 tiled true copies 5 symmetry h,8,30 refine_copies true
shapealign #2 to #1 copies 4 symmetry C4 symmetry_axis y symmetry_center 100,120,130
```

With symmetry, only one seed copy is searched and fitted. Candidates that
overlap their implied symmetry neighbours are rejected using the same
confidence-weighted overlap definition. The chosen seed is then expanded to
the requested `copies` count; there is no joint density averaging or repeated
global search for the other copies.

Use native ChimeraX point groups (`Cn`, `Dn`, `T`, `O`, `I`) or `h,rise,twist`.
Rise is in Angstroms and twist in degrees. The short helical form visits repeats
outward from the seed until the requested number of in-map copies is reached.
The native `h,rise,twist,count[,offset]` form restricts expansion to an explicit
repeat range, which must include zero. Point-group orientation suffixes such as
`I,222` use ChimeraX conventions.
Products such as `C2*h,8,30,5,-2` combine point and helical operators.

The default axis is reference-map Z through the grid center (grid index size/2).
Override it with
`symmetry_axis` and `symmetry_center`, in reference-map XYZ coordinates. This is
a geometric default, not an estimate from the density. Supply the actual axis
and center if the map is off-center or differently oriented.

The density score belongs to the seed only; expanded copies inherit it rather
than being independently scored. `score_fraction` does not limit symmetry
expansion. Each output copy must retain at least 99% of its confidence-weighted
fitting atoms inside the map. Fewer copies are returned if the operator range,
map bounds, or overlap checks prevent further expansion. Point groups produce
at most one complete symmetry orbit.

`refine_copies true` optionally performs a final local fit of each expanded
copy, including the seed, against the original map. This can relax exact
symmetry. Overlap is rechecked, reverting to the unrefined placement if needed.
A `radius_search` or `radius_refine` sphere restricts the seed fit, not expansion.
Both symmetry and multiple copies work with tiled search and optional bandpass
scoring. They add computation and currently support atomic-model/map fits only.

## Local regions

Center the view on the region of interest, then use:

```chimerax
shapealign #2 to #1 radius_search
shapealign #2 to #1 radius_search 60
shapealign #2 to #1 radius_refine 30
```

- `radius_search`: search a sphere in the reference map. Without a value, its
  radius is chosen from the moving model's size.
- `radius_refine`: refine the existing fit near the center of rotation.

Long termini can make the automatic radius large. If the CoFR crop exceeds its
memory limit, give an explicit radius, for example `radius_search 50`. For a
whole-map search, omit the radius and use `tiled true` instead.

## Large-map search and bandpass scoring

For a small atomic model in a large map, try the experimental tiled search:

```chimerax
shapealign #2 to #1 tiled true copy true
shapealign #2 to #1 tiled true bandpass true map_resolution 2.2 copy true
```

`tiled true` searches the whole reference map, or the `radius_search` sphere
if supplied. It is not otherwise restricted to the view center. Each
padded submap runs the same global orientation/translation search and native
refinement used for a regional search. Several completed fits from each submap
are pooled and compared against the original map. Padding allows fits to cross
tile boundaries. No symmetry information or previously fitted model is needed.
Spatially separated translation candidates are retained so one bright site
does not fill the entire orientation shortlist before finer-detail scoring.
The input map and surface are not cropped, zoned, duplicated, or restyled.

`bandpass true` ranks those candidates using finer detail than the current
fitting scale. Supply the map's estimated resolution with `map_resolution`
(Angstroms). The score is a model-local Pearson correlation of smooth Gaussian
difference bands, with the same confidence-weighted soft mask at every site.
For example, a 10 A fitting stage scores roughly the 10-5 A band. Scoring never
goes finer than `map_resolution`; at that limit the last available band is
reused. These are soft frequency weights, not sharp Fourier shells. They are
not independent validation data, probabilities, or guarantees of a correct fit.

Both options default to off and currently support moving atomic models with a
reference map only. `bandpass` can also be used with `radius_refine` or
`radius_search`. `tiled` is global and cannot be combined with `radius_refine`.
Without bandpass scoring, tiled candidates are ranked by density overlap on a
common map scale, not by separately normalized submap correlations. These
overlap scores are not correlation coefficients and should not be compared
numerically with non-tiled scores. Bandpass changes candidate ranking, not the
native refinement objective. `resolution` still controls the final fitting
scale, independently of the finer-band cutoff.
The explicit global search is not gated by `fallback_correlation`.

Map searches automatically parallelize sufficiently large batches of candidate
refinements using ChimeraX's native Python and available CPUs minus one.
Large FFTs use the same CPU budget. `workers 1` disables parallel work;
`workers N` sets an explicit maximum, bounded by available CPUs. Small fits stay
serial automatically to avoid startup overhead. Workers use extra memory and
temporary disk space for shared numerical arrays; model changes stay in the
main application.

Tiles are searched sequentially, with `x/y tiles complete` shown throughout
search and refinement. The search still examines all tiles; there is no early
exit based on a guessed quality threshold.

The GUI is serviced between numerical steps, at most ten times per second,
without adding waits or reducing search coverage. Single fits and parallel
atomic-model batches permit window/Log interaction; editing or closing a fitting
input cancels the pending fit. Serial batches and model/map batches repaint but
defer user input while their Undo transaction is open. An individual native
calculation can still briefly occupy the GUI thread.

Whole-map searches repeat the global search for each nonempty submap and can
take substantially longer than a single `radius_search` region.
`tile_size` optionally sets the approximate translation-core edge length in
Angstroms. Its default is model-sized, and it may be reduced to meet the FFT
limit. Smaller tiles reduce per-tile FFT work but increase overlap/repeated work.
Very large models or very fine scoring may still exceed the local grid limit;
restrict fitting atoms or use coarser settings if needed.

Inspect results, especially with repeated motifs, uncertain domain placement,
or weak finer-band agreement. Existing fitted copies are not excluded from a
new search, so repeated searches can return the same site.

## Prediction confidence

Known AlphaFold models, including local AlphaFoldDB CIFs with explicit database
and local pLDDT metadata, are weighted automatically. For other local predictions with
pLDDT in the B-factor field:

```chimerax
shapealign #2 to #1 confidence_moving plddt
shapealign #2 to #1 confidence_reference plddt
shapealign #2 to #1 confidence off
```

Do not use `plddt` for experimental B-factors.

## Options

| Option | Use |
| --- | --- |
| `resolution 8` | Set the fitting scale in Å; default is automatic |
| `objective gradient` | Try the experimental gradient-based objective |
| `tiled true` | Search a reference map with bounded, padded translation tiles |
| `bandpass true map_resolution 2.2` | Rank model/map candidates by finer-band agreement |
| `extra_starts 8` | Try more starting orientations |
| `workers 1` | Run serially |
| `workers 4` | Use four workers |
| `move false` | Report the fit without moving the model |
| `copy true` | Fit a new atomic-model copy; leave the original unchanged |
| `copies 6` | Fit up to six confidence-weighted, non-overlapping copies |
| `symmetry C6` | Fit one seed with implicit-neighbour clash checks, then expand |
| `refine_copies true` | Locally refine symmetry-expanded copies; default: off |
| `color_fit true` | Color moving atoms/cartoons by local density correlation |
| `correlation true` | Store correlation scores without changing colors |

Run `shapealign help` for examples or `usage shapealign` for all options.
Inspect the result: symmetry and conformational differences can produce
alternative fits.
