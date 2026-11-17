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

## Prediction confidence

Known AlphaFold models are weighted automatically. For local predictions with
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
| `extra_starts 8` | Try more starting orientations |
| `workers 1` | Run atomic-model batches serially |
| `workers 4` | Use four workers |
| `move false` | Report the fit without moving the model |
| `color_fit true` | Color moving atoms/cartoons by local density correlation |
| `correlation true` | Store correlation scores without changing colors |

Run `shapealign help` for examples or `usage shapealign` for all options.
Inspect the result: symmetry and conformational differences can produce
alternative fits.
