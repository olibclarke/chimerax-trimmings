# chimerax-trimmings
Useful aliases and startup settings for UCSF ChimeraX \(https://www.cgl.ucsf.edu/chimerax/). Add these to the "Startup" section of the "Preferences" pane and restart ChimeraX in order to have them take effect. _(Note: Some of these will only work with recent builds of ChimeraX)_. 

I would also recommend (if you are using a laptop) going through each tool in the "Tools" menu, launching it, right clicking on the relevant window and **unchecking** "Dockable tool" - this will ensure that each tool acts as a separate window and does not try to "snap" to the main GUI when you move it around, which I find undesirable on smaller screens (these settings will be persistent across restarts of ChimeraX).

This repository contains:

- `chimerax_trimmings.txt`: Startup commands, aliases, keybindings, and Shortcut buttonpanel setup.
- `chimerax_trimmings.py`: Companion ChimeraX Python script for a couple of things that are difficult or not possible to implement in ChimeraX command language.

## Installation

1. Make a local copy of `chimerax_startup.py` (somewhere permanent).
2. Edit the first line of `chimerax_trimmings.txt` so the `runscript` path points to your local copy of `chimerax_trimmings.py`.
3. In ChimeraX on macOS, open `UCSF ChimeraX -> Preferences... -> Startup` (I think under Favorites on Windows?)
4. Paste the contents of your edited `chimerax-trimmings.txt` into the startup commands box.
5. Restart ChimeraX.

## Companion Script

The Python script registers:

- `nextmodel`
- `prevmodel`
- `togglemaps`
- `togglemodels`

`togglemaps` and `togglemodels` toggle the currently displayed maps/models on/off. 

## Keybindings

- `F1`: Previous model
- `F2`: Next model
- `F3`: Center selection
- `F4`: View all
- `F5`: Toggle maps
- `F6`: Toggle models
- `F7`: Cycle lighting presets

## Button Panel

The startup file creates a `Shortcuts` button panel with quick access to commonly used tools/commands.

It currently includes buttons for:

- `Vol_Viewer`: open the Volume Viewer tool
- `Model_Panel`: open the Model Panel
- `Log`: open the Log viewer
- `default_disp`: Set some nice defaults for model/map display
- `map_sphere`: Zone maps to a 15 A sphere around the center of rotation
- `map_unsphere`: Undo map zoning
- `cofron`: Make sure the center of rotation is centered between the clip planes and show the center-of-rotation pivot marker
- `cofroff`: Hide the center-of-rotation pivot marker while keeping the CoFR centered between the clip planes.
- `cootmode`: Apply a Coot-style display preset for models/maps
- `mark_cofr`: Place a marker at the current center of rotation
- `togglemaps`: Hide shown maps, then restore the same set on the next use
- `togglemodels`: Hide shown atomic models, then restore the same set on the next use
- `reset_mouse_and_help`: Restore the default mouse bindings and open the mouse-mode help page
- `previous_model`: Step to the previous model in the Model Panel
- `next_model`: Step to the next model in the Model Panel

## Alias Reference

_**Note: for all of these aliases, on the Chimerax command line you can type "help `cofron`" or "usage `cofron`" for example to get a short description of usage of the indicated alias in the Log**_

### `cofron`

Centers the rotation point between the clip planes and shows the pivot marker.

Usage: `cofron`

```chimerax
alias cofron cofr centerofview showpivot 7,0.25
alias usage cofron synopsis "Centers the rotation point on the view and show the pivot"
```

### `cofroff`

Recenters the rotation point between the clip planes and hides the pivot marker.

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
alias usage symclip synopsis "Set symmetric clip planes around the center of rotation" $1 "Half-width"
```

### `selmodel`

Selects all atomic models.

Usage: `selmodel`

```chimerax
alias selmodel sel ##num_residues
alias usage selmodel synopsis "Select atomic models"
```

### `selmap`

Selects all non-atomic models such as maps.

Usage: `selmap`

```chimerax
alias selmap sel ~##num_residues
alias usage selmap synopsis "Select non-atomic models such as maps"
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

Shows a CA trace for proteins or phosphate backbone trace for nucleic acids, plus protein sidechains and nucleic-acid bases.

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

Marks selected residues as helix.

Usage: `helix <residue-spec>`

```chimerax
alias helix setattr $1 res is_helix true
alias usage helix synopsis "Mark residues as helix" $1 "residue-spec"
```

### `strand`

Marks selected residues as strand.

Usage: `strand <residue-spec>`

```chimerax
alias strand setattr $1 res is_strand true
alias usage strand synopsis "Mark residues as strand" $1 "residue-spec"
```

### `coil`

Clears helix and strand assignment for selected residues.

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

### `selside`

Selects protein sidechains and nucleic-acid bases for the given atom specification.

Usage: `selside <atom-spec>`

```chimerax
alias selside sel sidechain & $1 | ~backbone & nucleic & $1
alias usage selside synopsis "Select protein sidechains and nucleic-acid bases" $1 "atom-spec"
```

### `dispside`

Displays protein sidechains and nucleic-acid bases for the given atom specification.

Usage: `dispside <atom-spec>`

```chimerax
alias dispside sel sidechain & $1 | ~backbone & nucleic & $1; disp sel
alias usage dispside synopsis "Display protein sidechains and nucleic-acid bases" $1 "atom-spec"
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
- If you move the repository, update the `runscript` path in the first line of your Startup preferences panel.

## Example 
I prefer to operate ChimeraX with all panels undocked by default (except for this button panel!), as I find this is the best way to make use of screen real estate on a laptop (which is how I mostly work). So I have a button panel that allows me to quickly access the volume viewer, model panel and log, as well as other shortcuts that I find useful to have on hand. Here is what this looks like in practice:

<img width="1624" alt="image" src="https://github.com/olibclarke/chimeras-trimmings/assets/19766818/7f8f578a-48a3-4474-9526-72da8e7a4472">

And here is what it looks like in action:
https://www.dropbox.com/scl/fi/cdi295s8h7vc5ijf9wr54/chimerax_example_nav.mov?rlkey=mbqn4uk4yodhj1346thyf86oo&dl=0

