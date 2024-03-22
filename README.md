# chimerax-trimmings
Useful aliases and startup settings for UCSF ChimeraX \(https://www.cgl.ucsf.edu/chimerax/). Add these to the "Startup" section of the "Preferences" pane and restart ChimeraX in order to have them take effect. _(Note: Some of these will only work with the daily build)_. 

I would also recommend (if you are using a laptop) going through each tool in the "Tools" menu, launching it, right clicking on the relevant window and **unchecking** "Dockable tool" - this will ensure that each tool acts as a separate window and does not try to "snap" to the main GUI when you move it around, which I find undesirable on smaller screens (these settings will be persistent across restarts of ChimeraX).

Here is my current startup section for ChimeraX, broken down with descriptions of each section:

**Display settings:**
```
camera ortho
cofr centerofview
volume defaultvalues limitVoxelCount false voxelLimitForPlane 1000000000000000 voxelLimitForOpen 1000000000000000 saveSettings true
```
`camera ortho` sets the default camera mode to orthographic (which I prefer to perspective, and which allows for accurate scale bars). `cofr centerofview` sets the center of rotation midway between the two clip planes (`near` and `far`). The last line sets some defaults for displaying new volumes (this currently ony works in the daily build), setting the default step to 1 and switching off plane display of large volumes.


**Mouse mode settings**
```
mousemode right clip
mousemode alt right contour
mousemode alt left "translate selected models"
mousemode shift left "rotate selected models"
mousemode alt control left "pick blobs"
```
These make some alterations to the mouse modes:

Right-click-drag: Moves front clip plane

**Shift**-right-click-drag: Moves both clip planes

**Alt/option**-right-click-drag: Adjusts map threshold

**Alt/option**-left-click-drag: Translates the _selected_ model(s) in x-y plane

**Shift+Alt/option**-left-click-drag: Rotates the _selected_ model(s)

**Ctrl+Alt**-left-click: Picks & colors clicked blob (launches "Measure and Color Blobs" tool).

I would also suggest enabling "Mouse clipping enables _screen planes_" in the "clipping" section of preferences, rather than the default of altering _scene planes_, which in my experience can lead to unexpected results.

One other note - certain Mac accessibility settings may interfere with trackpad behavior in ChimeraX, so if you are seeing something weird, check the accessibility settings (e.g. the "use trackpad gesture to zoom" setting cn cause problems).


**Aliases**
```
alias cofron cofr centerofview showpivot 7,0.25
alias cofroff cofr centerofview showpivot false
```
`cofron` sets the center of rotation midway between the two clip planes, and adds a 3D marker for the center of rotation with orthogonal red, green and blue arrows (where R,G,B=X,Y,Z). `cofroff` does the same thing, but switches off the marker.  

```
alias symclip cofr centerofview; clip near -$1 far $1 position cofr
```
`symclip` sets the near and far clip planes symmetrically with respect to the center of rotation. So `symclip 5` would set the near clip plane 5 Å from the CoFR, and the far clip plane 5Å in the other direction, giving a 10 Å slab.  

```
alias cootmode set bgColor black; surface cap false; surface style solid; lighting flat; graphics silhouettes false; style stick; ~rib; color ##num_residues gold; color byhet ; disp;  ~disp @H*; style ions ball; style solvent ball; size ballscale 0.2;  size stickradius 0.07; transparency 70; cofr centerofview; clip near -10 far 10 position cofr; color ~##num_residues cornflower blue
alias cootmode_mesh surface cap false; surface style mesh; lighting flat; graphics silhouettes false; style stick; ~rib; color ##num_residues gold; color byhet ; disp;  ~disp @H*; style solvent ball; style ions ball; size ballscale 0.2;  size stickradius 0.07; cofr centerofview; clip near -10 far 10 position cofr; color ~##num_residues #3d60ffff; transparency 50
```
`cootmode` and `cootmode_mesh` give what I find to be pleasing and performant settings for inspecting atomic models in the context of density maps. By default, hydrogens are not displayed, as I find them irritating under most circumstances.  

```
alias ca_and_sidechains ~rib $1; ~surf $1; ~disp $1; disp @CA&protein&$1; disp @P&nucleic&$1; style $1 stick; disp sidechain&$1; disp ~backbone&nucleic&$1; size stickradius 0.1; size pseudobondradius 0.1
alias ca_trace ~rib $1; ~surf $1; ~disp $1; disp @CA&protein&$1; disp @P&nucleic&$1; style $1 stick; size stickradius 0.1; size pseudobondradius 0.1
```
`ca_and_sidechains` will display the selected model (executed as e.g. `ca_and_sidechains #1`) as a C-alpha (or phosphate for nucleic acids) backbone with attached sidechains/bases. `ca_trace` will do the same, just without the sidechains/bases.  

```
alias map_sphere_15 volume unzone ##~num_residues; sel; close #10000; marker #10000 position cofr; sel ~sel; volume zone ~##num_residues nearAtoms sel minimalBounds true range 15; close #10000
alias map_unsphere volume unzone ~##num_residues
```
`map_sphere_15` will limit display of all maps to a spherical 15Å zone around the center of rotation.  

```
alias default_mol_display ~disp; rib; rainbow chain palette RdYlBu-5; lighting soft
```
A nice default display setting for proteins.  

```
alias hidemaps surface unzone ##~num_residues; sel; close #10000; marker #10000 position cofr; sel ~sel; surface zone ~##num_residues nearAtoms sel distance 0; close #10000
alias showmaps surface unzone ~##num_residues
```
`hidemaps` and `showmaps` allow quick toggling of the display of the current maps, in order to view or interact with the atomic model underneath. Most useful bound to buttons (see below)  

```
alias caps_off surface cap false
alias caps_on surface cap true
```
I like having aliases to quickly switch surface caps on and off, as displaying caps can be useful for figures, but dramatically slows down rotation/translation of large maps (especially with `full` lighting).  

```
alias selbetween ks ri
```
Selects all residues (inclusive) between the selected residues. Temporary I think (until this is officially added as a selection mode in ChimeraX).  

```
alias helix setattr $1 res is_helix true
alias strand setattr $1 res is_strand true
alias coil setattr $1 res is_strand false; setattr $1 res is_helix false
```
Aliases for changing secondary structure assignment by selection. For example, `helix sel` will set the selected residues as helical, even if they were automatically assigned as strand or coil. Use with caution, and always check _why_ the automated assignment is failing.  

```
alias rock_movie cofr showpivot false; movie record; rock y 30; wait 136; movie encode ~/Desktop/rock_movie.mp4; stop
```
Makes a simple rocking movie, which can be seamlessly looped, e.g. for use in a presentation to show map/model fit. Hides center of rotation indicator if shown.

```
alias open_vseries "open browse vseries true"
```
Creates a volume series using files selected from the file browser GUI.

**Shortcut buttons**
```
buttonpanel Shortcuts rows 3 columns 5
buttonpanel Shortcuts add Vol_Viewer command "tool show 'Volume Viewer'"
buttonpanel Shortcuts add Model_Panel command "tool show Models"
buttonpanel Shortcuts add Log command "tool show Log"
buttonpanel Shortcuts add default_disp command "default_mol_display"
buttonpanel Shortcuts add map_sphere command "map_sphere_15"
buttonpanel Shortcuts add map_unsphere command "map_unsphere"
buttonpanel Shortcuts add cofron command "cofron"
buttonpanel Shortcuts add cofroff command "cofroff"
buttonpanel Shortcuts add cootmode command "cootmode"
buttonpanel Shortcuts add mark_cofr command "marker #20000 position cofr"
buttonpanel Shortcuts add hidemaps command "hidemaps"
buttonpanel Shortcuts add showmaps command "showmaps"
buttonpanel Shortcuts add reset_mouse command "mousemode alt right contour ; mousemode right clip ; mousemode alt left 'translate selected models' ; mousemode shift left 'rotate selected models' ; mousemode alt control left 'pick blobs'; help https://github.com/olibclarke/chimerax-trimmings/blob/main/default_mousemodes.md"
```

The first line will create a 3x5 panel of buttons that may be docked to a position of your choosing (I prefer the upper right - remember to right-click and check "Save tool position" & "Dockable tool" to make this persist), with useful shortcuts. Most are documented above; the only one that isn't is `mark cofr`, which places a marker at the center of rotation (useful for measurements of map features).

There are two additional buttons that I use **which will only work with the daily build**, when sequential model display controls are activated in the model panel:

```
buttonpanel Shortcuts add previous_model command "prevmodel"
buttonpanel Shortcuts add next_model command "nextmodel"
```
These nextmodel/prevmodel buttons will allow for quick switching between models, without needing to activate the model panel.

I prefer to operate ChimeraX with all panels undocked by default (except for this button panel!), as I find this is the best way to make use of screen real estate on a laptop (which is how I mostly work). So I have a button panel that allows me to quickly access the volume viewer, model panel and log, as well as other shortcuts that I find useful to have on hand. Here is what this looks like in practice:

<img width="1624" alt="image" src="https://github.com/olibclarke/chimeras-trimmings/assets/19766818/7f8f578a-48a3-4474-9526-72da8e7a4472">

And here is what it looks like in action:
https://www.dropbox.com/scl/fi/cdi295s8h7vc5ijf9wr54/chimerax_example_nav.mov?rlkey=mbqn4uk4yodhj1346thyf86oo&dl=0

