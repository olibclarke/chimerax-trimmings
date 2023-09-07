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
mousemode alt left "move picked models"
```
These make some alterations to the mouse modes, setting right-mouse-click-drag to move the front clip plane (shift-right-mouse-click-drag will move both clip planes in the same direction), setting alt-right-mouse-click-drag to adjust the contour level of the currently displayed map, and setting alt/option-left-click-drag to reposition the picked/selected model.


**Aliases**
```
alias cofron cofr centerofview showpivot true
alias cofroff cofr centerofview showpivot false
```
`cofron` sets the center of rotation midway between the two clip planes, and adds a 3D marker for the center of rotation with orthogonal red, green and blue arrows (where R,G,B=X,Y,Z). `cofroff` does the same thing, but switches off the marker.

```
alias symclip cofr centerofview; clip near -$1 far $1 position cofr
```
`symclip` sets the near and far clip planes symmetrically with respect to the center of rotation. So `symclip 5` would set the near clip plane 5 Å from the CoFR, and the far clip plane 5Å in the other direction, giving a 10 Å slab.

```
alias cootmode surface cap false; surface style solid; lighting flat; graphics silhouettes false; style stick; ~rib; color ##num_residues gold; color byhet ; disp;  ~disp @H*; style solvent ball; size ballscale 0.1;  size stickradius 0.07; transparency 70; symclip 10; color ##~num_residues cornflower blue
alias cootmode_mesh surface cap false; surface style mesh; lighting flat; graphics silhouettes false; style stick; ~rib; color ##num_residues gold; color byhet ; disp;  ~disp @H*; style solvent ball; size ballscale 0.1;  size stickradius 0.07; symclip 10; color ##~num_residues #3d60ffff; transparency 50;
```
`cootmode` and `cootmode_mesh` give what I find to be pleasing and performant settings for inspecting atomic models in the context of density maps. By default, hydrogens are not displayed, as I find them irritating under most circumstances.

```
alias map_sphere_15 surface unzone ##~num_residues; sel; close #10000; marker #10000 position cofr; sel ~sel; surface zone ##~num_residues nearAtoms sel distance 15; close #10000
alias map_unsphere surface unzone ##~num_residues
```
`map_sphere_15` will limit display of all maps to a spherical 15Å zone around the center of rotation.

```
alias default_mol_display ~disp; rib; rainbow chain palette RdYlBu-5; lighting soft
```
A nice default display setting for proteins.

```
alias hidemaps surface unzone ##~num_residues; sel; close #10000; marker #10000 position cofr; sel ~sel; surface zone ##~num_residues nearAtoms sel distance 0; close #10000
alias showmaps surface unzone ##~num_residues
```
`hidemaps` and `showmaps` allow quick toggling of the display of the current maps, in order to view or interact with the atomic model underneath. Most useful bound to buttons (see below)


**Shortcut buttons**
```
buttonpanel Shortcuts rows 3 columns 4
buttonpanel Shortcuts add Vol_Viewer command "tool show 'Volume Viewer'"
buttonpanel Shortcuts add Model_Panel command "tool show Models"
buttonpanel Shortcuts add Log command "tool show Log"
buttonpanel Shortcuts add default_disp command "default_mol_display"
buttonpanel Shortcuts add map_sphere command "map_sphere_15"
buttonpanel Shortcuts add map_unsphere command "map_unsphere"
buttonpanel Shortcuts add cofron command "cofron"
buttonpanel Shortcuts add cofroff command "cofroff"
buttonpanel Shortcuts add cootmode command "cootmode"
buttonpanel Shortcuts add mark_cofr command "marker #10000 position cofr"
buttonpanel Shortcuts add hidemaps command "hidemaps"
buttonpanel Shortcuts add showmaps command "showmaps"
```

The first line will create a 3x4 panel of buttons that may be docked to a position of your choosing (I prefer the upper right), with useful shortcuts. Most are documented above; the only one that isn't is `mark cofr`, which places a marker at the center of rotation (useful for measurements of map features). 

I prefer to operate ChimeraX with all panels undocked by default (except for this button panel!), as I find this is the best way to make use of screen real estate on a laptop (which is how I mostly work). So I have a button panel that allows me to quickly access the volume viewer, model panel and log, as well as other shortcuts that I find useful to have on hand. Here is what this looks like in practice:

<img width="1624" alt="image" src="https://github.com/olibclarke/chimeras-trimmings/assets/19766818/7f8f578a-48a3-4474-9526-72da8e7a4472">

And here is what it looks like in action:
https://www.dropbox.com/scl/fi/cdi295s8h7vc5ijf9wr54/chimerax_example_nav.mov?rlkey=mbqn4uk4yodhj1346thyf86oo&dl=0

