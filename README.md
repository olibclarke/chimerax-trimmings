# chimerax-trimmings
Useful aliases and startup settings for UCSF ChimeraX \(https://www.cgl.ucsf.edu/chimerax/). Add these to the "Startup" section of the "Preferences" pane and restart ChimeraX in order to have them take effect. _(Note: Some of these will only work with recent builds of ChimeraX)_. 

For more graphics space on a small screen, open a tool from the Tools menu, right-click its window, and uncheck **Dockable tool**. It will remain a separate window instead of snapping to the main GUI. This preference persists across ChimeraX restarts.

This repository contains:

- `chimerax_trimmings.txt`: Startup commands, aliases, keybindings, and Shortcut buttonpanel setup.
- `chimerax_trimmings.py`: Companion ChimeraX Python commands, including shape alignment of atomic models and density maps. See [SHAPE_ALIGN.md](SHAPE_ALIGN.md).
- `map_salami.py`: Experimental script that , if invoked via runscript, adds a command `map_salami`, which automates the generation of section-by-section map/model density fit figures (can run as `map_salami model #1 selection /A:100-220 map #2 zone 2 patches 3 movie true movie_resolution 1080p segment_size 40 image_height 2000 output_dir ~/Desktop/map_salami_panels`, everything except for map & model is optional; patches controls `MaxComponents` of surface, defaults to `None`; `movie_resolution` can be either 4k or 1080p, defaults to 1080p; `zone` is the surface zone clipping offset, defaults to 2 Å; `segment_size` is the residue range chunk to use, defaults to 40 but will expand to include entire helix/strand or short terminal segments if needed.)

## Installation

1. Put `chimerax_trimmings.txt`, `chimerax_trimmings.py`, and `map_salami.py` somewhere permanent on your machine.
2. Edit the first two lines of `chimerax_trimmings.txt` so the `runscript` paths point to your local copies of the two Python scripts. If upgrading, remove any older separate `runscript` line for `shape_align.py` from Startup preferences; its commands are now included in the main script.
3. Open ChimeraX Preferences and select Startup. On macOS, Preferences is in the UCSF ChimeraX application menu.
4. Paste the contents of your edited `chimerax_trimmings.txt` into the startup commands box.
5. Restart ChimeraX.

## Companion Script

The Python script registers:

- `nextmodel`
- `prevmodel`
- `togglemaps`
- `togglemodels`
- `reorient`
- `modelmotion`
- `exposedhydrophobics`
- `hydrophobicpatches`
- `patchbrowser`
- `shapealign`
- `domains`
- `domains cancel`
- `domainbrowser`
- `tmsegments`
- `tmsegments cancel`

`togglemaps` and `togglemodels` toggle the currently displayed maps/models on/off. `reorient` attempts to orient the view in a sensible direction - aligning the view to the shortest axis of an ellipsoid defined by the visible components - and repeated invocations of the command cycle through alignment to other axes of the ellipsoid. `modelmotion` activates & freezes selected models for motion, in order to mimick the behavior of legacy Chimera, described in more detail below.

## Fast Shape Alignment

Use `shapealign` to superimpose atomic models and density maps without sequence
alignment. Model/model, model/map, map/model, and map/map fitting use the same
syntax:

```chimerax
shapealign #2 to #1
shapealign #2 to #1 objective gradient
shapealign #2-5 to #1
shapealign #2/A to #1/B
shapealign #2 to #1 radius_search
shapealign #2 to #1 radius_search 60
shapealign #2 to #1 radius_refine 30
shapealign #2 to #1 color_fit true
```

The first input moves; `to` identifies one fixed reference. Multiple moving
inputs align independently. Atomic selections restrict fitting atoms but the
whole moving model is repositioned.
The reference is skipped automatically if it is included in the moving range.

Known AlphaFold database models are automatically weighted by pLDDT to reduce
the influence of uncertain tails. For local prediction files, use
`confidence_moving plddt`, `confidence_reference plddt`, or `confidence plddt`
for both inputs. These settings expect pLDDT in the B-factor field, not
experimental B-factors. `confidence off` disables weighting. The whole model
still moves rigidly; maps are unchanged.

Global atomic batches with at least four moving models automatically use the
available logical CPU count minus one, capped by the number of moving models.
Use `workers 1` to run serially or `workers N` to request a process count.
Map fits remain serial. Results retain one Undo action.

Unless a CoFR sphere is requested, maps use full density, independent of contours
and display regions. Map values,
grid metadata, styling, and visibility are unchanged; only the moving owner's
rigid transform changes.

The default `objective density` compares density overlap. Try the experimental
`objective gradient` when bulk overlap gives an implausible orientation.
`radius_search` globally searches a reference-map sphere around the current
center of rotation, keeping the moving selection intact. Omit its value to use
1.5 times the moving input's enclosing radius. `radius_refine` locally
refines an existing alignment using both inputs' density inside that sphere.
Both radii are in Angstroms and cannot be combined. `color_fit true` colors moving atomic
models by local density agreement; it is not available for moving maps.

Type `shapealign` or `shapealign help` to open the Log with a short, formatted
summary. Click **All options** or run `usage shapealign` for complete syntax.
The command is included in `chimerax_trimmings.py` and supports Undo/Redo.
See [Shape Alignment](SHAPE_ALIGN.md) for installation, all options, examples,
and practical limitations.

## Domain Coloring and Annotation

`domains` colors protein residues by domain/family and opens a dockable
**Protein Domains** browser. It is included in `chimerax_trimmings.py`.

```chimerax
domains #1
domains #1 source interpro
domains #1-3
domains #1/A color false
domainbrowser #2
domains cancel
```

**Sources.** Pfam sends each full chain sequence to the
[EMBL-EBI HMMER API](https://hmmer-web-docs.readthedocs.io/en/latest/api.html)
for a new search using fixed, curated family cutoffs. InterPro sends sequence
checksums to the [InterPro Matches API](https://github.com/ProteinsWebTeam/interpro-matches-api)
and retrieves existing matches from databases including SUPERFAMILY, CATH-Gene3D,
Pfam, SMART, and CDD. It can provide broader coverage, but requires the exact
sequence to be indexed; use Pfam for novel or modified sequences.

Matches describe domains, repeats, families, and superfamilies; they may overlap
or leave regions unannotated. Boundaries follow database alignments, including
missing residues and insertion codes, rather than partitioning the structure
geometrically. There is no adjustable confidence threshold. Matching annotations
share colors across models; unassigned target residues are gray. Only residues
in the requested selection receive colors and attributes.

**Browser.** Results persist across searches. Choose a model from the **Model**
dropdown or use **< / >** to cycle through analyzed models. Each model remembers
its current domain. **Show only selected model** hides other atomic structures;
Undo restores visibility. Click a domain or use **Previous / Next** to select
its residues. **Auto-center on domain**, **Focus domain**, and **Restore view**
control the camera. Switching models alone keeps the current view and selection.

Details load automatically: available function descriptions, InterPro/GO
annotations, database links, and representative or example PDB structures.
These describe the family; functions and specificity can differ between members.
**Reload details** retries metadata retrieval without repeating the search.

| Option | Purpose | Default |
| --- | --- | --- |
| `source` | `pfam` search or `interpro` lookup | `pfam` |
| `color` | Color atoms and cartoons; attributes are always assigned | `true` |
| `include_overlaps` | Include competing Pfam assignments | `false` |
| `refresh` | Repeat the search instead of using session-cached results | `false` |
| `gui` | Open or refresh the browser | `true` |
| `wait` | Wait for completion before subsequent script commands | `false` |
| `timeout` | Maximum lookup time in seconds | `600` |

Searches run in the background in graphical sessions; use `wait true` in scripts
that depend on the result. Colors, attributes, and browser results support
Undo/Redo. Reanalyzing a chain replaces its results while retaining other models.
Results are cached per source for the current session; rerun `domains` after
loading a saved session to rebuild the browser.

Residue attributes `domain_name` and `domain_accession` hold the primary match;
`domain_names` and `domain_accessions` retain overlapping matches separated by
` | `. For example: `select #1 & ::domain_accession=PF00071`.


## Transmembrane Helices and Beta Strands

`tmsegments` runs [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM) on BioLib and
annotates both predicted **TM helices** (`M`, orange) and **TM beta strands**
(`B`, cyan). Other residues are gray. Both atom and cartoon colors are updated;
use `color false` to assign attributes while retaining existing colors. Predictions
do not change the model's secondary-structure assignments or coordinates.

```chimerax
tmsegments #1
tmsegments #1/A
tmsegments #1-3 wait true
tmsegments #1 color false
tmsegments #1 refresh true timeout 1800
tmsegments cancel
```

The atom specification identifies **complete chains** to analyze: touching one
residue requests its whole protein chain. Full chain sequences, where available,
are sent to BioLib; coordinates stay local. Identical sequences are predicted
once and results are cached for the current ChimeraX session. `refresh true`
requests fresh predictions. GUI submissions run in the background; `wait true`
waits for completion and is useful in scripts. Headless sessions always wait.
While a remote job is pending, the status bar shows a running message with an
elapsed-seconds counter, updated about once per second. This includes time queued
on BioLib. Completion replaces it with the elapsed time and counts of TM helices
and TM beta strands. Cancellation, timeout, and failure also update the status bar.
The default timeout is 1200 seconds. Cancellation/timeout stops the local client;
an already submitted remote job may continue on BioLib.

### BioLib setup

The command is included in `chimerax_trimmings.py`; no extra startup script is
needed. The BioLib client runs in a separate Python environment. On macOS, create
the default environment in Terminal (using an installed Python 3):

```bash
python3 -m venv ~/.local/share/chimerax-trimmings/deeptmhmm
~/.local/share/chimerax-trimmings/deeptmhmm/bin/python -m pip install pybiolib
```

To use another environment, pass its Python executable explicitly:

```chimerax
tmsegments #1 python /absolute/path/to/venv/bin/python
```

Alternatively set `CHIMERAX_DEEPTMHMM_PYTHON` in the environment from which
ChimeraX is launched. If BioLib requires authentication, log in once from Terminal
with that environment's `biolib login`, or supply `BIOLIB_TOKEN` in ChimeraX's
launch environment. Tokens are not command arguments. See the
[BioLib Python documentation](https://biolib.com/docs/using-applications/python/).
The client interface was checked against pybiolib 1.4.647.

### Import and inspect predictions

An existing DeepTMHMM `predicted_topologies.3line` file can be imported without
installing the BioLib client or submitting a remote job:

```chimerax
tmsegments #1 results /absolute/path/predicted_topologies.3line
select #1 & ::tm_type=helix
select #1 & ::tm_type=strand
select #1 & ::tm_topology=S
select #1/A & ::tm_segment=2
```

The file must contain an **exact sequence match** for each requested chain;
identifiers need not match model or chain names. Extra sequences in a batch file
are allowed. Imported predictions take precedence over cached remote results
for that run and do not replace the remote cache.

| Residue attribute | Meaning |
| --- | --- |
| `tm_topology` | Original DeepTMHMM label: `M`, `B`, `S`, `I`, `O`, or `P` |
| `tm_type` | `helix`, `strand`, or an empty string for other residues |
| `tm_segment` | One-based TM segment number within the chain; `0` outside TM segments |

Signal peptides (`S`) retain their own topology label and are excluded from TM
segment counts. `I`, `O`, and `P` retain DeepTMHMM's inside/outside/periplasm labels;
these are sequence predictions, not an estimated membrane plane in 3D.

The Log shows a clickable list of both segment types, including sequence ranges
and their corresponding structure residue specifications. Missing coordinates,
residue-number offsets, and insertion codes are respected. A segment with no
modeled residues is still reported. Unknown/nonstandard sequence letters are
rejected rather than silently deleted or replaced. If the full polymer sequence
is unavailable, the helper warns before using ChimeraX's inferred chain sequence.

Annotations and coloring support Undo/Redo. Reanalyzing a chain replaces its TM
annotations, including clearing earlier assignments when no TM segments are
predicted; other chains' records remain available. Predictions for a chain that
was deleted or whose sequence/residue mapping changed during submission are
skipped. The combined records are available as `session._trimmings_tm_records`.
Attributes/records are session annotations; this helper does not add a custom
mmCIF export or persist the prediction cache across ChimeraX sessions.


## Exposed Hydrophobic Residues

`exposedhydrophobics` replaces the selection with whole protein residues whose
relative sidechain solvent accessibility exceeds a fractional cutoff.

```chimerax
exposedhydrophobics #1
exposedhydrophobics #1/A minimum_exposure 0.4
exposedhydrophobics #1 residue_names LEU,ILE,VAL,PHE probe_radius 1.4
exposedhydrophobics #1 minimum_exposure 0.25 minimum_area 10
```

| Argument | Meaning | Default |
| --- | --- | --- |
| Atom specification | Candidate residues, e.g. `#1`, `#1/A`, or `sel` | Required |
| `minimum_exposure` | Relative sidechain exposure must exceed this fraction, from 0 to 1 | `0.25` (25%) |
| `minimum_area` | Optional additional sidechain SASA cutoff, in square Angstroms | `0` |
| `residue_names` | Comma-separated residue names to consider | `ALA,VAL,ILE,LEU,MET,PHE,TRP,PRO` |
| `probe_radius` | Solvent probe radius in Angstroms | `1.4` |

The calculation uses all non-water atoms in each candidate's model, including
hidden chains, nucleic acids, ligands, and ions. Restricting candidates to one
chain does not remove burial by other chains in that model. Separate models are
calculated independently, so overlapping alternative structures do not bury each
other. Existing atom radii and any hydrogens present are used. Sidechain area is
summed using ChimeraX's `sideonly` definition, excluding backbone atoms and CA.
Glycine therefore has zero sidechain area.

Relative exposure is the sidechain SASA in the full model divided by its sidechain
SASA with only that residue present. The reference retains the residue's own
backbone and conformation and uses the same atom radii, hydrogens, and probe radius
as the full-model calculation. This compares exposed fractions rather than imposing
the same absolute area on large and small sidechains. It is an isolated-residue
reference, not a published maximum-accessibility table or tripeptide reference.
Residues with zero reference area have zero relative exposure. Incomplete
sidechains are normalized using only the atoms present, so a highly exposed
fragment does not imply that the missing sidechain would also be exposed.

The command assigns three attributes to every protein residue in the models it
calculates: `sidechain_sasa` (square Angstroms), `sidechain_sasa_reference` (square
Angstroms), and `relative_sidechain_sasa` (0 to 1). They support related selections
without recalculating, for example:

```chimerax
select #1 & :ARG,LYS,ASP,GLU & ::relative_sidechain_sasa>0.25
```

Rerun `exposedhydrophobics` after changing coordinates, radii, model composition,
or the probe radius to refresh those values. It does not create a surface or
change model visibility, coloring, or coordinates. Selection changes are undoable.

## Binder-Accessible Hydrophobic Patches

`hydrophobicpatches` screens for exposed hydrophobic clusters that a globular
binder, approximated as a rigid sphere, could contact. It selects the reported
patches and, by default, gives each patch a distinct atom/cartoon color.

```chimerax
hydrophobicpatches #1 binder_radius 15
hydrophobicpatches #1/A binder_mass 15 context #2
hydrophobicpatches #1 binder_radius 15 samples 128 highlight false
select #1 & ::hydrophobic_patch_rank=1
```

Supply exactly one of `binder_radius` (Angstroms) and `binder_mass` (kDa).
Mass is converted to an equivalent sphere at an assumed protein density of
1.35 g/cm^3; e.g. 15 kDa corresponds to about 16.4 Angstroms radius.

| Option | Meaning | Default |
| --- | --- | --- |
| Atom specification | Candidate residues; the full containing model remains the steric environment | Required |
| `binder_radius` | Binder sphere radius in Angstroms | Supply this or `binder_mass` |
| `binder_mass` | Approximate globular binder mass in kDa | Supply this or `binder_radius` |
| `context` | Additional blocker atoms from separate models, e.g. glycans or assembly partners | None |
| `minimum_exposure` | Minimum relative sidechain exposure with the solvent probe | `0.25` |
| `probe_radius` | Solvent probe radius in Angstroms, distinct from the binder radius | `1.4` |
| `minimum_residues` | Minimum hydrophobic residues in a patch | `3` |
| `contact_distance` | Maximum atom-surface to binder-surface gap in Angstroms | `3` |
| `link_distance` | Maximum distance between exposed sidechain atom centers linking residues | `6` |
| `samples` | Deterministic sphere-placement directions per exposed heavy sidechain atom | `64` |
| `max_patches` | Maximum number of distinct, non-overlapping patches reported | `10` |
| `approach` | Require a straight, collision-free approach from outside | `true` |
| `highlight` | Color patch atoms and cartoons | `true` |
| `gui` | Automatically open or refresh the patch browser in GUI sessions | `true` |
| `residue_names` | Hydrophobic residue names | `ALA,VAL,ILE,LEU,MET,PHE,TRP,PRO` |

The screen first applies the same relative sidechain-exposure definition as
`exposedhydrophobics`, using the solvent probe (1.4 Angstroms by default). It then samples sphere
placements tangent to exposed sidechain atoms. Every reported patch must have
one shared clash-free placement near all its reported residues, and its residues
must form a connected cluster at `link_distance`. With `approach true`, the sphere
must also have an unobstructed straight outward path to infinity. This is stricter
than testing whether individual residues are independently accessible.

All non-water atoms in the containing model contribute to exposure and steric
occlusion, including hidden chains, ligands, and modeled glycans. Each candidate
model is screened independently; `context` adds blockers from separate models.
Patches may span chains within a model, but are not combined across models.

Patches are ranked by the summed water-probe SASA of hydrophobic sidechain atoms
near the sampled binder sphere. This is nearby exposed area, not predicted buried
interface area or binding energy. The highest-ranked patch is kept first, then
lower-ranked overlapping alternatives are omitted. Residue attribute
`hydrophobic_patch_rank` is the rank, or zero for unassigned protein residues in
the analyzed models. Placement centers and approach directions in scene coordinates
are available in `session._trimmings_hydrophobic_patches`.

Use these as candidate epitopes for inspection, not predictions of affinity,
specificity, or binder design success. A sphere omits binder shape, sidechain
flexibility, missing glycans, and unmodeled membrane/assembly context. Finite
direction sampling can miss accessible placements; increase `samples` to refine
the screen. A failed straight-approach test does not rule out a curved approach.
Incomplete target sidechains retain the exposure limitations described above.

Color and selection changes are undoable together; undo a highlight before
comparing another set of parameters. `highlight false` leaves colors untouched.
The command uses NumPy and SciPy supplied with ChimeraX, creates no model copies
or probe models, and leaves coordinates, map contours, and display masks unchanged.

### Patch Browser

In a graphical ChimeraX session, a search automatically opens a dockable
**Hydrophobic Patches** panel. Its list-and-navigation layout follows ResNotes,
but ResNotes does not need to be installed.

- Click a patch in the ranked list, or use **Previous** / **Next**, to select it.
  Navigation wraps at either end.
- **Auto-center on patch** fits the patch and adjusts local clipping. Uncheck it
  to change selection while keeping the current camera.
- **View from binder approach** orients focused views along the sampled approach
  direction. Uncheck it to retain the current orientation when centering.
- **Focus patch** centers the current patch manually.
- **Restore view** restores the camera, clipping, and rotation center captured
  when the panel opened. It does not restore or change model positions.
- **Refresh list** loads the latest search results. New searches refresh the
  existing panel, and closed-model results are removed automatically.

The details show the residue specification, nearby hydrophobic area, and binder
sphere radius. Focusing follows a moved target's current transform; rescreen after
structural edits or changes in the target's position relative to blocker models.
Hidden models are left hidden, and the panel does not change display masks.

Reopen the latest results with `patchbrowser`. Suppress automatic opening with:

```chimerax
hydrophobicpatches #1 binder_mass 15 gui false
```

Without a GUI, searches remain command-only. Browser state is not saved in
session files; rerun the search after reopening a session.

## Modelmotion

The `modelmotion` commands provide legacy Chimera-style activation independently
of atom/model selection. Structures, maps, and independent surfaces are supported.
Groups expand to their movable members; surfaces owned by a structure or map move
with that owner. Specify the owner rather than its surface submodel.

E.g.
```python
modelmotion active #1
modelmotion unfreeze #3
modelmotion freeze #1
modelmotion status
modelmotion all
modelmotion off
```

`active` replaces the movement set, `unfreeze` adds models, and `freeze` removes
models. If `freeze` is the first command, all other current models remain active.
Newly opened models are active by default; existing frozen models stay frozen.
An empty active set prevents rotation/translation; clearing selection has no effect.

Mouse and trackpad input work simultaneously and use the same active/frozen model
set. With the Trimmings startup bindings and ChimeraX's standard middle-button and
wheel bindings, mouse controls are:

| Mouse input | Behavior |
| --- | --- |
| Left drag | Rotate active models together |
| Middle-button drag | Translate active models in the screen plane |
| Control + middle-button drag | Translate active models along screen Z |
| Scroll wheel | Keep the existing whole-view zoom binding |
| Right drag | Keep the Trimmings clipping control |
| Alt + right drag | Keep the Trimmings contour control |

Custom mouse bindings are respected: existing model rotation/translation bindings
are redirected to the active set when `modelmotion` is enabled. Selected-model
movement bindings then use the active set independently of selection.

| Trackpad gesture | Movement |
| --- | --- |
| Two-finger swipe | Rotate the active models together |
| Two-finger twist | Rotate active models around the screen-normal axis |
| Three-finger swipe | Translate active models in the screen plane |
| Control + three-finger swipe | Translate active models along screen Z |
| Pinch | Keep the existing zoom binding |

Trackpad gestures must be enabled in ChimeraX and available to it from macOS.
Rotation uses the current ChimeraX center of rotation, including the pivot set by
`cofron`, `cofroff`, or an explicit `cofr` position. Mouse bindings that currently
rotate/translate models are redirected to the same active set; selection, clipping,
and contour bindings retain their existing behavior. Trackpad movement is undoable
after a short pause. `modelmotion off` restores the previous mouse and trackpad
bindings changed on activation.
Freezing applies to these movement modes; commands such as `fitmap` can still move
models. Activation is temporary session state and is not saved in session files.

## Keybindings

- `F1`: Previous model
- `F2`: Next model
- `F3`: Center selection
- `F4`: View all
- `F5`: Toggle maps
- `F6`: Toggle models
- `F7`: Cycle lighting presets
- `F8`: Cycle model display options
- `F9`: Previous map (Ignores models)
- `F10`: Next map (Ignores models)

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
- If you move the scripts, update both `runscript` paths in your Startup preferences.

## Example 
The shortcut button panel provides quick access to the Volume Viewer, Model Panel, Log, and display commands. It can stay docked while other tool windows are undocked to leave more room for graphics:

<img width="1624" alt="image" src="https://github.com/olibclarke/chimeras-trimmings/assets/19766818/7f8f578a-48a3-4474-9526-72da8e7a4472">

And here is what it looks like in action:
https://www.dropbox.com/scl/fi/cdi295s8h7vc5ijf9wr54/chimerax_example_nav.mov?rlkey=mbqn4uk4yodhj1346thyf86oo&dl=0
