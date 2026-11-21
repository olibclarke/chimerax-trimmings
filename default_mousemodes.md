_**Mouse modes:**_

**Right-click-drag**: Moves front clip plane

**Shift**-right-click-drag: Moves both clip planes

**Alt/option**-right-click-drag: Adjusts map threshold

**Alt/option**-left-click-drag: Translates the selected model(s) in x-y plane

**Alt/option**-left-click-**Ctrl**-drag: Translates the selected model(s) along z. Order important! Must press Ctrl after clicking, and z-translation only happens with vertical mouse motion.

**Shift+Alt/option**-left-click-drag: Rotates the selected model(s)

**Ctrl+Alt**-left-click: Picks & colors clicked blob (launches "Measure and Color Blobs" tool).

**Ctrl**-left-click: Picks the frontmost atom, bond, cartoon residue, or atom-linked shape through transparent contours and surfaces. Opaque surfaces block the click.

**Ctrl**-left-drag: Uses ChimeraX's native rectangle selection for atoms, bonds, cartoon residues, and atom-linked shapes, including those behind other objects. It respects clipping and does not select maps or surface patches.

Neither gesture highlights transparent molecular surface patches.

**Ctrl+Shift**-left-click or drag: Toggles the same selection.

**Ctrl+Alt+Shift**-left-click or drag: Uses ChimeraX's standard selection to toggle objects, including maps and surfaces.

### Restore standard ChimeraX selection

Run these commands in ChimeraX to restore standard object picking for Ctrl-click and Ctrl-drag, including maps and surfaces.

```chimerax
mousemode control left select
mousemode control shift left select
mousemode control alt shift left none
```

To keep this behavior across restarts, add these commands to the end of your Startup preferences. Reloading Trimmings or clicking `reset_mouse_and_help` re-enables visible-structure picking for the current session.

_**Function keys:**_

**F1** - previous model

**F2** - next model

**F3** - center selection

**F4** - view all

**F5** - Toggle maps

**F6** - Toggle models

**F7** - Cycle lighting presets

**F8** - Cycle model display options

**F9** - Previous map (Ignores models)

**F10** - Next map (Ignores models)

**For a description of chimerax-trimmings as a whole, see here:** [Readme](README.md)
