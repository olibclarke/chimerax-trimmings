# ChimeraX Trimmings

Shortcuts and tools for model and map work in UCSF ChimeraX.

## Install

1. Download this repository and keep it in a permanent folder.
2. Edit the two `runscript` paths at the top of [chimerax_trimmings.txt](chimerax_trimmings.txt).
3. Paste that file into ChimeraX **Preferences → Startup**, then restart ChimeraX.

All domain tools are included in `chimerax_trimmings.py`.

## Local domain library

1. Install [Foldseek](https://github.com/steineggerlab/foldseek#installation).
   On macOS: `brew install foldseek`.
2. Download the [domain library](https://github.com/olibclarke/chimerax-trimmings/releases/download/domains-2026-11-17/trimmings-ecod-cdd-v295.2.tar.gz).
   Extract it and put the `ecod-cdd` folder inside `.domain_data` in this repository.
3. In ChimeraX:

```chimerax
domains #1 source foldseek
```

For a library stored elsewhere, add `library "/path/to/ecod-cdd"`.
The library contains ECOD/CDD references; some protein families remain unrepresented.

## Domains and sequence features

```chimerax
domains #1
domains #1 source interpro
domains #1 source foldseek
domains #1-3
domains #1 source features
domainbrowser
domains cancel
```

- Expand domain groups to see individual copies.
- Use the model arrows to switch models, or **Previous / Next** to browse domains.
- **Nearby / Selected** identifies domains near the view center or selected residues.
- **Fetch & align domain** loads an exemplar; **Show full context** reveals its full structure.
- Use `color false` to preserve colors, `features false` to skip motifs, or
  `include_tentative false` to hide weaker Pfam matches.

Pfam searches sequences online unless a local database is configured. InterPro
requires an exact database sequence. Foldseek searches locally; descriptions and
exemplars are downloaded. Structural matches suggest a family, not a proven function.

## Shape alignment

```chimerax
shapealign #2 to #1
shapealign #2-5 to #1
shapealign #2 to #1 radius_search
shapealign #2 to #1 color_fit true
```

Aligns models or maps by shape. The first input moves; `to` stays fixed.
See [shape alignment options](SHAPE_ALIGN.md).

## Transmembrane segments

```chimerax
tmsegments #1
tmsegments #1 results /path/predicted_topologies.3line
tmsegments cancel
```

Predicts TM helices and beta strands with [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM).
To run predictions, install its BioLib client once in Terminal:

```bash
python3 -m venv ~/.local/share/chimerax-trimmings/deeptmhmm
~/.local/share/chimerax-trimmings/deeptmhmm/bin/python -m pip install pybiolib
```

Importing an existing `.3line` result does not require the client.

## Other tools

| Command | Use |
| --- | --- |
| `nextmodel` / `prevmodel` | Switch models |
| `nextitem` / `previtem` | Switch models or maps |
| `togglemodels` / `togglemaps` | Hide/restore visible models or maps |
| `cyclelighting` / `cyclemodeldisplay` | Cycle display presets |
| `reorient` | Reorient visible objects |
| `modelmotion` | Control which models can move |
| `exposedhydrophobics` | Find exposed hydrophobic residues |
| `hydrophobicpatches` / `patchbrowser` | Find and browse hydrophobic patches |
| `map_salami` | Make map/model-fit panels |

Use `usage <command>` or `usage aliases` for options.
The **Shortcuts** panel provides common display and navigation controls.
See [mouse defaults](default_mousemodes.md) for mouse bindings.

## Licenses

Original code: [MIT](LICENSE). Data attribution: [DATA_NOTICE.md](DATA_NOTICE.md).
[Third-party terms and citations](THIRD_PARTY.md) apply. PROSITE and DeepTMHMM
have separate commercial-use requirements.
