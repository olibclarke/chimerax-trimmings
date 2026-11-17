# Licenses and acknowledgments

Original code: [MIT](LICENSE), copyright 2026 Oliver Clarke.
Data retains its upstream terms: [DATA_NOTICE.md](DATA_NOTICE.md).

## Software

| Component | License | Citation |
| --- | --- | --- |
| [ChimeraX](https://www.rbvi.ucsf.edu/chimerax/docs/licensing.html) | Separate noncommercial/commercial application licenses | [Pettersen et al.](https://doi.org/10.1002/pro.3943) |
| [Foldseek](https://github.com/steineggerlab/foldseek/blob/master/LICENSE.md) | GPL v3 | [van Kempen et al.](https://doi.org/10.1038/s41587-023-01773-0) |
| [HMMER](https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE) | BSD-3-Clause and upstream notices | [Eddy](https://doi.org/10.1371/journal.pcbi.1002195) |
| [Gemmi](https://github.com/project-gemmi/gemmi/blob/master/LICENSE.txt) | MPL-2.0 | [Wojdyr](https://doi.org/10.21105/joss.04200) |
| [NumPy](https://numpy.org/doc/stable/license.html) | BSD-3-Clause and distribution notices | [Harris et al.](https://doi.org/10.1038/s41586-020-2649-2) |
| [SciPy](https://scipy.org/about/) | BSD-3-Clause and distribution notices | [Virtanen et al.](https://doi.org/10.1038/s41592-019-0686-2) |

Dependency programs and environments are not bundled.
The SEG implementation is ported from [NCBI's public-domain code](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/src/algo/blast/core/blast_seg.c);
credit Ilya Dondoshansky and Wootton & Federhen.

## Optional services

- [PROSITE](https://prosite.expasy.org/prosite_license.html): CC BY-NC-ND 4.0.
  Commercial use requires separate permission. Use `features false` to skip
  automatic PROSITE scans without that permission. Definitions are downloaded
  unchanged and are not included in the library release.
- [DeepTMHMM](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/): separate
  academic/commercial terms. The MIT [BioLib client](https://pypi.org/project/pybiolib/)
  does not license the predictor. [Citation](https://doi.org/10.1101/2022.04.08.487609).

## Data citations

[ECOD](https://doi.org/10.1093/nar/gkag1132) ·
[Pfam](https://doi.org/10.1093/nar/gkaa913) ·
[InterPro](https://doi.org/10.1093/nar/gkac993) ·
[AlphaFold DB](https://doi.org/10.1093/nar/gkab1061) ·
[AlphaFold](https://doi.org/10.1038/s41586-021-03819-2) ·
[SIFTS](https://doi.org/10.1093/nar/gky1114) ·
[UniProt](https://www.uniprot.org/help/publications) ·
[CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) ·
[wwPDB](https://www.wwpdb.org/about/citing).

Cite original structure papers when discussing individual exemplars.
