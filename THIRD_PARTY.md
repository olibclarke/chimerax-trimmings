# Licenses and acknowledgments

Original Trimmings code is [MIT licensed](LICENSE), copyright 2026 Oliver Clarke.
Database contents retain their upstream terms; see [DATA_NOTICE.md](DATA_NOTICE.md),
which is included in generated library archives. Terms checked 17 November 2026.

## Data

| Resource | Terms / acknowledgment |
| --- | --- |
| [ECOD](https://doi.org/10.1093/nar/gkag1132) | CC BY 4.0; retain release, domain UID, attribution and change notices. |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/faq), [UniProt](https://www.uniprot.org/help/license) | CC BY 4.0; retain source accessions and attribution. |
| [Pfam/InterPro](https://interpro-documentation.readthedocs.io/en/latest/license.html) | CC0 1.0 for InterPro, Pfam, PRINTS and SFLD downloadable data; acknowledge the resources. Other member databases retain their own terms. |
| [wwPDB](https://www.wwpdb.org/about/usage-policies), [PDBe/SIFTS](https://www.ebi.ac.uk/pdbe/about/public-data-access-statement) | CC0; credit the resources and original structure authors. |
| [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) | NCBI and original member-database terms; do not apply a blanket CC0 claim to imported annotations. |

Archives retain source metadata, coverage reports and processing provenance.
Cropping, C-alpha conversion, coordinate translation, reannotation and database
encoding are modifications, not official upstream releases. The published Pfam
structural snapshot is partial: 46,654 of 52,484 folded-scope families; its full
sequence-profile bundle contains all 53,275 profiles in Pfam 39.0.

## Software

| Component | License / use |
| --- | --- |
| [Foldseek](https://github.com/steineggerlab/foldseek/blob/master/LICENSE.md) | GPL v3; separately installed executable. |
| [HMMER](https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE) | BSD-3-Clause and upstream component notices; separately installed executable. |
| [Gemmi](https://github.com/project-gemmi/gemmi/blob/master/LICENSE.txt) | MPL-2.0; installed builder dependency. |
| [NumPy](https://numpy.org/doc/stable/license.html), [SciPy](https://scipy.org/about/) | BSD-3-Clause project licenses; installed numerical dependencies. |
| [curl](https://curl.se/docs/copyright.html), Python | curl license; Python/PSF terms; installed tools. |
| [NCBI SEG](https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/master/src/algo/blast/core/blast_seg.c) | Python port of the public-domain NCBI implementation; credit Ilya Dondoshansky and Wootton & Federhen. |

No dependency binaries or environments are included in the data archives.
Using Foldseek to generate data does not itself make that data GPL; see the
[GNU output FAQ](https://www.gnu.org/licenses/gpl-faq.html#GPLOutput).
Redistributing dependency binaries or copying their code requires compliance
with their actual licenses and notices.

## Application and optional features

- [ChimeraX](https://www.rbvi.ucsf.edu/chimerax/docs/licensing.html) has separate
  noncommercial/commercial application licenses; it is not redistributed here.
- [PROSITE](https://prosite.expasy.org/prosite_license.html) definitions are
  CC BY-NC-ND 4.0. Commercial use requires separate permission. The helper retains
  the downloaded original file unchanged; definitions and caches are excluded
  from releases. Automatic sequence features include PROSITE; use `features false`
  without the necessary permission. Pfam motifs and the local SEG implementation
  are separate from PROSITE.
- [DeepTMHMM](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/) has
  separate academic/commercial terms. The MIT [pybiolib](https://pypi.org/project/pybiolib/)
  client does not license the predictor. We bundle no predictor weights/container.

## Cite what you use

- ECOD: [Schaeffer et al., doi:10.1093/nar/gkag1132](https://doi.org/10.1093/nar/gkag1132).
- Foldseek: [van Kempen et al., 2023](https://doi.org/10.1038/s41587-023-01773-0).
- HMMER: [Eddy, 2011](https://doi.org/10.1371/journal.pcbi.1002195).
- Pfam: [Mistry et al., 2021](https://doi.org/10.1093/nar/gkaa913).
- InterPro: [Paysan-Lafosse et al., 2023](https://doi.org/10.1093/nar/gkac993).
- AlphaFold DB: [Varadi et al., 2022](https://doi.org/10.1093/nar/gkab1061);
  AlphaFold method: [Jumper et al., 2021](https://doi.org/10.1038/s41586-021-03819-2).
- SIFTS: [Dana et al., 2019](https://doi.org/10.1093/nar/gky1114).
- [wwPDB citation guidance](https://www.wwpdb.org/about/citing), individual
  structure papers, [UniProt](https://www.uniprot.org/help/publications) and
  [CDD/member resources](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml), as applicable.
- ChimeraX: [Pettersen et al., 2021](https://doi.org/10.1002/pro.3943).
- Gemmi: [Wojdyr, 2022](https://doi.org/10.21105/joss.04200).
- NumPy: [Harris et al., 2020](https://doi.org/10.1038/s41586-020-2649-2);
  SciPy: [Virtanen et al., 2020](https://doi.org/10.1038/s41592-019-0686-2).
- [PROSITE references](https://prosite.expasy.org/prosite_details.html) and
  [DeepTMHMM, Hallgren et al., 2022](https://doi.org/10.1101/2022.04.08.487609), when enabled.
