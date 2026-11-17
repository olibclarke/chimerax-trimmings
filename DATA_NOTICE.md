# Domain-library data attribution

Trimmings software is MIT licensed. Database contents retain their upstream
terms; the MIT license does not relicense third-party data.

| Source | Credit and terms |
| --- | --- |
| [ECOD](http://prodata.swmed.edu/ecod/) | ECOD team, UT Southwestern/UCF; [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). [Reference](https://doi.org/10.1093/nar/gkag1132). |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/) | Google DeepMind and EMBL-EBI; [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). [Reference](https://doi.org/10.1093/nar/gkab1061). |
| [UniProt](https://www.uniprot.org/) | UniProt Consortium; [CC BY 4.0](https://www.uniprot.org/help/license). |
| [Pfam/InterPro](https://www.ebi.ac.uk/interpro/) | EMBL-EBI and contributing resource teams; [CC0 1.0 and attribution guidance](https://interpro-documentation.readthedocs.io/en/latest/license.html). |
| [PDB](https://www.wwpdb.org/) | wwPDB and original structure depositors; [CC0 1.0](https://www.wwpdb.org/about/usage-policies). Retain structure IDs and primary-paper credits. |
| [PDBe/SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) | PDBe/UniProt; [PDBe CC0 policy](https://www.ebi.ac.uk/pdbe/about/public-data-access-statement). [Reference](https://doi.org/10.1093/nar/gky1114). |
| [NCBI CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml), where included | NCBI and original member-database providers; [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/). Imported material retains its original provider terms. |

Changes by Trimmings include domain selection/cropping, C-alpha coordinate
conversion/translation, sequence-profile validation/reannotation, metadata joins
and Foldseek database encoding. This is a modified collection, not an official
upstream release; no endorsement is implied. Source accessions, domain UIDs,
versions, provenance and coverage limits are recorded in the supplied manifest,
reference metadata and coverage report. Sources apply to the records actually
included in each bundle. Preserve these notices when redistributing it.

Built and searched with [Foldseek](https://doi.org/10.1038/s41587-023-01773-0),
[HMMER](https://doi.org/10.1371/journal.pcbi.1002195) and
[Gemmi](https://doi.org/10.21105/joss.04200). No program binaries, PROSITE
definitions or DeepTMHMM weights are included. See `THIRD_PARTY.md` for software
and optional-service terms.
