# Trimmings library data attribution

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

## Ligand Search archive

| Source | Credit and terms |
| --- | --- |
| [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) | wwPDB and original depositors; [CC0 1.0](https://www.wwpdb.org/about/usage-policies). CCD coordinates, connectivity and descriptors are selected and re-encoded as a search index. |
| [RCSB PDB](https://www.rcsb.org/pages/usage-policy) | RCSB PDB and wwPDB; occurrence counts are derived from the public PDB entry/component data. |
| [ChEBI](https://www.ebi.ac.uk/chebi/about) | EMBL-EBI ChEBI team; [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Trimmings selects metabolite records, joins them to CCD or PubChem identities, and derives filters from chemical class and role annotations. ChEBI IDs and the ontology version are retained. |
| [PubChem3D](https://pubchem.ncbi.nlm.nih.gov/docs/pubchem3d) | NCBI/NLM PubChem; [NCBI data-use policy](https://www.ncbi.nlm.nih.gov/home/about/policies/). The non-CCD coordinates are computed PubChem3D conformers, not depositor-supplied experimental coordinates. Trimmings verifies full InChIKeys, retains PubChem CIDs, and re-encodes coordinates and connectivity for fitting. |
| [GNPS2 ClassyFire2](https://apps.gnps2.org/classyfire) / ChemOnt | Credit GNPS2 for the computed compound-to-class assignments and [Djoumbou Feunang et al.](https://doi.org/10.1186/s13321-016-0174-y) for ChemOnt. The index includes class IDs, short labels, evidence and engine/rule versions, not classifier code, rules, ontology definitions or standalone taxonomy files. No upstream data or service is relicensed by Trimmings. |

This archive is a modified collection, not an upstream release. Keep this notice
with redistributed copies. Its manifest records source links and snapshot counts;
upstream identifiers remain in the indexes and are included in hit exports where available.

## Domain libraries

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
