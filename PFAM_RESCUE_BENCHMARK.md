# Pfam sensitivity experiment — 14 November 2026

The useful improvements in these controls were **newer Pfam profiles and padded
region searches**. Disabling HMMER's screening filters with `--max` added no
accepted matches and took longer than ordinary targeted region searches.
The subsequent [local-Pfam integration](PFAM_LIBRARY.md#use-the-profiles-in-chimerax)
uses ordinary padded retries, with an online fallback for systems without a local
installation. The timings below describe the original controlled experiment.

## Results

| Case | Current online service | Standard local Pfam 39.0 | Additional accepted matches from padded retries | Additional benefit from `--max` |
| --- | --- | --- | --- | --- |
| A0AAE4MDN0, 643 residues | No reported matches | PF39420 / DUF14543, 304–452 | None | None |
| A0A1V4UJD3, 548 residues | TrkA_N 65–155; RyR 460–528 | Same matches plus PF39420 / DUF14543, 227–366 | None | None |
| Large RyR subunit, 4,249 modeled sequence positions | Previously recorded 17 curated assignments | 17 assignments | EF-hand_7 / PF13499, author residues **4075–4132** | None beyond ordinary retries |

Two fresh API searches confirmed that the online service advertises **Pfam 37.2**.
The local sequence-profile database is **Pfam 39.0**, with all **53,275** models.
PF39420 is a sequence-family assignment in the second RCK-like structural region
of each small protein; its annotated function remains unknown. It is not a
functional annotation of an RCK domain.

The first RCK-like region and RyR-like region of A0AAE4MDN0 remain below Pfam
gathering thresholds. The weaker RyR occurrence in A0A1V4UJD3 is not rescued.
A separate control targeting Ion_trans / PF00520 on A0AAE4MDN0 residues 1–120
also produced no accepted match, including with `--max`. These negative results
do not establish that the proposed folds are absent.

The large RyR EF-hand occupies modeled sequence positions **3639–3696** in this
structure, which map to author residues **4075–4132**. Its profile coverage is
**91%**. The domain score is **23.5 bits** in the full chain, below the **30-bit**
gathering cutoff, and **31.9 bits** in the automatically proposed 88-residue
window (sequence positions 3622–3709). The cropped independent E-value is
**3.8 × 10⁻⁷**, calculated using the complete profile count. The aligned interval
and composition-bias score are the same: the gain comes from changing query
context, rather than removing screening filters or lowering the acceptance score.
An independent Foldseek match already supports an EF-hand in this region.

## Runtime

Median wall times from three sequential repetitions, HMMER 3.4, four CPU workers,
with the filesystem cache warm:

| Case | Complete Pfam scan | Ordinary targeted region retries | Targeted region retries with `--max` |
| --- | ---: | ---: | ---: |
| A0AAE4MDN0 | 0.44 s | 0.010 s | 0.027 s |
| A0A1V4UJD3 | 0.43 s | 0.015 s | 0.128 s |
| Large RyR subunit | 0.50 s | 0.033 s | 0.237 s |

Candidate-profile extraction and indexing add approximately **0.035, 0.105,
and 0.065 seconds**, respectively. These measurements exclude ChimeraX rendering,
network requests, and the prior Foldseek searches supplying candidate families
and structural intervals. They describe sequence scanning, not a complete
Foldseek annotation run. The FKBP and calmodulin representatives were also checked;
`--max` added no accepted matches. Only one representative of each repeated
sequence in the twelve-chain complex was tested.

## Method and limitations

- The baseline uses the complete local database with `hmmscan --cut_ga`.
- Candidate panels contain explicit Pfam mappings from qualifying, previously
  computed Foldseek matches, plus all families accepted by the baseline scan.
  The panels contain 8, 33, and 19 models for the two small proteins and RyR.
- Retry windows cover unassigned stretches of at least 20 residues and structural
  intervals with less than 80% coverage by baseline sequence annotations, padded
  by 12 residues at each end. The experiment uses 4, 6, and 22 windows, respectively.
- Standard versus `--max` searches use the same panels and query windows. Both
  sequence and domain gathering cutoffs remain fixed. Composition-bias score
  correction remains enabled. Reported accepted coordinates and domain scores
  were identical between the corresponding standard and `--max` runs.
- `-Z 53275` retains the complete profile count for independent E-values when
  searching a subset. Conditional domain E-values are recorded separately.
- Local tables precede Pfam clan competition: overlapping family matches are not
  independent domains. Calmodulin's competing EF-hand profiles are therefore not
  counted as additional domain discoveries.
- Below-cutoff reporting was examined separately and never counted as an accepted
  rescue. Repeated crops add searches and change score context; fixed gathering
  thresholds do not guarantee an unchanged false-positive rate. Structural
  corroboration and validation on a broader control set remain useful.
- This is a small control set and a targeted strategy. A family omitted by the
  initial candidate selection cannot be recovered by its targeted retries.

The next practical improvement is a local Pfam 39.0 option with ordinary padded
retries where structural evidence suggests a missed domain. These controls do
not justify adding a user-facing `--max` sensitivity setting.

## Reproduce

The benchmark accepts general case definitions; it contains no RyR-specific
search logic. Case definitions, extracted candidate profiles, query FASTAs, all
HMMER domain tables, command logs, and structured results are retained under
`.domain_data/diagnostics/pfam-rescue/`. Repetitions are in `repeat-2/` and
`repeat-3/`; `service_comparison.json` records the fresh API comparison and
database versions. Input definitions contain sequences, original residue maps,
and the candidate families and intervals used for the experiment.

```bash
python3 benchmark_pfam_rescue.py \
  --input .domain_data/diagnostics/pfam-rescue/inputs.json \
  --profiles .domain_data/pfam-complete/profiles/Pfam39.0/Pfam-A.hmm \
  --metadata .domain_data/pfam-audit/v39.0/Pfam-A.hmm.dat.gz \
  --models 53275 --threads 4 \
  --out .domain_data/diagnostics/pfam-rescue/recheck
```
