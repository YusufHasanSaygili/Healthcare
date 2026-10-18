# Model card — all coding types v2

For research and educational use only; not for clinical diagnosis or treatment decisions.

## Model details

- **Name:** `all_types_v2.onnx` (13.9 MB, one ONNX graph, CPU inference with ONNX Runtime)
- **Built:** 2026-10-18, as a personal portfolio project
- **Task:** binary classification of germline coding variants in GRCh38, pathogenic (1) vs benign (0)
- **Inputs:** 55 `f_` features computed in this repository, the ESM-2 blocks for the variant's protein position, and
  the variant type. No precomputed pathogenicity score of any kind (CADD, REVEL, AlphaMissense, SpliceAI, SIFT,
  PolyPhen-2, dbNSFP, …) is an input, label, filter or benchmark.
- **Components:**
  - LightGBM trained jointly on all coding types (40-trial sweep, selected by the mean per-type validation ROC-AUC)
  - Dual tower network: protein tower (ESM-2 650M embedding window and log-likelihood features) and measurement tower,
    fused by an MLP head; variant type one-hot; the protein tower is masked for rows without a protein position and
    randomly during training; the seed-0 model is used, 3 seeds were trained to measure variance
  - Per-type route chosen on validation: average of both models for missense, in-frame deletion, in-frame insertion
    and synonymous; LightGBM alone for stop gained, frameshift, canonical splice and the small types
  - Per-type Platt calibration fit on validation (small types pooled)
- **Output:** calibrated probability of pathogenicity and a three-way decision: Likely pathogenic (≥ 0.90),
  Uncertain, Likely benign (≤ 0.10). Explanations report the change in calibrated log-odds when one feature group
  is replaced by values from 32 seeded training variants.

## Intended use

- Exploring how far public raw data and a protein language model go without other predictors' scores
- Teaching and portfolio demonstration, run locally on `127.0.0.1`

**Not intended for** clinical diagnosis, treatment decisions, carrier screening, or any decision about a person.
The output is not an ACMG/AMP classification.

## Training data

| Source | Version | Role |
|---|---|---|
| ClinVar | 2025-10-06 | Labels: Pathogenic / Likely pathogenic → 1, Benign / Likely benign → 0; germline; ≥ 1 review star; VUS and conflicting excluded |
| GRCh38 (Ensembl 115), MANE Select v1.5 | — | Our annotator: variant type, codon and amino acid change, NMD, splice side |
| gnomAD v4.1 exomes | — | Allele frequency and counts; site coverage separates "not observed" from "not measured" |
| gnomAD v4.1 constraint | — | LOEUF, missense Z |
| UCSC phyloP100way, phastCons100way | hg38 | Conservation at the site (and across the span for indels) |
| UniProtKB/Swiss-Prot | 2026_03 | Curated regions mapped to MANE proteins with identical sequence (90.4% of proteins) |
| ESM-2 `esm2_t33_650M_UR50D` | pinned revision | Wild-type marginals, log-likelihood ratio, entropy, embeddings |

Training rows (missense + other coding types): 589,925 in train genes, 126,322 in validation genes. Synonymous
variants are not subsampled. Variants scored in more than one MANE gene are dropped.

## Evaluation design

- **Gene-disjoint split:** all 19,229 MANE genes assigned to train / validation / test (≈ 70/15/15, fixed seed);
  genes with overlapping coding sequence or a shared labeled position form one group. Split file pinned by hash.
- **Locked test:** test and temporal rows are stored separately and readable only by final-evaluation code. The v2
  test set was opened once; the models did not change afterwards (hashes in the report).
- **Temporal test:** variants first added to ClinVar between 2025-10-06 and 2026-10-04 and absent from training;
  reported for all genes and separately for test genes only.
- **Selection:** routes, calibration, thresholds and ablation models were all fixed on validation before the test
  was opened.
- Metrics per variant type with 95% gene-level bootstrap confidence intervals.

## Results (gene-disjoint test)

| Type | n | Minority class | ROC-AUC [95% CI] | PR-AUC | MCC | ECE | Uncertain |
|---|---|---|---|---|---|---|---|
| Missense | 27,500 | 8,003 P | 0.990 [0.989–0.992] | 0.976 | 0.909 | 0.002 | 11.2% |
| Stop gained | 8,531 | 67 B | 0.972 [0.954–0.989] | 0.9997 | 0.776 | 0.001 | 0.3% |
| Canonical splice | 6,293 | 91 B | 0.950 [0.926–0.973] | 0.999 | 0.630 | 0.005 | 0.5% |
| Frameshift | 14,176 | 138 B | 0.902 [0.841–0.959] | 0.998 | 0.785 | 0.001 | 0.5% |
| In-frame deletion | 490 | 191 B | 0.990 [0.982–0.997] | 0.992 | 0.928 | 0.023 | 10.8% |
| In-frame insertion | 213 | 87 P | 0.985 [0.967–0.998] | 0.978 | 0.783 | 0.096 | 22.1% |
| Synonymous | 98,509 | 81 P | 0.889 [0.836–0.936] | 0.290 | 0.314 | 0.000 | 0.1% |

Start lost (475), stop lost (55) and in-frame delins (35) have too few of one class to score.

**Temporal test, ROC-AUC:** missense 0.989 (test genes only 0.988), frameshift 0.934 (0.822), stop gained 0.883,
canonical splice 0.909, in-frame deletion 0.985 (0.983), synonymous 0.858. Missense ECE rises to 0.032.

## Ablations (gene-disjoint test, ROC-AUC)

| Type | v2 | LightGBM | Dual tower | LightGBM without allele frequency | LightGBM without ESM |
|---|---|---|---|---|---|
| Missense | 0.990 | 0.989 | 0.988 | 0.961 | 0.984 |
| Stop gained | 0.972 | 0.972 | 0.964 | 0.827 | 0.974 |
| Canonical splice | 0.950 | 0.950 | 0.883 | 0.840 | 0.940 |
| Frameshift | 0.902 | 0.902 | 0.898 | 0.721 | 0.900 |
| In-frame deletion | 0.990 | 0.989 | 0.988 | 0.948 | 0.988 |
| In-frame insertion | 0.985 | 0.971 | 0.988 | 0.894 | 0.971 |
| Synonymous | 0.889 | 0.881 | 0.864 | 0.831 | 0.867 |

The Phase 1 missense-only model scored 0.990 on the same test genes ([`missense_v1.md`](missense_v1.md)).

## Known limitations and risks

- **Allele frequency circularity.** Many ClinVar benign labels, especially for loss-of-function types, were
  assigned because the variant is common. The no-frequency column is the better estimate of what the model learned
  beyond that rule; for frameshift it falls to 0.72.
- **Label noise.** ClinVar contains misclassified variants; the model inherits them.
- **Ancestry bias.** gnomAD and ClinVar over-represent European ancestry, so frequency evidence is weaker elsewhere.
- **Small classes.** Benign loss-of-function and pathogenic synonymous variants are rare; their intervals are wide
  and three types could not be scored.
- **Calibration drift** on newer ClinVar submissions (overconfident between 0.3 and 0.8 for missense).
- **Scope.** GRCh38, MANE Select only, coding SNVs and small indels. Deep intronic, non-canonical splice, UTR,
  regulatory, mitochondrial, structural and somatic variants are out of scope and are rejected by the service.
- **Single transcript.** Effects on non-MANE isoforms are ignored.

## Reproducibility

Source versions and checksums: [`config/sources.toml`](../config/sources.toml). Steps: [README](../README.md#reproduce-from-scratch).
Detailed numbers: [`all_types_v2.md`](all_types_v2.md), [`all_types_v2.json`](all_types_v2.json),
[`joint_lightgbm.json`](joint_lightgbm.json), [`joint_dual_tower.json`](joint_dual_tower.json),
[`onnx_export_v2.json`](onnx_export_v2.json).
