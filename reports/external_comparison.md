# External comparison: AlphaMissense and REVEL

For research and educational use only; not for clinical diagnosis or treatment decisions.

AlphaMissense and REVEL scores were downloaded only to evaluate them next to our model on the same locked missense variants. They are never a model input, feature, label or filter, and our models did not change (model changed: False).

Read the numbers with these caveats:

- Our v2 model uses gnomAD allele frequency as an input. AlphaMissense and REVEL do not, but AlphaMissense learned its benign signal from human and primate population frequency, so it is frequency-informed through its training labels. The frequency-free row removes our direct use of it.
- Our test genes were never seen by our models. REVEL (trained on HGMD and rare ESP variants, from predictors trained on HGMD, SwissVar and ClinVar) and AlphaMissense (calibrated on ClinVar) may have seen these variants or genes; gene-disjointness does not apply to them.
- Both tools predate the temporal labels, but since the 2022 ClinGen guidance submitters use REVEL and AlphaMissense as ACMG PP3/BP4 evidence, so some newer ClinVar labels partly rest on these very scores. Both effects favour the external tools.
- When a tool has no score on our MANE transcript, its highest score on another transcript is used; the count is given under each set.

## Gene-disjoint test

26,522 of 27,500 variants scored by all four (7,802 pathogenic, 4,329 genes). Scores from a non-MANE transcript: AlphaMissense 7,136, REVEL 3,557.

| Scorer | ROC-AUC [95% gene CI] | PR-AUC [95% gene CI] |
|---|---|---|
| Ours: v2 (all features) | 0.9903 [0.9882–0.9919] | 0.9770 [0.9720–0.9808] |
| Ours: LightGBM without allele frequency | 0.9602 [0.9550–0.9639] | 0.9059 [0.8923–0.9147] |
| AlphaMissense (2023) | 0.9721 [0.9693–0.9748] | 0.9375 [0.9289–0.9448] |
| REVEL (2016) | 0.9765 [0.9740–0.9790] | 0.9464 [0.9375–0.9539] |

| Three-way call | Decided | Likely pathogenic correct | Likely benign correct | Pathogenic found | Benign found |
|---|---|---|---|---|---|
| Ours: v2 (0.10 / 0.90) | 88.8% | 97.9% | 99.1% | 78.5% | 91.5% |
| AlphaMissense (0.34 / 0.564) | 95.0% | 90.1% | 95.8% | 81.0% | 93.0% |

## Temporal (added to ClinVar after October 2025), all genes

17,671 of 18,197 variants scored by all four (6,010 pathogenic, 5,925 genes). Scores from a non-MANE transcript: AlphaMissense 4,865, REVEL 2,722.

| Scorer | ROC-AUC [95% gene CI] | PR-AUC [95% gene CI] |
|---|---|---|
| Ours: v2 (all features) | 0.9886 [0.9865–0.9903] | 0.9793 [0.9743–0.9833] |
| Ours: LightGBM without allele frequency | 0.9726 [0.9692–0.9757] | 0.9520 [0.9403–0.9608] |
| AlphaMissense (2023) | 0.9646 [0.9602–0.9697] | 0.9326 [0.9180–0.9440] |
| REVEL (2016) | 0.9731 [0.9685–0.9765] | 0.9521 [0.9396–0.9619] |

| Three-way call | Decided | Likely pathogenic correct | Likely benign correct | Pathogenic found | Benign found |
|---|---|---|---|---|---|
| Ours: v2 (0.10 / 0.90) | 85.9% | 96.4% | 99.3% | 86.4% | 83.4% |
| AlphaMissense (0.34 / 0.564) | 94.6% | 88.1% | 95.5% | 85.0% | 89.4% |

## Temporal, test genes only

2,800 of 2,919 variants scored by all four (888 pathogenic, 1,536 genes). Scores from a non-MANE transcript: AlphaMissense 738, REVEL 407.

| Scorer | ROC-AUC [95% gene CI] | PR-AUC [95% gene CI] |
|---|---|---|
| Ours: v2 (all features) | 0.9878 [0.9833–0.9916] | 0.9722 [0.9599–0.9809] |
| Ours: LightGBM without allele frequency | 0.9659 [0.9578–0.9722] | 0.9230 [0.9006–0.9405] |
| AlphaMissense (2023) | 0.9714 [0.9651–0.9775] | 0.9415 [0.9289–0.9534] |
| REVEL (2016) | 0.9740 [0.9679–0.9795] | 0.9436 [0.9289–0.9563] |

| Three-way call | Decided | Likely pathogenic correct | Likely benign correct | Pathogenic found | Benign found |
|---|---|---|---|---|---|
| Ours: v2 (0.10 / 0.90) | 86.2% | 95.8% | 99.6% | 82.9% | 85.8% |
| AlphaMissense (0.34 / 0.564) | 95.4% | 90.4% | 95.3% | 83.5% | 92.3% |
