# All coding types v2 — final evaluation

For research and educational use only; not for clinical diagnosis or treatment decisions.

Test opened once. Model changed after opening: False.

Routes chosen on validation: canonical_splice → lightgbm, frameshift → lightgbm, inframe_deletion → average, inframe_indel → lightgbm, inframe_insertion → average, missense → average, start_lost → lightgbm, stop_gained → lightgbm, stop_lost → lightgbm, synonymous → average

Temporal variants by gene set: test_gene: 20,856, train_gene: 88,360, validation_gene: 21,484; 13,209 share a training position.

Every table is per variant type; no number mixes types. Comparison columns are ROC-AUC of the joint LightGBM, the joint dual tower, and LightGBM trained without allele frequency or without ESM features.

## Gene-disjoint test (156,277 variants)

| Type | n | Pathogenic | ROC-AUC [95% CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE | Uncertain | LightGBM | Dual tower | LightGBM no-AF | LightGBM no-ESM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| canonical_splice | 6293 | 6202 | 0.9497 [0.9263–0.9733] | 0.9991 | 0.6296 | 0.7625 | 0.0083 | 0.0053 | 0.5% | 0.9497 | 0.8828 | 0.8395 | 0.9400 |
| frameshift | 14176 | 14038 | 0.9022 [0.8414–0.9592] | 0.9981 | 0.7846 | 0.8474 | 0.0035 | 0.0012 | 0.5% | 0.9022 | 0.8982 | 0.7209 | 0.9004 |
| inframe_deletion | 490 | 299 | 0.9902 [0.9823–0.9968] | 0.9924 | 0.9276 | 0.9659 | 0.0287 | 0.0232 | 10.8% | 0.9890 | 0.9875 | 0.9482 | 0.9877 |
| inframe_indel | 35 | 29 | too few of one class | | | | | | | | | | |
| inframe_insertion | 213 | 87 | 0.9851 [0.9674–0.9977] | 0.9782 | 0.7832 | 0.8714 | 0.0778 | 0.0958 | 22.1% | 0.9706 | 0.9876 | 0.8939 | 0.9711 |
| missense | 27500 | 8003 | 0.9902 [0.9886–0.9922] | 0.9763 | 0.9086 | 0.9537 | 0.0291 | 0.0021 | 11.2% | 0.9889 | 0.9883 | 0.9607 | 0.9842 |
| start_lost | 475 | 457 | too few of one class | | | | | | | | | | |
| stop_gained | 8531 | 8464 | 0.9719 [0.9537–0.9885] | 0.9997 | 0.7760 | 0.8355 | 0.0027 | 0.0005 | 0.3% | 0.9719 | 0.9643 | 0.8273 | 0.9735 |
| stop_lost | 55 | 42 | too few of one class | | | | | | | | | | |
| synonymous | 98509 | 81 | 0.8886 [0.8362–0.9356] | 0.2896 | 0.3142 | 0.5494 | 0.0007 | 0.0004 | 0.1% | 0.8806 | 0.8643 | 0.8310 | 0.8672 |

## Temporal (added to ClinVar after October 2025), all genes (130,700 variants)

| Type | n | Pathogenic | ROC-AUC [95% CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE | Uncertain | LightGBM | Dual tower | LightGBM no-AF | LightGBM no-ESM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| canonical_splice | 6563 | 6506 | 0.9087 [0.8309–0.9546] | 0.9986 | 0.4385 | 0.6484 | 0.0066 | 0.0040 | 0.4% | 0.9087 | 0.8082 | 0.7343 | 0.8992 |
| frameshift | 20585 | 20397 | 0.9340 [0.8931–0.9597] | 0.9992 | 0.6617 | 0.7418 | 0.0045 | 0.0030 | 0.4% | 0.9340 | 0.8471 | 0.7778 | 0.9277 |
| inframe_deletion | 534 | 290 | 0.9852 [0.9734–0.9937] | 0.9833 | 0.8967 | 0.9436 | 0.0436 | 0.0395 | 9.0% | 0.9825 | 0.9804 | 0.9496 | 0.9818 |
| inframe_indel | 59 | 47 | too few of one class | | | | | | | | | | |
| inframe_insertion | 250 | 68 | 0.9783 [0.9576–0.9947] | 0.9576 | 0.8457 | 0.9081 | 0.0463 | 0.0303 | 14.8% | 0.9707 | 0.9459 | 0.9181 | 0.9754 |
| missense | 18197 | 6116 | 0.9886 [0.9870–0.9900] | 0.9791 | 0.8874 | 0.9514 | 0.0397 | 0.0323 | 14.0% | 0.9879 | 0.9853 | 0.9727 | 0.9834 |
| start_lost | 251 | 244 | too few of one class | | | | | | | | | | |
| stop_gained | 9346 | 9278 | 0.8827 [0.8357–0.9178] | 0.9989 | 0.4375 | 0.6247 | 0.0054 | 0.0043 | 0.3% | 0.8827 | 0.7696 | 0.7533 | 0.8608 |
| stop_lost | 44 | 36 | too few of one class | | | | | | | | | | |
| synonymous | 74871 | 61 | 0.8582 [0.7964–0.9125] | 0.2184 | 0.2418 | 0.5410 | 0.0007 | 0.0004 | 0.1% | 0.8551 | 0.8342 | 0.8154 | 0.8580 |

## Temporal, test genes only (never seen in training or model choice) (20,856 variants)

| Type | n | Pathogenic | ROC-AUC [95% CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE | Uncertain | LightGBM | Dual tower | LightGBM no-AF | LightGBM no-ESM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| canonical_splice | 1130 | 1121 | too few of one class | | | | | | | | | | |
| frameshift | 2848 | 2824 | 0.8224 [0.6102–0.9893] | 0.9969 | 0.6776 | 0.7498 | 0.0041 | 0.0024 | 0.3% | 0.8224 | 0.8406 | 0.5546 | 0.8394 |
| inframe_deletion | 75 | 45 | 0.9830 [0.9476–1.0000] | 0.9864 | 0.8921 | 0.9333 | 0.0382 | 0.0458 | 8.0% | 0.9807 | 0.9833 | 0.9793 | 0.9778 |
| inframe_indel | 7 | 7 | too few of one class | | | | | | | | | | |
| inframe_insertion | 35 | 13 | too few of one class | | | | | | | | | | |
| missense | 2919 | 919 | 0.9882 [0.9848–0.9918] | 0.9728 | 0.8947 | 0.9570 | 0.0366 | 0.0289 | 13.5% | 0.9868 | 0.9852 | 0.9668 | 0.9816 |
| start_lost | 82 | 80 | too few of one class | | | | | | | | | | |
| stop_gained | 1405 | 1393 | too few of one class | | | | | | | | | | |
| stop_lost | 12 | 11 | too few of one class | | | | | | | | | | |
| synonymous | 12343 | 11 | too few of one class | | | | | | | | | | |

## Temporal, validation genes (used for model choice) (21,484 variants)

| Type | n | Pathogenic | ROC-AUC [95% CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE | Uncertain | LightGBM | Dual tower | LightGBM no-AF | LightGBM no-ESM |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| canonical_splice | 1152 | 1140 | too few of one class | | | | | | | | | | |
| frameshift | 3033 | 3002 | 0.9443 [0.8848–0.9931] | 0.9993 | 0.7091 | 0.8058 | 0.0038 | 0.0020 | 0.5% | 0.9443 | 0.9356 | 0.7305 | 0.9372 |
| inframe_deletion | 82 | 43 | 0.9911 [0.9756–1.0000] | 0.9919 | 0.9061 | 0.9487 | 0.0356 | 0.0342 | 11.0% | 0.9839 | 0.9911 | 0.9100 | 0.9797 |
| inframe_indel | 9 | 6 | too few of one class | | | | | | | | | | |
| inframe_insertion | 47 | 12 | too few of one class | | | | | | | | | | |
| missense | 2982 | 1033 | 0.9872 [0.9824–0.9913] | 0.9771 | 0.8909 | 0.9504 | 0.0395 | 0.0203 | 14.1% | 0.9852 | 0.9848 | 0.9633 | 0.9796 |
| start_lost | 63 | 62 | too few of one class | | | | | | | | | | |
| stop_gained | 1526 | 1511 | too few of one class | | | | | | | | | | |
| stop_lost | 11 | 9 | too few of one class | | | | | | | | | | |
| synonymous | 12579 | 11 | too few of one class | | | | | | | | | | |

Phase 1 missense model on the same test genes: ROC-AUC 0.9902.

Allele frequency: for loss-of-function types most ClinVar benign labels follow population frequency, so the no-AF column is the honest measure of what the model learned beyond it.
