# Missense model v1 — final evaluation

Variant type: **missense SNVs only** (GRCh38, MANE Select). For research and educational use only; not for clinical diagnosis or treatment decisions.

- Model: Missense ensemble (LightGBM + dual tower, Platt), ONNX (`missense_v1.onnx`)
- Test set opened in this report for the first time; model changed after testing: **no**
- Thresholds and calibration were fit on the validation set only (likely pathogenic ≥ 0.90, likely benign ≤ 0.10).
- Ablations: retrained with one seed from the frozen hyperparameters on training data only, before the locked sets were read; the full model is the exported ONNX (float32 trees); ensembles of ablated members are uncalibrated, so compare their ROC/PR, not Brier/ECE.
- Temporal sets: temporal_train_gene and temporal_validation_gene share genes with data the model saw (see shares_position_with_train per set); temporal_test_gene is the gene-disjoint temporal estimate.

## test

27501 variants, 4493 genes, 8003 pathogenic.

| Model | ROC-AUC [95% gene CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE |
|---|---|---|---|---|---|---|
| ensemble_calibrated | 0.9902 [0.9884–0.9921] | 0.9765 | 0.9086 | 0.9527 | 0.0296 | 0.0044 |
| lightgbm | 0.9895 [0.9878–0.9917] | 0.9743 | 0.8908 | 0.9356 | 0.0344 | 0.0276 |
| dual_tower | 0.9888 [0.9869–0.9908] | 0.9730 | 0.8952 | 0.9427 | 0.0328 | 0.0093 |
| no_frequency/lightgbm | 0.9616 [0.9580–0.9655] | 0.9090 | 0.7612 | 0.8739 | 0.0698 | 0.0230 |
| no_frequency/dual_tower | 0.9643 [0.9607–0.9676] | 0.9170 | 0.7856 | 0.8927 | 0.0657 | 0.0176 |
| no_frequency/ensemble_uncalibrated | 0.9687 [0.9657–0.9720] | 0.9263 | 0.7974 | 0.8967 | 0.0617 | 0.0210 |
| no_esm/lightgbm | 0.9849 [0.9827–0.9872] | 0.9630 | 0.8624 | 0.9195 | 0.0429 | 0.0314 |
| no_esm/dual_tower | 0.9832 [0.9808–0.9857] | 0.9590 | 0.8459 | 0.9066 | 0.0443 | 0.0330 |
| no_esm/ensemble_uncalibrated | 0.9851 [0.9829–0.9873] | 0.9641 | 0.8679 | 0.9228 | 0.0417 | 0.0299 |

| Decision | Share | Pathogenic fraction |
|---|---|---|
| likely_benign | 0.650 | 0.009 |
| uncertain | 0.122 | 0.508 |
| likely_pathogenic | 0.228 | 0.980 |

## temporal_all

18205 variants, 6116 genes, 6117 pathogenic.

| Model | ROC-AUC [95% gene CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE |
|---|---|---|---|---|---|---|
| ensemble_calibrated | 0.9887 [0.9869–0.9904] | 0.9797 | 0.8873 | 0.9508 | 0.0393 | 0.0333 |
| lightgbm | 0.9881 [0.9863–0.9896] | 0.9778 | 0.8932 | 0.9495 | 0.0368 | 0.0105 |
| dual_tower | 0.9869 [0.9849–0.9888] | 0.9762 | 0.8839 | 0.9459 | 0.0409 | 0.0325 |
| no_frequency/lightgbm | 0.9720 [0.9687–0.9753] | 0.9500 | 0.8183 | 0.9013 | 0.0600 | 0.0259 |
| no_frequency/dual_tower | 0.9756 [0.9720–0.9788] | 0.9579 | 0.8485 | 0.9220 | 0.0522 | 0.0101 |
| no_frequency/ensemble_uncalibrated | 0.9773 [0.9742–0.9801] | 0.9606 | 0.8511 | 0.9212 | 0.0513 | 0.0224 |
| no_esm/lightgbm | 0.9841 [0.9819–0.9860] | 0.9689 | 0.8740 | 0.9402 | 0.0434 | 0.0154 |
| no_esm/dual_tower | 0.9794 [0.9768–0.9820] | 0.9604 | 0.8537 | 0.9258 | 0.0519 | 0.0374 |
| no_esm/ensemble_uncalibrated | 0.9832 [0.9809–0.9854] | 0.9676 | 0.8729 | 0.9396 | 0.0453 | 0.0283 |

| Decision | Share | Pathogenic fraction |
|---|---|---|
| likely_benign | 0.550 | 0.006 |
| uncertain | 0.151 | 0.292 |
| likely_pathogenic | 0.299 | 0.967 |

## temporal_test_gene

2919 variants, 1593 genes, 919 pathogenic.

| Model | ROC-AUC [95% gene CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE |
|---|---|---|---|---|---|---|
| ensemble_calibrated | 0.9882 [0.9843–0.9914] | 0.9732 | 0.8911 | 0.9545 | 0.0364 | 0.0277 |
| lightgbm | 0.9863 [0.9821–0.9899] | 0.9659 | 0.8937 | 0.9498 | 0.0363 | 0.0126 |
| dual_tower | 0.9873 [0.9832–0.9907] | 0.9724 | 0.8856 | 0.9473 | 0.0390 | 0.0284 |
| no_frequency/lightgbm | 0.9673 [0.9603–0.9731] | 0.9260 | 0.7749 | 0.8782 | 0.0668 | 0.0360 |
| no_frequency/dual_tower | 0.9716 [0.9659–0.9768] | 0.9421 | 0.8192 | 0.9083 | 0.0587 | 0.0148 |
| no_frequency/ensemble_uncalibrated | 0.9737 [0.9679–0.9787] | 0.9434 | 0.8233 | 0.9083 | 0.0573 | 0.0261 |
| no_esm/lightgbm | 0.9814 [0.9769–0.9862] | 0.9510 | 0.8683 | 0.9386 | 0.0443 | 0.0133 |
| no_esm/dual_tower | 0.9804 [0.9754–0.9854] | 0.9521 | 0.8659 | 0.9333 | 0.0478 | 0.0419 |
| no_esm/ensemble_uncalibrated | 0.9818 [0.9773–0.9869] | 0.9541 | 0.8768 | 0.9425 | 0.0440 | 0.0289 |

| Decision | Share | Pathogenic fraction |
|---|---|---|
| likely_benign | 0.588 | 0.004 |
| uncertain | 0.144 | 0.375 |
| likely_pathogenic | 0.268 | 0.966 |

## temporal_train_gene

12304 variants, 2918 genes, 4165 pathogenic.

| Model | ROC-AUC [95% gene CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE |
|---|---|---|---|---|---|---|
| ensemble_calibrated | 0.9894 [0.9872–0.9912] | 0.9817 | 0.8869 | 0.9511 | 0.0398 | 0.0377 |
| lightgbm | 0.9890 [0.9868–0.9909] | 0.9807 | 0.8962 | 0.9518 | 0.0360 | 0.0169 |
| dual_tower | 0.9874 [0.9846–0.9894] | 0.9775 | 0.8848 | 0.9469 | 0.0412 | 0.0360 |
| no_frequency/lightgbm | 0.9746 [0.9697–0.9786] | 0.9568 | 0.8343 | 0.9098 | 0.0559 | 0.0233 |
| no_frequency/dual_tower | 0.9779 [0.9737–0.9812] | 0.9640 | 0.8598 | 0.9277 | 0.0482 | 0.0103 |
| no_frequency/ensemble_uncalibrated | 0.9794 [0.9755–0.9825] | 0.9662 | 0.8648 | 0.9278 | 0.0476 | 0.0223 |
| no_esm/lightgbm | 0.9858 [0.9826–0.9883] | 0.9738 | 0.8804 | 0.9445 | 0.0418 | 0.0229 |
| no_esm/dual_tower | 0.9792 [0.9754–0.9820] | 0.9618 | 0.8523 | 0.9258 | 0.0530 | 0.0399 |
| no_esm/ensemble_uncalibrated | 0.9842 [0.9806–0.9867] | 0.9711 | 0.8764 | 0.9424 | 0.0449 | 0.0333 |

| Decision | Share | Pathogenic fraction |
|---|---|---|
| likely_benign | 0.540 | 0.006 |
| uncertain | 0.154 | 0.257 |
| likely_pathogenic | 0.306 | 0.967 |

## temporal_validation_gene

2982 variants, 1605 genes, 1033 pathogenic.

| Model | ROC-AUC [95% gene CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE |
|---|---|---|---|---|---|---|
| ensemble_calibrated | 0.9867 [0.9816–0.9907] | 0.9772 | 0.8857 | 0.9463 | 0.0400 | 0.0206 |
| lightgbm | 0.9863 [0.9817–0.9901] | 0.9748 | 0.8806 | 0.9398 | 0.0404 | 0.0102 |
| dual_tower | 0.9852 [0.9794–0.9895] | 0.9746 | 0.8782 | 0.9404 | 0.0418 | 0.0234 |
| no_frequency/lightgbm | 0.9646 [0.9555–0.9711] | 0.9370 | 0.7925 | 0.8880 | 0.0701 | 0.0361 |
| no_frequency/dual_tower | 0.9679 [0.9596–0.9746] | 0.9413 | 0.8294 | 0.9115 | 0.0623 | 0.0162 |
| no_frequency/ensemble_uncalibrated | 0.9704 [0.9623–0.9764] | 0.9477 | 0.8210 | 0.9064 | 0.0607 | 0.0249 |
| no_esm/lightgbm | 0.9801 [0.9744–0.9852] | 0.9596 | 0.8535 | 0.9247 | 0.0493 | 0.0148 |
| no_esm/dual_tower | 0.9793 [0.9731–0.9844] | 0.9609 | 0.8478 | 0.9190 | 0.0514 | 0.0314 |
| no_esm/ensemble_uncalibrated | 0.9804 [0.9746–0.9857] | 0.9618 | 0.8550 | 0.9254 | 0.0483 | 0.0256 |

| Decision | Share | Pathogenic fraction |
|---|---|---|
| likely_benign | 0.556 | 0.008 |
| uncertain | 0.145 | 0.369 |
| likely_pathogenic | 0.298 | 0.966 |

## Figures

- `figures/missense_v1_curves_test.png`, `figures/missense_v1_curves_temporal_test_gene.png`, `figures/missense_v1_reliability.png`
## Interpretation

- **Gene-disjoint test (headline for missense):** ROC-AUC 0.9902 [0.9884–0.9921], PR-AUC 0.9765, ECE 0.0044;
  12.2% of variants are left Uncertain, Likely pathogenic calls are 98.0% pathogenic, Likely benign calls 0.9%.
  Test metrics match validation (0.9901), so validation-based selection did not overfit noticeably.
- **Temporal (variants first added to ClinVar after October 2025):** discrimination holds (ROC-AUC 0.9887; 0.9882 on
  the gene-disjoint `temporal_test_gene` subset), but calibration drifts (ECE 0.033): in the 0.3–0.8 range the
  observed pathogenic fraction is clearly below the predicted probability. Likely pathogenic precision drops from
  98.0% to 96.6–96.7%, still above the 90% target. The calibrator is not refit on temporal data, because the model
  must not change after testing.
- **Allele frequency ablation:** without any frequency feature, test ROC-AUC falls from 0.9902 to 0.9687
  (ensemble) and PR-AUC from 0.9765 to 0.9263. Frequency carries real signal, and part of it may be
  circular (some ClinVar benign labels rely on frequency).
- **ESM ablation:** without ESM-2 inputs, test ROC-AUC 0.9851 vs 0.9902 and PR-AUC 0.9641 vs 0.9765. The protein
  language model helps most when frequency is unavailable (dual tower no-AF 0.9643 vs LightGBM no-AF 0.9616).
- **LightGBM vs dual tower:** close on discrimination (test 0.9895 vs 0.9888). Their calibrated average has the best
  ROC-AUC and PR-AUC on every set. On test the dual tower is better calibrated than LightGBM (ECE 0.0093 vs 0.0276);
  on the temporal sets raw LightGBM is better calibrated than the final model (ECE 0.0105 vs 0.0333).
