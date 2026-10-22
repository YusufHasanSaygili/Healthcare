# External comparison: AlphaMissense and REVEL

For research and educational use only; not for clinical diagnosis or treatment decisions.

AlphaMissense and REVEL scores were downloaded only to evaluate them next to our model on the same locked missense variants. They are never a model input, feature, label or filter, and our models did not change (model changed: False).

Read the numbers with these caveats:

- Our v2 model uses gnomAD allele frequency as an input. AlphaMissense and REVEL do not, but AlphaMissense learned its benign signal from human and primate population frequency, so it is frequency-informed through its training labels. The frequency-free row removes our direct use of it.
- Our test genes were never seen by our models. REVEL (trained on HGMD and rare ESP variants, from predictors trained on HGMD, SwissVar and ClinVar) and AlphaMissense (calibrated on ClinVar) may have seen these variants or genes; gene-disjointness does not apply to them.
- Both tools predate the temporal labels, but since the 2022 ClinGen guidance submitters use REVEL and AlphaMissense as ACMG PP3/BP4 evidence, so some newer ClinVar labels partly rest on these very scores. Both effects favour the external tools.
- When a tool has no score on our MANE transcript, its highest score on another transcript is used; the count is given under each set.
- dbNSFP tables use each tool's rankscore (oriented so higher means more damaging, taken from its most damaging transcript). Several tools there, such as ClinPred, MetaRNN and BayesDel with AF, use allele frequency and were trained on ClinVar labels, so on the gene-disjoint test set they may have seen our test variants. The temporal sets guard against that only if a tool's training data predates October 2025, which we could not verify tool by tool.

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

### Gene-disjoint test: every predictor in dbNSFP 5.4a

Each tool is scored on the variants it covers, and our model on exactly those variants. The difference is ours minus the tool, with a paired 95% gene-bootstrap interval.

| Tool | Variants | Tool ROC-AUC | Ours v2 | Difference [95% CI] | Ours without AF |
|---|---|---|---|---|---|
| MetaRNN | 27,236 | 0.9927 | 0.9902 | -0.0025 [-0.0042, -0.0013] | 0.9607 |
| ClinPred | 27,144 | 0.9919 | 0.9902 | -0.0017 [-0.0030, -0.0005] | 0.9605 |
| BayesDel (with AF) | 27,317 | 0.9901 | 0.9902 | +0.0001 [-0.0012, +0.0018] | 0.9606 |
| MutationTaster | 25,666 | 0.9867 | 0.9909 | +0.0042 [+0.0021, +0.0062] | 0.9616 |
| PHACTboost | 26,220 | 0.9845 | 0.9902 | +0.0057 [+0.0038, +0.0074] | 0.9600 |
| MutScore | 27,121 | 0.9791 | 0.9903 | +0.0111 [+0.0091, +0.0133] | 0.9607 |
| REVEL | 26,190 | 0.9760 | 0.9903 | +0.0143 [+0.0117, +0.0168] | 0.9608 |
| BayesDel (no AF) | 27,317 | 0.9756 | 0.9902 | +0.0146 [+0.0120, +0.0177] | 0.9606 |
| VEST4 | 26,525 | 0.9725 | 0.9903 | +0.0178 [+0.0155, +0.0201] | 0.9614 |
| VARITY_R | 25,698 | 0.9725 | 0.9903 | +0.0178 [+0.0152, +0.0211] | 0.9597 |
| AlphaMissense | 27,379 | 0.9712 | 0.9902 | +0.0190 [+0.0163, +0.0223] | 0.9605 |
| VARITY_ER | 25,698 | 0.9679 | 0.9903 | +0.0224 [+0.0196, +0.0250] | 0.9597 |
| gMVP | 25,888 | 0.9672 | 0.9904 | +0.0232 [+0.0205, +0.0257] | 0.9610 |
| MutPred2 | 27,473 | 0.9633 | 0.9902 | +0.0269 [+0.0240, +0.0299] | 0.9607 |
| Eigen | 26,052 | 0.9543 | 0.9902 | +0.0359 [+0.0329, +0.0388] | 0.9609 |
| CADD | 27,500 | 0.9526 | 0.9902 | +0.0376 [+0.0346, +0.0411] | 0.9607 |
| DEOGEN2 | 26,503 | 0.9514 | 0.9902 | +0.0388 [+0.0327, +0.0450] | 0.9599 |
| MutFormer | 27,115 | 0.9509 | 0.9903 | +0.0394 [+0.0340, +0.0444] | 0.9608 |
| M-CAP | 22,881 | 0.9475 | 0.9879 | +0.0404 [+0.0348, +0.0448] | 0.9661 |
| MVP | 25,561 | 0.9469 | 0.9897 | +0.0428 [+0.0374, +0.0480] | 0.9616 |
| MisFit | 26,918 | 0.9467 | 0.9902 | +0.0435 [+0.0387, +0.0499] | 0.9604 |
| MetaLR | 27,148 | 0.9443 | 0.9902 | +0.0459 [+0.0408, +0.0529] | 0.9605 |
| MetaSVM | 27,148 | 0.9431 | 0.9902 | +0.0471 [+0.0415, +0.0553] | 0.9605 |
| PolyPhen-2 HVAR | 25,985 | 0.9424 | 0.9907 | +0.0483 [+0.0439, +0.0520] | 0.9603 |
| FATHMM-XF | 25,976 | 0.9415 | 0.9902 | +0.0487 [+0.0433, +0.0529] | 0.9611 |
| GPN-MSA | 27,500 | 0.9393 | 0.9902 | +0.0510 [+0.0461, +0.0554] | 0.9607 |
| popEVE | 26,051 | 0.9389 | 0.9912 | +0.0523 [+0.0462, +0.0592] | 0.9614 |
| ESM1b | 27,499 | 0.9372 | 0.9902 | +0.0530 [+0.0459, +0.0603] | 0.9607 |
| PROVEAN | 25,892 | 0.9356 | 0.9903 | +0.0547 [+0.0506, +0.0588] | 0.9609 |
| SIFT | 25,580 | 0.9311 | 0.9904 | +0.0593 [+0.0547, +0.0635] | 0.9611 |
| PolyPhen-2 HDIV | 25,985 | 0.9297 | 0.9907 | +0.0609 [+0.0551, +0.0659] | 0.9603 |
| MutationAssessor | 24,623 | 0.9292 | 0.9907 | +0.0614 [+0.0548, +0.0686] | 0.9601 |
| SIFT4G | 26,551 | 0.9222 | 0.9903 | +0.0681 [+0.0625, +0.0739] | 0.9615 |
| LIST-S2 | 26,997 | 0.9167 | 0.9905 | +0.0738 [+0.0672, +0.0812] | 0.9605 |
| DANN | 27,366 | 0.9032 | 0.9902 | +0.0870 [+0.0818, +0.0928] | 0.9606 |
| PrimateAI | 26,744 | 0.8985 | 0.9903 | +0.0917 [+0.0838, +0.1009] | 0.9604 |
| MPC | 24,231 | 0.8869 | 0.9905 | +0.1036 [+0.0915, +0.1149] | 0.9608 |

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

### Temporal (added to ClinVar after October 2025), all genes: every predictor in dbNSFP 5.4a

Each tool is scored on the variants it covers, and our model on exactly those variants. The difference is ours minus the tool, with a paired 95% gene-bootstrap interval.

| Tool | Variants | Tool ROC-AUC | Ours v2 | Difference [95% CI] | Ours without AF |
|---|---|---|---|---|---|
| MetaRNN | 18,004 | 0.9847 | 0.9886 | +0.0039 [+0.0021, +0.0057] | 0.9727 |
| ClinPred | 18,009 | 0.9810 | 0.9886 | +0.0076 [+0.0055, +0.0094] | 0.9726 |
| BayesDel (with AF) | 18,101 | 0.9807 | 0.9886 | +0.0079 [+0.0053, +0.0107] | 0.9726 |
| MutScore | 17,923 | 0.9792 | 0.9886 | +0.0095 [+0.0069, +0.0123] | 0.9727 |
| MutationTaster | 16,845 | 0.9771 | 0.9891 | +0.0120 [+0.0093, +0.0147] | 0.9733 |
| PHACTboost | 16,852 | 0.9748 | 0.9886 | +0.0139 [+0.0113, +0.0170] | 0.9724 |
| REVEL | 17,370 | 0.9734 | 0.9889 | +0.0154 [+0.0120, +0.0188] | 0.9730 |
| gMVP | 16,847 | 0.9724 | 0.9885 | +0.0162 [+0.0138, +0.0194] | 0.9727 |
| BayesDel (no AF) | 18,101 | 0.9723 | 0.9886 | +0.0163 [+0.0134, +0.0195] | 0.9726 |
| VARITY_R | 16,272 | 0.9712 | 0.9881 | +0.0169 [+0.0145, +0.0202] | 0.9715 |
| VARITY_ER | 16,272 | 0.9693 | 0.9881 | +0.0188 [+0.0163, +0.0217] | 0.9715 |
| MutPred2 | 18,189 | 0.9677 | 0.9886 | +0.0209 [+0.0183, +0.0236] | 0.9727 |
| VEST4 | 17,648 | 0.9663 | 0.9887 | +0.0223 [+0.0194, +0.0254] | 0.9728 |
| AlphaMissense | 17,855 | 0.9660 | 0.9887 | +0.0227 [+0.0199, +0.0263] | 0.9727 |
| MVP | 17,549 | 0.9541 | 0.9884 | +0.0343 [+0.0282, +0.0411] | 0.9726 |
| DEOGEN2 | 17,297 | 0.9537 | 0.9888 | +0.0351 [+0.0308, +0.0404] | 0.9729 |
| MisFit | 17,667 | 0.9499 | 0.9886 | +0.0386 [+0.0327, +0.0451] | 0.9726 |
| Eigen | 16,767 | 0.9481 | 0.9883 | +0.0402 [+0.0327, +0.0502] | 0.9727 |
| MetaLR | 18,013 | 0.9426 | 0.9886 | +0.0460 [+0.0402, +0.0544] | 0.9726 |
| MetaSVM | 18,013 | 0.9423 | 0.9886 | +0.0463 [+0.0394, +0.0542] | 0.9726 |
| ESM1b | 18,193 | 0.9416 | 0.9886 | +0.0470 [+0.0414, +0.0530] | 0.9727 |
| M-CAP | 17,389 | 0.9412 | 0.9881 | +0.0469 [+0.0409, +0.0545] | 0.9731 |
| PolyPhen-2 HVAR | 16,847 | 0.9370 | 0.9887 | +0.0517 [+0.0440, +0.0595] | 0.9729 |
| FATHMM-XF | 16,669 | 0.9361 | 0.9883 | +0.0522 [+0.0428, +0.0639] | 0.9726 |
| CADD | 18,197 | 0.9358 | 0.9886 | +0.0528 [+0.0472, +0.0578] | 0.9727 |
| PROVEAN | 17,202 | 0.9352 | 0.9891 | +0.0539 [+0.0465, +0.0611] | 0.9732 |
| MutationAssessor | 16,221 | 0.9347 | 0.9887 | +0.0540 [+0.0456, +0.0621] | 0.9729 |
| popEVE | 14,457 | 0.9335 | 0.9894 | +0.0558 [+0.0464, +0.0640] | 0.9733 |
| GPN-MSA | 18,197 | 0.9334 | 0.9886 | +0.0552 [+0.0477, +0.0634] | 0.9727 |
| MutFormer | 17,965 | 0.9269 | 0.9886 | +0.0617 [+0.0459, +0.0861] | 0.9727 |
| SIFT4G | 16,936 | 0.9200 | 0.9887 | +0.0687 [+0.0629, +0.0753] | 0.9727 |
| SIFT | 17,079 | 0.9187 | 0.9891 | +0.0704 [+0.0631, +0.0769] | 0.9732 |
| PolyPhen-2 HDIV | 16,847 | 0.9179 | 0.9887 | +0.0708 [+0.0616, +0.0792] | 0.9729 |
| LIST-S2 | 17,698 | 0.9000 | 0.9884 | +0.0884 [+0.0809, +0.0982] | 0.9720 |
| PrimateAI | 17,666 | 0.8894 | 0.9885 | +0.0991 [+0.0864, +0.1095] | 0.9725 |
| DANN | 18,124 | 0.8705 | 0.9886 | +0.1181 [+0.1086, +0.1258] | 0.9726 |
| MPC | 15,455 | 0.8536 | 0.9893 | +0.1357 [+0.1117, +0.1593] | 0.9736 |

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

### Temporal, test genes only: every predictor in dbNSFP 5.4a

Each tool is scored on the variants it covers, and our model on exactly those variants. The difference is ours minus the tool, with a paired 95% gene-bootstrap interval.

| Tool | Variants | Tool ROC-AUC | Ours v2 | Difference [95% CI] | Ours without AF |
|---|---|---|---|---|---|
| MetaRNN | 2,887 | 0.9887 | 0.9880 | -0.0006 [-0.0038, +0.0027] | 0.9663 |
| ClinPred | 2,880 | 0.9885 | 0.9880 | -0.0004 [-0.0042, +0.0024] | 0.9665 |
| BayesDel (with AF) | 2,895 | 0.9844 | 0.9881 | +0.0037 [-0.0002, +0.0076] | 0.9668 |
| MutScore | 2,877 | 0.9825 | 0.9880 | +0.0055 [+0.0017, +0.0095] | 0.9667 |
| PHACTboost | 2,774 | 0.9780 | 0.9876 | +0.0096 [+0.0051, +0.0142] | 0.9657 |
| MutationTaster | 2,669 | 0.9773 | 0.9892 | +0.0120 [+0.0066, +0.0171] | 0.9679 |
| VARITY_R | 2,715 | 0.9759 | 0.9874 | +0.0115 [+0.0077, +0.0152] | 0.9655 |
| REVEL | 2,768 | 0.9732 | 0.9882 | +0.0150 [+0.0092, +0.0215] | 0.9668 |
| BayesDel (no AF) | 2,895 | 0.9732 | 0.9881 | +0.0149 [+0.0094, +0.0199] | 0.9668 |
| VARITY_ER | 2,715 | 0.9732 | 0.9874 | +0.0142 [+0.0099, +0.0186] | 0.9655 |
| VEST4 | 2,813 | 0.9731 | 0.9886 | +0.0155 [+0.0111, +0.0199] | 0.9676 |
| AlphaMissense | 2,910 | 0.9716 | 0.9882 | +0.0165 [+0.0117, +0.0216] | 0.9666 |
| MutPred2 | 2,918 | 0.9687 | 0.9881 | +0.0195 [+0.0145, +0.0251] | 0.9668 |
| gMVP | 2,748 | 0.9652 | 0.9877 | +0.0225 [+0.0150, +0.0294] | 0.9665 |
| Eigen | 2,750 | 0.9618 | 0.9881 | +0.0262 [+0.0199, +0.0340] | 0.9671 |
| DEOGEN2 | 2,811 | 0.9588 | 0.9879 | +0.0291 [+0.0220, +0.0359] | 0.9664 |
| GPN-MSA | 2,919 | 0.9556 | 0.9882 | +0.0326 [+0.0252, +0.0407] | 0.9668 |
| CADD | 2,919 | 0.9538 | 0.9882 | +0.0344 [+0.0273, +0.0419] | 0.9668 |
| PolyPhen-2 HVAR | 2,766 | 0.9535 | 0.9876 | +0.0341 [+0.0252, +0.0408] | 0.9662 |
| MisFit | 2,844 | 0.9508 | 0.9878 | +0.0370 [+0.0286, +0.0474] | 0.9660 |
| FATHMM-XF | 2,746 | 0.9508 | 0.9880 | +0.0371 [+0.0299, +0.0452] | 0.9669 |
| MVP | 2,769 | 0.9473 | 0.9881 | +0.0409 [+0.0319, +0.0497] | 0.9679 |
| PROVEAN | 2,737 | 0.9437 | 0.9888 | +0.0451 [+0.0376, +0.0532] | 0.9672 |
| popEVE | 2,775 | 0.9431 | 0.9880 | +0.0449 [+0.0362, +0.0548] | 0.9671 |
| M-CAP | 2,721 | 0.9426 | 0.9871 | +0.0445 [+0.0352, +0.0564] | 0.9687 |
| ESM1b | 2,919 | 0.9406 | 0.9882 | +0.0475 [+0.0368, +0.0594] | 0.9668 |
| SIFT | 2,693 | 0.9394 | 0.9891 | +0.0498 [+0.0400, +0.0600] | 0.9672 |
| MetaSVM | 2,881 | 0.9394 | 0.9880 | +0.0486 [+0.0378, +0.0619] | 0.9665 |
| PolyPhen-2 HDIV | 2,766 | 0.9391 | 0.9876 | +0.0485 [+0.0369, +0.0578] | 0.9662 |
| MetaLR | 2,881 | 0.9375 | 0.9880 | +0.0505 [+0.0406, +0.0625] | 0.9665 |
| MutFormer | 2,872 | 0.9361 | 0.9881 | +0.0520 [+0.0423, +0.0626] | 0.9670 |
| SIFT4G | 2,813 | 0.9318 | 0.9886 | +0.0568 [+0.0470, +0.0673] | 0.9675 |
| MutationAssessor | 2,571 | 0.9301 | 0.9877 | +0.0577 [+0.0454, +0.0719] | 0.9654 |
| LIST-S2 | 2,859 | 0.9286 | 0.9879 | +0.0592 [+0.0480, +0.0721] | 0.9658 |
| DANN | 2,900 | 0.9104 | 0.9881 | +0.0777 [+0.0665, +0.0892] | 0.9668 |
| PrimateAI | 2,818 | 0.9037 | 0.9881 | +0.0845 [+0.0747, +0.0972] | 0.9665 |
| MPC | 2,511 | 0.8775 | 0.9891 | +0.1116 [+0.0881, +0.1391] | 0.9676 |
