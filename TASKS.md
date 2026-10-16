# Task List

The project is split into small end-to-end tasks. Each task is finished as **code + tests + summary output** and ends
with a concrete result the user can try on their own machine.

**Current task:** 8.2b ESM features for in-frame indels

Rules:
- Tasks are done in order. The next task does not start until the current one is finished.
- A task is not finished until its **Acceptance** items are met.
- When the user confirms the **Try it** step: tick the box, update "Current task", and add at most 3 lines to
  "Completion Notes" at the end.
- No features from later tasks are added while working on a task.
- If a task turns out larger than expected, split it into subtasks after asking the user.
- Long-running steps (downloads, ESM, training) are first validated on a small sample (`--chrom 22`, `--limit`),
  then run in full in the background; the assistant reads only the summary JSON.
- Tasks marked 🔍 invoke the **leakage-auditor**.

---

## Phase 0 — Foundation

- [x] **0.1 Project skeleton, environment and Graphify**
  - `git init`, `.gitignore` (list in CLAUDE.md Section 14)
  - `uv` project: Python 3.12, `pyproject.toml`, `src/healthcare/` package skeleton, ruff, pytest
  - PyTorch installed with CUDA
  - `uv run healthcare doctor`: Python/PyTorch/CUDA versions, GPU name and memory, free space on drive A:
  - `config/sources.toml` skeleton, `README.md` draft (English + Turkish sections, disclaimer included)
  - **Graphify installation (inside this task, never postponed):** `uv tool install graphifyy` →
    `graphify install --project` → `graphify claude install` → `graphify hook install` → `graphify update .`
  - Acceptance: `doctor` sees the GPU with CUDA; `pytest` and `ruff check` are clean; the Graphify graph is built and
    excludes `data/` and `artifacts/`
  - Try it: `uv run healthcare doctor` shows the RTX 4070 Ti and "CUDA: available"

- [x] **0.2 Shared infrastructure**
  - Download helper: resume, checksum verification, skip if present, `--force`
  - File logging (`logs/`), small summary JSON writer, path settings, `--chrom` / `--limit` sample mode
  - Acceptance: an interrupted download resumes; a summary JSON is produced; tests pass
  - Try it: a small source file (e.g. the MANE summary table) downloads, is skipped on the second run, and the summary
    JSON is read

## Phase 1 — Data and Our Own Variant Annotator

- [x] **1.1 Reference genome and MANE Select**
  - GRCh38 primary assembly FASTA (pyfaidx index), MANE Select GTF/GFF and protein FASTA (versions pinned)
  - Transcript table (Parquet): gene, strand, exons, CDS boundaries
  - Acceptance: protein translated from the CDS matches the MANE protein FASTA for ≥99% of transcripts; mismatches
    listed in the summary
  - Try it: the summary JSON shows transcript count and translation agreement; a sample gene's (e.g. BRCA1) protein
    sequence is verified

- [x] **1.2 ClinVar labels** 🔍
  - ClinVar GRCh38 VCF: October 2025 release for training + October 2026 release for the temporal test
    (variants first added in between); variants whose label changed between the two releases go to a separate
    reclassification table, never into the test set
  - Labeling rules from CLAUDE.md Section 7; REF check against the reference genome; deduplication
  - VUS records in a separate table; HGVS protein and molecular consequence fields kept in separate columns for
    validation only
  - Acceptance: unit tests for labeling rules; summary JSON with class × stars × type distribution
  - Try it: the summary shows how many pathogenic / benign / VUS records exist and their type distribution

- [x] **1.3 Our own variant annotator: SNVs**
  - Genomic SNV → MANE transcript → codon → amino acid change and protein position
  - SNV types: missense, synonymous, stop-gained, stop-lost, start-lost, canonical splice (±1/2), intronic / out of scope
  - Unit tests for both strands, codon boundaries and exon edges
  - Acceptance: missense agreement with the ClinVar HGVS protein field ≥99%; mismatching examples written to the
    summary and reviewed; `/code-review` on the task diff passes
  - Try it: the agreement report is read; a few known variants are run through the CLI and show the expected amino
    acid change

- [x] **1.4a gnomAD allele frequencies and gene constraint** (exomes; raw kept on A:)
  - gnomAD v4.x exome (and genome if needed) sites VCFs streamed chromosome by chromosome → Parquet
    (chrom, pos, ref, alt, AC, AN, AF, nhomalt, grpmax AF)
  - Gene constraint table (LOEUF, missense Z)
  - Raw files deleted only with the user's approval and `--delete-raw`
  - Acceptance: validated on `--chrom 22`, then a full background run; disk usage in the summary
  - Try it: a known variant's frequency comes back from a DuckDB query in milliseconds

- [x] **1.4b Exome coverage signal** (replaces targeted genomes, see O.5)
  - Exome allele number at every MANE CDS ±20 bp position (`coverage.parquet`), so "not observed" can be told apart
    from "not measured"
  - Acceptance: coverage built by `build-gnomad`; low-coverage share of labeled missense measured
  - Try it: the share of labeled missense variants at low-coverage positions is read from the summary

- [x] **1.5 Evolutionary conservation**
  - UCSC hg38 phyloP100way (and phastCons100way) wigFix files → Parquet keeping only MANE CDS ±20 bp positions
  - Acceptance: chromosome 22 validation; known-position test; disk usage in the summary
  - Try it: a position's phyloP value is queried

- [x] **1.6 UniProt protein region features**
  - Swiss-Prot curated features (domain, active site, binding site, transmembrane, signal peptide, disulfide, etc.)
  - Mapping to MANE proteins (direct when sequences are identical; otherwise excluded and counted)
  - Acceptance: mapped protein ratio in the summary; unit test
  - Try it: a sample protein's region map is shown as a table

## Phase 2 — Dataset and Split

- [x] **2.1 Missense feature table** 🔍
  - Labeled missense variants + frequency, conservation, gene constraint, UniProt regions, amino acid physicochemistry
    (Grantham, BLOSUM62, hydrophobicity, volume, charge change), relative protein position
  - No ESM yet
  - Missing gnomAD frequency handled explicitly (fill 0 for not observed and/or a flag); the no-allele-frequency
    ablation drops every frequency-derived column (`af`, `ac`, `an`, `nhomalt`, `af_grpmax`, `grpmax`, `filter_pass`,
    missingness flags)
  - Acceptance: per-feature missing-value rate in the summary; no column comes from a forbidden source or from ClinVar
    classification fields
  - Try it: the feature table schema and missing-value report are read

- [x] **2.2 Gene-disjoint split and temporal test** 🔍
  - Gene-disjoint train / validation / test with a fixed seed (about 70 / 15 / 15, class balance considered)
  - Temporal test: variants added to ClinVar after the training release and absent from training
  - Split on the MANE gene assigned by our annotator (never ClinVar's GENEINFO); genes sharing any coding position
    are merged into one split group (connected components); test that no `(chrom,pos)` appears in two sets
  - Watch gene memorization: `f_protein_length` + `f_gene_*` identify most genes; compare train/val gap and an
    ablation without gene-level features in 3.1
  - Temporal metrics also reported on the gene-disjoint subset (temporal variants whose genes are in validation/test),
    since most temporal genes also appear in training
  - Split file + hash; automated tests verifying gene and variant overlap
  - Acceptance: overlap tests pass; no 🔴 findings from the leakage-auditor
  - Try it: the split summary shows gene / variant / class counts per set and overlap = 0

## Phase 3 — Baseline Model

- [x] **3.1 LightGBM baseline**
  - Train on the training set, early stopping on validation; metrics **on the validation set only**
  - Small hyperparameter sweep selected on validation only
  - Feature importance with shap; ablation without allele frequency
  - Acceptance: validation metrics and importance plot under `reports/`; test set not opened
  - Try it: validation ROC-AUC / PR-AUC / MCC and top features are read

## Phase 4 — ESM-2 (Protein Tower)

- [x] **4.1 ESM-2 setup and single-protein trial**
  - `facebook/esm2_t33_650M_UR50D` on GPU; wild-type probabilities and per-residue representations for one protein
  - Windowing strategy for proteins longer than 1022 amino acids
  - Acceptance: runtime and GPU memory measured; long-protein test passes
  - Try it: one protein's output and compute time appear in the summary

- [x] **4.2 Precompute for all MANE proteins**
  - Resumable background batch job; storage format and size (fp16, dimensionality reduction if needed) decided at the
    start of the task (ask the user if a new library is needed)
  - Acceptance: processed protein count, error count and total disk usage in the summary
  - Try it: ESM output for any protein is read from disk instantly

- [x] **4.3 ESM features and contribution** 🔍
  - Log-likelihood ratio we compute (wild-type vs mutant amino acid), position entropy, representation around the position
  - Scalar ESM features added to LightGBM to measure their contribution on validation
  - Acceptance: with/without ESM validation comparison under `reports/`
  - Try it: ESM's metric contribution is shown in a table

## Phase 5 — Dual Tower Model

- [x] **5.1 Dual tower model**
  - PyTorch: protein tower (ESM representation) + measurement tower + fusion head; tower dropout; variant type input
  - Overnight sweep (12–16 h budget, measured first on 3–4 trial runs): random or Bayesian search (Optuna only if
    the user approves the library) over hidden size (e.g. d128–d512), depth, dropout, tower-dropout rate, ESM input
    form, learning rate, batch size; one seed per config, all runs appended to one results table
  - Confirmation round: top ~10 configs re-trained with 5 seeds; winner = best validation mean with low seed variance
    (guards against winner's curse); optional gene-grouped 5-fold CV inside training data
  - Early stopping on validation
  - Acceptance: validation metrics reported side by side with LightGBM; `/code-review` on the task diff passes
  - Try it: dual tower vs LightGBM validation table is read

- [x] **5.2 Calibration and the "Uncertain" band**
  - Model = average of LightGBM and dual tower probabilities (decided 2026-10-13, see CLAUDE.md Decisions)
  - Calibration on validation (temperature scaling or isotonic); Likely pathogenic / Uncertain / Likely benign thresholds
    at 90% precision on each side (ACMG-style); report the share of variants left Uncertain
  - Reliability plot, Brier, ECE
  - Acceptance: thresholds and calibration come from the validation set only
  - Try it: calibration plot and chosen thresholds are read

- [x] **5.3 ONNX export**
  - Dual tower + calibration → ONNX; inference with ONNX Runtime
  - Acceptance: PyTorch / ONNX parity test (within tolerance); single-variant CPU latency measured
  - Try it: sample prediction and latency via ONNX

## Phase 6 — Final Evaluation (Missense)

- [x] **6.1 Opening the test set** 🔍 (before the test set is opened)
  - Gene-disjoint test + temporal test metrics, gene-level bootstrap confidence intervals
  - Ablations: without allele frequency, without ESM, LightGBM vs dual tower
  - ROC, PR and calibration plots; `reports/missense_v1.md`
  - Acceptance: no model changes after testing (if any, stated explicitly in the report)
  - Try it: the final report is read

## Phase 7 — Prediction Service and Web Interface

- [x] **7.1 Prediction pipeline and CLI**
  - Input: genomic (`17:43063931 G>A`) or protein (`BRCA1 p.Arg1699Trp`, mapped to genomic coordinates by our annotator)
  - DuckDB feature lookup → ONNX prediction → explanation (features most responsible for the decision)
  - `uv run healthcare predict ...`; clear message for out-of-scope variants
  - Acceptance: end-to-end tests with known variants; single query in milliseconds
  - Try it: a few known pathogenic and benign variants are queried from the CLI

- [x] **7.2 FastAPI**
  - `/api/health`, `/api/predict`, `/api/predict/batch` (size- and row-limited VCF upload), `/api/model-card`
  - pydantic schemas; `127.0.0.1` only; CORS limited to the local Vite origin; uploads never written to disk
  - Acceptance: API tests (valid / invalid input, size limit); `/code-review` on the task diff passes
  - Try it: a single-variant query from the `/docs` page

- [x] **7.3 Web skeleton**
  - `web/`: React + TypeScript + Vite; API types via `openapi-typescript`; API client
  - Single-variant query page (no design yet)
  - Acceptance: `tsc --noEmit` and `vite build` clean; a backend schema change produces a type error
  - Try it: a variant is entered in the browser and the result appears

- [x] **7.4 Design pass**
  - Design direction and token system with UI UX Pro Max (color, typography, spacing; dark mode)
  - Result card: probability, three-way decision (including "Uncertain"), chart of the features driving the decision,
    disclaimer
  - English UI by default with an optional Turkish language toggle (ask before adding an i18n library)
  - Ask the user before adding a charting library
  - Finish with a web-design-guidelines audit; follow vercel-react-best-practices
  - Acceptance: no critical audit findings remain
  - Try it: the interface is reviewed in light and dark mode, in English and Turkish

- [x] **7.5 Batch query and model card pages**
  - VCF upload → results table → CSV download
  - Model card: data sources and versions, per-type metrics, baseline comparison, limitations
  - Acceptance: clean error messages for large / invalid files
  - Try it: a sample VCF is uploaded and the results are downloaded as CSV

## Phase 8 — Expansion (Other Coding Variant Types)

- [x] **8.1 Extend the annotator to indels**
  - Left-normalization, frameshift, in-frame insertion / deletion, canonical splice effects of indels
  - Acceptance: per-type agreement report against ClinVar HGVS; unit tests; `/code-review` on the task diff passes
  - Try it: known frameshift and in-frame variants are run through the CLI

- [x] **8.2a Type-specific features** 🔍
  - mRNA with UTRs from the reference genome + MANE exons (no new source)
  - Truncation position, fraction of protein retained, novel residues, NMD rule (last exon / 55 nt before the last
    exon-exon junction), stop-loss extension, next in-frame Met for start loss, splice site side/position and the
    frame of the skipped exon, distance to the nearest exon boundary, indel length
  - One table for all Phase 2 types (synonymous included; sampling is decided in 8.3); the Phase 1 missense table is
    not touched; a variant scored in more than one gene is dropped
  - Acceptance: per-type missing-value report; leakage-auditor has no 🔴 findings
  - Try it: the per-type feature report is read

- [ ] **8.2b ESM features for in-frame indels**
  - ESM-based features for in-frame indels; protein tower masked for types without them
  - Acceptance: per-type missing-value report including ESM columns
  - Try it: the feature report is read

- [ ] **8.3 Joint model training** 🔍
  - Dual tower retrained or fine-tuned on all types (including missense); LightGBM baseline updated
  - The Phase 1 gene split is reused unchanged
  - Acceptance: per-type validation metrics; missense performance does not drop meaningfully versus Phase 1
  - Try it: the per-type validation table is read

- [ ] **8.4 Final evaluation v2 and integration** 🔍 (before the test set is opened)
  - Per-type test + temporal test metrics, ablations, `reports/all_types_v2.md`
  - New ONNX; prediction pipeline and interface support the new types and indel input
  - Acceptance: rules from 6.1; interface and API tests pass
  - Try it: missense, nonsense, frameshift and synonymous examples are queried in the interface

## Phase 9 — Delivery

- [ ] **9.1 Final README and model card**
  - Single README with English and Turkish sections: results, architecture diagram, screenshots, data sources and
    licenses, reproduction steps, limitations, disclaimer
  - Try it: the README steps can be followed end to end on a clean clone

- [ ] **9.2 Demo pack and release**
  - Full data tables are too large, so a small demo pack for a selected gene subset (ONNX model + lookup tables for
    those genes) that fits GitHub Releases size limits
  - Going public is the user's decision
  - Try it: the interface runs from the demo pack without downloading the full data

## Optional (only if the user asks)

- [ ] **O.1 VEP cross-check in Docker:** compare our annotator's output against VEP's consequence field on a sample.
  VEP score plugins are not used.
- [ ] **O.2 Masked-marginal ESM score:** compute the log-likelihood ratio by masking each variant position and compare
  with the wild-type method.
- [ ] **O.3 ESM-2 3B comparison:** measure the larger model's contribution and cost.
- [ ] **O.4 Our own small masked protein model:** compare ESM with a small protein language model trained from scratch.
- [ ] **O.5 Targeted gnomAD genome frequencies:** fetch genome records for MANE CDS ±20 bp with the tabix reader
  (`data/tabix.py`) and in-window bisection over BGZF blocks (16 kb leaf bins are ~3 MB each, so whole-bin fetching
  is not viable). Only worth it if results are poor at low-coverage positions (<1% of labeled missense).

---

## Completion Notes

(At most 3 lines per finished task: date, task, key result / decision.)
- 2026-10-11 · 0.1 · Python 3.12 (uv), PyTorch 2.14.1+cu130 on RTX 4070 Ti; `healthcare doctor` CLI; Graphify installed (AST-only `graphify update .`, hooks active).
- 2026-10-11 · 0.2 · Stdlib downloader (resume, sha256, skip, `--force`), `healthcare download`, summaries in `logs/summaries/`; first source MANE v1.5 summary pinned. `--chrom`/`--limit` deferred to 1.1.
- 2026-10-11 · 1.1 · 19,229 MANE v1.5 transcripts on 1–22/X/Y, CDS translation match 99.995% (non-ATG starts and selenocysteine accepted, flagged in `translation_note`); 60 genes annotated only on fix/alt patches are excluded.
- 2026-10-11 · 1.2 · ClinVar 2025-10-06: 1.52M labels (282k P / 1.24M B), 2.07M VUS; temporal (new in 2026-10-04): 186,768 (44k P); 3,392 VUS resolved, 28 flips. Leakage audit PASS; metadata columns prefixed `clinvar_`/`validation_`.
- 2026-10-11 · 1.2 · Fixed `.gitignore`/`.graphifyignore` (`data/` was hiding `src/healthcare/data/`; `mane.py` was missing from the 1.1 commit).
- 2026-10-11 · 1.3 · SNV annotator: missense HGVS agreement 175,203/175,203 (100%), consequence agreement 99.999%; missense train set 183,375 (53,371 P / 130,004 B, 14,934 genes). `/code-review`: 8 findings, all fixed.
- 2026-10-11 · 1.4a · gnomAD v4.1 exomes: 183.7M variants → 2.2 GB Parquet (24 chr, md5-verified), DuckDB lookup ~40 ms; 17,481 MANE constraint rows (no `mis_pphen`). Missense: 91% of benign vs 33% of pathogenic seen in gnomAD. Leakage audit PASS.
- 2026-10-12 · 1.4b · Targeted genomes dropped: tabix leaf bins (16 kb ≈ 3 MB) make it ~full download; exome covers 99.96% of 41.1M CDS±20 positions, only 0.7% B / 0.3% P labeled missense at low coverage. `coverage.parquet` added instead; genomes moved to O.5.
- 2026-10-12 · 1.5 · UCSC phyloP100way + phastCons100way for 41.1M CDS±20 positions (99.99% filled) → 203 MB; lookup ~8 ms. Labeled missense median phyloP: benign 0.90 vs pathogenic 7.64.
- 2026-10-12 · 1.6 · UniProt 2026_03: 17,384/19,229 MANE proteins matched by exact sequence (90.4%), 461k curated features (VARIANT/MUTAGEN and purely automatic evidence excluded). Key-region share: benign 17.5% vs pathogenic 36.1%.
- 2026-10-12 · 2.1 · Missense features: train 183,343 (53,354 P, 14,934 genes), temporal 18,205, VUS 1.88M; 40 `f_` columns. Multi-gene variants dropped (32 train). Leakage audit PASS; absence from gnomAD is a strong signal (63% P vs 8% B).
- 2026-10-12 · 2.2 · Split pinned (sha256): train 128,340 / val 27,502 / test 27,501 variants, 29.1% P each; all 19,229 MANE genes assigned, CDS-overlapping genes grouped (32 merged groups), 0 overlaps. Test + temporal only in `splits/locked/`, readable via `load_split(allow_locked=True)`. Temporal test_gene subset: 2,919. Leakage audit PASS (5 findings fixed).
- 2026-10-12 · 3.1 · LightGBM baseline (30-trial sweep, early stopping on 10% held-out train genes): validation ROC-AUC 0.985 [0.982–0.987], PR-AUC 0.963, MCC 0.873. No-AF 0.944, no-gene-level 0.983 (no gene memorization). Test not opened.
- 2026-10-12 · 4.1 · ESM-2 650M (pinned revision, safetensors, fp16) on RTX 4070 Ti: BRCA1 1,863 aa in 0.12 s, TTN 35,991 aa (70 windows) in 2.6 s, peak 1.35 GiB. Wild-type-marginal LLR: R1699W −13.9; open-split BRCA1 median LLR benign −3.8 vs pathogenic −11.2.
- 2026-10-12 · 4.2 · ESM-2 store: 19,228 proteins (11.1M residues), 0 errors, 11 min; full fp16 embeddings 28.5 GB + log-probs 0.45 GB per protein `.npy`, memory-mapped read ~9 ms.
- 2026-10-12 · 4.3 · 5 `f_esm_*` features (45 total). LightGBM validation ROC-AUC 0.985→0.990, PR-AUC 0.963→0.976; no-AF 0.944→0.962 (PR +0.045); ESM-only 0.904. Split pin now hashes (variant, gene, label) identity; assignments unchanged. Leakage audit PASS.
- 2026-10-13 · 5.1 · Optuna sweep (40% random + 60% TPE, no pruning): 4,645 trials in 7.7 h (stopped by a CUDA driver error at 04:09; now handled). Winner #2636 (d512, window3 ESM): validation ROC-AUC 0.9890±0.0003 vs LightGBM 0.9899; no-AF 0.9637 vs 0.9617; ECE 0.016 vs 0.024. Sweep best 0.9903 shrank to 0.9890 over 5 seeds (winner's curse). `/code-review`: 7/8 fixed.
- 2026-10-13 · 5.2 · Platt on ensemble logit (out-of-fold ECE 0.0149→0.0042, Brier 0.0299→0.0288). Decision on calibrated probability ≥0.90 / ≤0.10 (set-level 90% precision made the bands overlap). Out-of-fold validation: 23.2% LP (98.1% P), 65.1% LB (0.9% P), 11.7% Uncertain.
- 2026-10-13 · 5.3 · Single `missense_v1.onnx` (14 MB): LightGBM + dual tower + Platt, preprocessing in-graph. Validation parity: dual tower 7e-7, final 7.5e-3 (float32 trees vs f_delta_* noise splits), 0 decision changes, ECE 0.0037. CPU single variant median 0.48 ms.
- 2026-10-14 · 6.1 · Test opened once (leakage audit PASS; ONNX sha256 pinned in report). Gene-disjoint test ROC-AUC 0.9902 [0.9884–0.9921], PR-AUC 0.9765, ECE 0.0044, 12.2% Uncertain; temporal 0.9887 (gene-disjoint subset 0.9882) but ECE 0.033 (overconfident in 0.3–0.8). No-AF 0.9687, no-ESM 0.9851. Model unchanged.
- 2026-10-16 · 7.1 · `healthcare predict` (genomic or `GENE p.Xxx123Yyy`) reuses the training feature functions; `join_sources` now prefilters by chrom/pos range (one variant 2.5 s → ~120 ms). `serve-check`: 0 feature mismatches and 0 probability diff on 100 validation variants, median 155 ms with explanation. `serve/` imports no torch.
- 2026-10-16 · 7.2 · `healthcare serve`: FastAPI on 127.0.0.1 (health, predict, plain-text VCF batch ≤2 MB / 500 variants in memory, model card), CORS for the Vite origin only; 30 API/serve tests (httpx dev-only). Serving now keeps nested-gene variants like training (416 validation variants were wrongly rejected); 500-variant batch 14.6 s. `/code-review`: 3 fixed, 2 deferred.
- 2026-10-16 · 7.3 · `web/`: React 19 + Vite 8 + TypeScript 5.9 (pinned: openapi-typescript needs TS 5); types from `web/openapi.json` (`healthcare openapi` → `npm run gen:api`), a pytest drift check fails until regenerated; renaming a backend field breaks `tsc`. Single-variant page verified in the browser against the live API.
- 2026-10-16 · 7.4 · Red underground-lab HUD theme modelled on the RE remakes and Requiem cutscenes (own SVG: layered cables, a glossy black symbiote-like goo mass, red/purple/green specimen tanks with ringed metal caps and under-lights holding a pulsing veined tissue mass on feed tubes, a wall-anchored goo web and a glowing serum sample with virions, a CRT monitor that reacts to the decision; no game IP), light/dark, EN/TR via a typed dictionary (no i18n or chart library): decision badge with icon, banded probability gauge, diverging log-odds chart with a plain-language summary, measurements. System fonts only (offline). 12 Vitest tests; web-design-guidelines audit: no critical findings, 4 minor fixed.
- 2026-10-16 · 7.5 · Hash-routed pages (no router library): Batch (VCF) uploads the file as-is to `/api/predict/batch` after client checks (empty, >2 MB, compressed), shows a decision table, line errors and a CSV download (formula-safe cells); API 413/415/422 map to localized messages. Model card page lists thresholds, test/temporal metrics with bootstrap CIs, baseline/ablations, collapsed data versions and limitations. `web/public/sample.vcf` covers all three decisions, out of scope and a line error; 18 Vitest tests.
- 2026-10-16 · 8.1 · Indel annotator (`annotate/indel.py`): left-normalization, frameshift (HGVS short `fs`), in-frame del/dup/ins/delins with the protein 3' rule, stop gained/lost/retained, start lost, canonical splice; MNVs stay out of scope. ClinVar train: 192,880 indels, all already left-normalized; ClinVar-equivalent protein change (edits applied to the MANE protein) frameshift 99.97%, in-frame del 99.87%, ins 99.29%, delins 100% (exact strings 58–100%: ClinVar mixes ins/dup and repeat notation); indel consequence agreement 94.4%. `/code-review`: 5/6 fixed.
- 2026-10-16 · 8.2a · `build-phase2-features`: 9 Phase 2 types, 51 `f_` columns (18 type-specific: truncation, NMD 55 nt rule via mRNA built from GRCh38 + MANE exons, stop-loss extension, next Met, splice side/skipped-exon frame, boundary distance, indel length; span phyloP). Train 845,024 rows (73% synonymous), temporal 112,503, VUS 71,778; 303 multi-gene dropped; type-specific missing <1% except partial stop-loss deletions. Benign frameshifts escape NMD 48% vs 15%. Missense table unchanged (serve-check 0 diffs). Leakage audit PASS: 181 missense variants also scored as Phase 2 in another gene, to be settled in 8.3; type alone nearly separates labels, so report per type.
