# Task List

The project is split into small end-to-end tasks. Each task is finished as **code + tests + summary output** and ends
with a concrete result the user can try on their own machine.

**Current task:** 1.3 Our own variant annotator: SNVs

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

- [ ] **1.3 Our own variant annotator: SNVs**
  - Genomic SNV → MANE transcript → codon → amino acid change and protein position
  - SNV types: missense, synonymous, stop-gained, stop-lost, start-lost, canonical splice (±1/2), intronic / out of scope
  - Unit tests for both strands, codon boundaries and exon edges
  - Acceptance: missense agreement with the ClinVar HGVS protein field ≥99%; mismatching examples written to the
    summary and reviewed; `/code-review` on the task diff passes
  - Try it: the agreement report is read; a few known variants are run through the CLI and show the expected amino
    acid change

- [ ] **1.4 gnomAD allele frequencies and gene constraint**
  - gnomAD v4.x exome (and genome if needed) sites VCFs streamed chromosome by chromosome → Parquet
    (chrom, pos, ref, alt, AC, AN, AF, nhomalt, grpmax AF)
  - Gene constraint table (LOEUF, missense Z)
  - Raw files deleted only with the user's approval and `--delete-raw`
  - Acceptance: validated on `--chrom 22`, then a full background run; disk usage in the summary
  - Try it: a known variant's frequency comes back from a DuckDB query in milliseconds

- [ ] **1.5 Evolutionary conservation**
  - UCSC hg38 phyloP100way (and phastCons100way) wigFix files → Parquet keeping only MANE CDS ±20 bp positions
  - Acceptance: chromosome 22 validation; known-position test; disk usage in the summary
  - Try it: a position's phyloP value is queried

- [ ] **1.6 UniProt protein region features**
  - Swiss-Prot curated features (domain, active site, binding site, transmembrane, signal peptide, disulfide, etc.)
  - Mapping to MANE proteins (direct when sequences are identical; otherwise excluded and counted)
  - Acceptance: mapped protein ratio in the summary; unit test
  - Try it: a sample protein's region map is shown as a table

## Phase 2 — Dataset and Split

- [ ] **2.1 Missense feature table** 🔍
  - Labeled missense variants + frequency, conservation, gene constraint, UniProt regions, amino acid physicochemistry
    (Grantham, BLOSUM62, hydrophobicity, volume, charge change), relative protein position
  - No ESM yet
  - Acceptance: per-feature missing-value rate in the summary; no column comes from a forbidden source or from ClinVar
    classification fields
  - Try it: the feature table schema and missing-value report are read

- [ ] **2.2 Gene-disjoint split and temporal test** 🔍
  - Gene-disjoint train / validation / test with a fixed seed (about 70 / 15 / 15, class balance considered)
  - Temporal test: variants added to ClinVar after the training release and absent from training
  - Split on the MANE gene assigned by our annotator (never ClinVar's GENEINFO); variants mapping to several genes
    are handled explicitly
  - Temporal metrics also reported on the gene-disjoint subset (temporal variants whose genes are in validation/test),
    since most temporal genes also appear in training
  - Split file + hash; automated tests verifying gene and variant overlap
  - Acceptance: overlap tests pass; no 🔴 findings from the leakage-auditor
  - Try it: the split summary shows gene / variant / class counts per set and overlap = 0

## Phase 3 — Baseline Model

- [ ] **3.1 LightGBM baseline**
  - Train on the training set, early stopping on validation; metrics **on the validation set only**
  - Small hyperparameter sweep selected on validation only
  - Feature importance with shap; ablation without allele frequency
  - Acceptance: validation metrics and importance plot under `reports/`; test set not opened
  - Try it: validation ROC-AUC / PR-AUC / MCC and top features are read

## Phase 4 — ESM-2 (Protein Tower)

- [ ] **4.1 ESM-2 setup and single-protein trial**
  - `facebook/esm2_t33_650M_UR50D` on GPU; wild-type probabilities and per-residue representations for one protein
  - Windowing strategy for proteins longer than 1022 amino acids
  - Acceptance: runtime and GPU memory measured; long-protein test passes
  - Try it: one protein's output and compute time appear in the summary

- [ ] **4.2 Precompute for all MANE proteins**
  - Resumable background batch job; storage format and size (fp16, dimensionality reduction if needed) decided at the
    start of the task (ask the user if a new library is needed)
  - Acceptance: processed protein count, error count and total disk usage in the summary
  - Try it: ESM output for any protein is read from disk instantly

- [ ] **4.3 ESM features and contribution** 🔍
  - Log-likelihood ratio we compute (wild-type vs mutant amino acid), position entropy, representation around the position
  - Scalar ESM features added to LightGBM to measure their contribution on validation
  - Acceptance: with/without ESM validation comparison under `reports/`
  - Try it: ESM's metric contribution is shown in a table

## Phase 5 — Dual Tower Model

- [ ] **5.1 Dual tower model**
  - PyTorch: protein tower (ESM representation) + measurement tower + fusion head; tower dropout; variant type input
  - Hyperparameter sweep (hidden size e.g. d128/d256/d384, depth, dropout, tower-dropout rate) × several seeds,
    all logged to one results table; the winner is chosen on validation only and reported as mean ± std over seeds
  - Early stopping on validation
  - Acceptance: validation metrics reported side by side with LightGBM; `/code-review` on the task diff passes
  - Try it: dual tower vs LightGBM validation table is read

- [ ] **5.2 Calibration and the "Uncertain" band**
  - Calibration on validation (temperature scaling or isotonic); Likely pathogenic / Uncertain / Likely benign thresholds
  - Reliability plot, Brier, ECE
  - Acceptance: thresholds and calibration come from the validation set only
  - Try it: calibration plot and chosen thresholds are read

- [ ] **5.3 ONNX export**
  - Dual tower + calibration → ONNX; inference with ONNX Runtime
  - Acceptance: PyTorch / ONNX parity test (within tolerance); single-variant CPU latency measured
  - Try it: sample prediction and latency via ONNX

## Phase 6 — Final Evaluation (Missense)

- [ ] **6.1 Opening the test set** 🔍 (before the test set is opened)
  - Gene-disjoint test + temporal test metrics, gene-level bootstrap confidence intervals
  - Ablations: without allele frequency, without ESM, LightGBM vs dual tower
  - ROC, PR and calibration plots; `reports/missense_v1.md`
  - Acceptance: no model changes after testing (if any, stated explicitly in the report)
  - Try it: the final report is read

## Phase 7 — Prediction Service and Web Interface

- [ ] **7.1 Prediction pipeline and CLI**
  - Input: genomic (`17:43045712 G>A`) or protein (`BRCA1 p.Arg1699Trp`, mapped to genomic coordinates by our annotator)
  - DuckDB feature lookup → ONNX prediction → explanation (features most responsible for the decision)
  - `uv run healthcare predict ...`; clear message for out-of-scope variants
  - Acceptance: end-to-end tests with known variants; single query in milliseconds
  - Try it: a few known pathogenic and benign variants are queried from the CLI

- [ ] **7.2 FastAPI**
  - `/api/health`, `/api/predict`, `/api/predict/batch` (size- and row-limited VCF upload), `/api/model-card`
  - pydantic schemas; `127.0.0.1` only; CORS limited to the local Vite origin; uploads never written to disk
  - Acceptance: API tests (valid / invalid input, size limit); `/code-review` on the task diff passes
  - Try it: a single-variant query from the `/docs` page

- [ ] **7.3 Web skeleton**
  - `web/`: React + TypeScript + Vite; API types via `openapi-typescript`; API client
  - Single-variant query page (no design yet)
  - Acceptance: `tsc --noEmit` and `vite build` clean; a backend schema change produces a type error
  - Try it: a variant is entered in the browser and the result appears

- [ ] **7.4 Design pass**
  - Design direction and token system with UI UX Pro Max (color, typography, spacing; dark mode)
  - Result card: probability, three-way decision (including "Uncertain"), chart of the features driving the decision,
    disclaimer
  - English UI by default with an optional Turkish language toggle (ask before adding an i18n library)
  - Ask the user before adding a charting library
  - Finish with a web-design-guidelines audit; follow vercel-react-best-practices
  - Acceptance: no critical audit findings remain
  - Try it: the interface is reviewed in light and dark mode, in English and Turkish

- [ ] **7.5 Batch query and model card pages**
  - VCF upload → results table → CSV download
  - Model card: data sources and versions, per-type metrics, baseline comparison, limitations
  - Acceptance: clean error messages for large / invalid files
  - Try it: a sample VCF is uploaded and the results are downloaded as CSV

## Phase 8 — Expansion (Other Coding Variant Types)

- [ ] **8.1 Extend the annotator to indels**
  - Left-normalization, frameshift, in-frame insertion / deletion, canonical splice effects of indels
  - Acceptance: per-type agreement report against ClinVar HGVS; unit tests; `/code-review` on the task diff passes
  - Try it: known frameshift and in-frame variants are run through the CLI

- [ ] **8.2 Type-specific features** 🔍
  - Truncation position, fraction of protein lost, NMD rule (last exon / 50–55 nucleotide rule), splice distance
  - ESM-based features for in-frame indels; protein tower masked for types without them
  - Acceptance: per-type missing-value report
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

---

## Completion Notes

(At most 3 lines per finished task: date, task, key result / decision.)
- 2026-10-11 · 0.1 · Python 3.12 (uv), PyTorch 2.14.1+cu130 on RTX 4070 Ti; `healthcare doctor` CLI; Graphify installed (AST-only `graphify update .`, hooks active).
- 2026-10-11 · 0.2 · Stdlib downloader (resume, sha256, skip, `--force`), `healthcare download`, summaries in `logs/summaries/`; first source MANE v1.5 summary pinned. `--chrom`/`--limit` deferred to 1.1.
- 2026-10-11 · 1.1 · 19,229 MANE v1.5 transcripts on 1–22/X/Y, CDS translation match 99.995% (non-ATG starts and selenocysteine accepted, flagged in `translation_note`); 60 genes annotated only on fix/alt patches are excluded.
- 2026-10-11 · 1.2 · ClinVar 2025-10-06: 1.52M labels (282k P / 1.24M B), 2.07M VUS; temporal (new in 2026-10-04): 186,768 (44k P); 3,392 VUS resolved, 28 flips. Leakage audit PASS; metadata columns prefixed `clinvar_`/`validation_`.
- 2026-10-11 · 1.2 · Fixed `.gitignore`/`.graphifyignore` (`data/` was hiding `src/healthcare/data/`; `mane.py` was missing from the 1.1 commit).
