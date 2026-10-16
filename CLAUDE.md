# Healthcare — Genetic Variant Pathogenicity Prediction — Project Instructions

This file defines how the AI assistant must work in this repository.
Read this file and the task list at the start of every session. Treat them as the authoritative project context;
do not rediscover the project or re-litigate what is documented here.

@TASKS.md

Short notes on completed tasks live in the "Completion Notes" section at the end of `TASKS.md`.
There is no separate code map file; the **Graphify** graph fills that role (see Section 6).

---

# 1. Project Definition

- **What:** A from-scratch AI model that predicts whether germline genetic variants are **Pathogenic** or
  **Benign**, plus a local web interface that uses it.
- **Why:** Personal CV / portfolio project. The problem is the university-level problem of the TEKNOFEST 2026
  Artificial Intelligence in Healthcare Competition, but the team is **not competing**. The competition
  specification, report templates and format are out of scope and must not be used as references.
- **Data:** Public datasets only. The model and every feature are produced by us.
- **Product:** A model that runs on the user's machine (Windows 11) plus a web interface served on `localhost`.
  No cloud, no deployment, no mobile app, no user accounts.
- **Repository:** GitHub `Healthcare`. Created as private; the decision to make it public belongs to the user.
- **Disclaimer:** For research and educational use only; not for clinical diagnosis or treatment decisions.
  This statement must be visible in the README and in the interface.

---

# 2. MOST CRITICAL RULE: No Precomputed Pathogenicity Scores

The model learns its decision **only from raw measurements and features we compute ourselves**.
Using another model's "this variant is harmful" verdict is copying someone else's answer sheet and invalidates
the project's from-scratch claim.

## Forbidden (as input, feature, label, pre-filter, ensemble member or benchmark — never downloaded, never used)

CADD, REVEL, AlphaMissense, SpliceAI, PrimateAI / PrimateAI-3D, EVE, precomputed ESM1b/ESM variant scores,
SIFT, PolyPhen-2, MutationTaster, MutationAssessor, PROVEAN, FATHMM, MetaSVM, MetaLR, MetaRNN, BayesDel, ClinPred,
VEST, DANN, Eigen, M-CAP, MVP, gMVP, VARITY, and every other variant-effect predictor.
**dbNSFP** (the database aggregating these scores) is not used at all. Score-adding Ensembl VEP plugins are not used.

## Allowed sources

| Source | Purpose | Note |
|---|---|---|
| **ClinVar** (GRCh38 VCF, pinned release) | Labels + variant identity | Classification, review status, submitter count and similar fields **are never features**. HGVS / molecular consequence fields are used only to validate our own annotator. |
| **GRCh38 reference genome** (primary assembly FASTA) | Reference base, codons | |
| **MANE Select** (GTF/GFF + protein FASTA) | One canonical transcript per gene | |
| **gnomAD v4.x** allele frequencies | Population frequency | Count statistic, not a verdict. See the circularity warning below. |
| **gnomAD gene constraint** (LOEUF, missense Z) | Gene tolerance to damage | Count-based statistic; agreed with the user. |
| **UCSC phyloP / phastCons** (hg38, wigFix) | Evolutionary conservation | Statistic derived from cross-species alignments. |
| **UniProt Swiss-Prot** curated features | Domains, active sites, binding sites, transmembrane regions, etc. | Curated (manual) annotations only. |
| **ESM-2** (`facebook/esm2_t33_650M_UR50D`) | Protein tower | **We run it locally ourselves.** Pretrained weights only; a masked language model never trained on variant labels. |
| Amino acid physicochemical tables (Grantham, BLOSUM62, hydrophobicity, volume, charge) | Measurement tower | |

- Ask the user **before** adding any data source not listed here. Check whether it is a predictor trained on
  labeled variant data or derived from ClinVar.
- **gnomAD circularity warning:** some ClinVar Benign labels were assigned precisely because a variant is common
  in the population, so allele frequency can make the model look stronger than it is. Every final report also
  includes an ablation trained **without allele frequency**.

---

# 3. Scope

- Genome build: **GRCh38 only**. Chromosome names are normalized to `1..22, X, Y` (no `chr` prefix).
- Transcripts: **MANE Select only** (one transcript per gene).
- **Phase 1:** missense SNVs (the project's showcase).
- **Phase 2:** synonymous, nonsense (stop-gained), stop-lost, start-lost, canonical splice (intron ±1/2),
  frameshift, in-frame insertions/deletions.
- **Permanently out of scope:** deep intronic, non-canonical splice region, UTR, regulatory, intergenic,
  mitochondrial, structural variants / CNVs, somatic variants. The README states why.

---

# 4. Architecture

## Pipeline

```
download (data/raw) → process (data/interim) → our own annotator (annotate) → features
     → ESM precompute (esm) → gene-disjoint split (split) → training (LightGBM baseline + dual tower)
     → calibration → evaluation (reports/) → ONNX export → FastAPI → React interface
```

## Model

- **Single model, dual tower:**
  - **Protein tower:** ESM-2 outputs (representation around the variant position + scalars we compute, such as
    the log-likelihood ratio).
  - **Measurement tower:** conservation, allele frequency, gene constraint, amino acid physicochemistry,
    position within the protein, UniProt features, variant type, truncation position, etc.
  - The towers are fused into a single pathogenicity probability.
- **Tower dropout:** the protein tower is randomly masked during training so the model learns to decide without
  protein information. This lets Phase 2 types (synonymous, nonsense, etc.) run through the same model.
- **Variant type** is a model input; type-specific features go into type-specific input blocks.
- **Phase 2 retraining:** when new types are added the model is trained or fine-tuned **together with the missense
  data** (training only on new types causes forgetting). The result is still a single ONNX file.
- **Baseline:** LightGBM on the same features, reported side by side with the dual tower in every report.
- **Output:** calibrated pathogenicity probability + a three-way decision: *Likely pathogenic / Uncertain /
  Likely benign*. Thresholds are chosen on the validation set only.

## Inference

- ESM outputs are precomputed **once** per MANE protein and stored on disk.
- Frequency and conservation data live in **Parquet** files partitioned by chromosome and sorted by position,
  queried with **DuckDB**. Target: a single-variant query in milliseconds.
- The model runs on **ONNX Runtime**; the serving side (`serve/`) does not depend on PyTorch.
- FastAPI binds to `127.0.0.1` only.

## Module boundaries (`src/healthcare/`)

| Module | Responsibility | Rule |
|---|---|---|
| `data/` | Downloads, raw files → Parquet | Raw files are never modified. |
| `annotate/` | Genomic coordinate → MANE transcript → codon → amino acid change, variant type | Pure functions; no I/O; unit tested. Our own annotator; no VEP dependency. |
| `features/` | Feature table generation | Any label-derived statistic is computed from the training set only. |
| `esm/` | ESM-2 precompute and readers | GPU; resumable. |
| `split/` | Gene-disjoint and temporal splits | The split file is pinned and never silently regenerated. |
| `models/` | LightGBM baseline, dual tower, calibration, ONNX export | Never reads raw files, only `data/processed`. |
| `evaluate/` | Metrics, ablations, plots, reports | The test set is touched only in final evaluation tasks. |
| `serve/` | Prediction pipeline + FastAPI | Never imports PyTorch. |
| `cli.py` | `uv run healthcare <command>` entry point | Stays thin; logic lives in modules. |

Modules do not touch each other's internal file formats directly; they communicate through Parquet schemas and
function interfaces.

## Folder layout

```
HealthPproject/                 (GitHub repository name: Healthcare)
├── CLAUDE.md, TASKS.md
├── README.md                   → single README, English section + Turkish section
├── pyproject.toml, uv.lock
├── config/sources.toml         → data source URLs, versions, checksums (committed)
├── src/healthcare/             → Python package (modules above)
├── tests/                      → pytest; only small sample data under tests/fixtures/
├── web/                        → React + TypeScript + Vite interface
├── reports/                    → metric summaries, plots, model card (committed, small files)
├── data/                       → raw/, interim/, processed/ (never committed)
├── artifacts/                  → model weights, ONNX, ESM outputs (never committed)
└── logs/                       → long-running job logs (never committed)
```

The project lives on drive A: (Samsung 990 Pro, ~1 TB free); `data/` and `artifacts/` are inside the project.

---

# 5. Technology

| Layer | Technology |
|---|---|
| Language / environment | Python 3.12, **uv** |
| Deep learning | PyTorch (CUDA) — RTX 4070 Ti 12 GB, i5-13600KF |
| Protein language model | ESM-2 650M via Hugging Face `transformers` |
| Data processing | polars, pyarrow, DuckDB, pyfaidx |
| Classical ML | LightGBM, scikit-learn, shap (baseline explanations) |
| Model export / inference | onnx, onnxruntime, onnxmltools (LightGBM → ONNX, export only) |
| Backend | FastAPI, uvicorn, pydantic |
| Frontend | React + TypeScript + Vite; types generated from the FastAPI schema with `openapi-typescript` |
| Testing / quality | pytest, ruff; frontend `tsc --noEmit`, `vite build`, Vitest |

- Ask the user before adding any library not listed here, and explain why.
- **Windows compatibility:** avoid libraries without Windows wheels (pysam, pyBigWig, cyvcf2, etc.).
  Parse VCF / wigFix with polars or plain Python. Docker is only for optional cross-checks.
- Do not use a technology merely because it is listed; add it only when the current task requires it.
- Do not suggest paid services or accounts.

---

# 6. Plugins and Skills

Use **only** the following. Do not use similarly named forks or versions.
Where they conflict, **the rules in this file take precedence.**

| Plugin | Source | Purpose | Rule |
|---|---|---|---|
| **Context7** | https://context7.com (project `.mcp.json`) | Current docs for PyTorch, transformers, polars, DuckDB, LightGBM, ONNX Runtime, FastAPI, React, Vite | Check current usage with Context7 before using a library API. |
| **ponytail** | https://github.com/DietrichGebert/ponytail | Avoid unnecessary code and libraries | Use `full` mode, never `ultra`. Never cut validation, error handling or tests. |
| **Graphify** | https://github.com/Graphify-Labs/graphify (PyPI: `graphifyy`) | Code knowledge graph; query the graph instead of reading files one by one | Installed and used per the **mandatory rule** below. |
| **UI UX Pro Max** | https://github.com/nextlevelbuilder/ui-ux-pro-max-skill | Interface color, typography, layout, UX decisions | Design decisions only. Ask before adding any extra UI library. |
| **vercel-react-best-practices** | https://github.com/vercel-labs/agent-skills | React performance and correctness rules | Apply rules that fit React + Vite; skip Next.js-specific rules (RSC, App Router, server actions, Next cache). |
| **web-design-guidelines** | https://github.com/vercel-labs/agent-skills | Accessibility / UX audit | Use as an **audit checklist** at the end of the design pass. |
| **vercel-composition-patterns** | https://github.com/vercel-labs/agent-skills | React component composition patterns | Apply to React + Vite components from 7.3 on; skip Next.js-specific guidance. |
| **building-components** | https://github.com/vercel/components.build | Accessible, composable, themeable components | Component code only; color, typography and layout decisions stay with UI UX Pro Max. |
| **before-and-after** | https://github.com/vercel-labs/before-and-after | Presents before/after screenshots of a page change | Used when the user restyles the site (7.4+). Screenshots come from the app's built-in browser pane, not agent-browser; ask before installing any headless browser. |
| **no-ai-slop** | https://github.com/petergyang/no-ai-slop | Edits prose to remove AI-writing patterns | Writing only: interface copy, README, model card, reports. Not a design skill. |

All three are installed when task 7.3 starts, not before.

**Not used (do not suggest again):** get-shit-done / GSD (conflicts with our task workflow, token-heavy),
codebase-memory-mcp (duplicates Graphify), Supabase (no database), Strix (login-free localhost app, requires a paid
LLM), agent-reach (social media access not needed), Playwright CLI and Vercel agent-browser (the app's built-in
browser pane is enough), Taste Skill and awesome-design-md (conflict with UI UX Pro Max), Vercel deploy/CLI/autoship
(no Vercel deployment), Next.js, React Native and Turborepo skills (the web app is React + Vite, single package).

## Graphify — MANDATORY RULE

1. Graphify is installed **at the earliest possible moment**: in Task 0.1, as soon as `pyproject.toml`, the package
   skeleton and `.gitignore` exist, within the same task. It is never postponed to a later task.
2. Installation order:
   ```
   uv tool install graphifyy
   graphify install --project          # if not detected: graphify install --project --platform windows
   graphify claude install             # hook that nudges the assistant to the graph before searching
   graphify hook install               # automatic rebuild on git commit / checkout
   ```
   Then build the first graph with `graphify update .` (AST only, no LLM cost). `/graphify .` inside Claude Code adds
   a semantic pass over docs and costs model tokens; use it only when explicitly needed.
3. `data/`, `artifacts/`, `logs/`, `web/node_modules/` and `graphify-out/` are excluded from the graph
   (`.gitignore` is honored automatically; use `.graphifyignore` if needed). The 100+ GB data must never be indexed.
4. Once installed, every question about the code goes to `graphify query` / `graphify path` / `graphify explain`
   **first**, followed only by targeted reads of the needed file.
5. Run `graphify update .` at the end of every task (in case the git hook missed changes).
6. `graphify-out/` is never committed.

## Install commands (for reinstalling)

```
# ponytail and UI UX Pro Max are enabled in the project .claude/settings.json. If needed (two separate messages):
/plugin marketplace add DietrichGebert/ponytail
/plugin install ponytail@ponytail
/plugin marketplace add nextlevelbuilder/ui-ux-pro-max-skill
/plugin install ui-ux-pro-max@ui-ux-pro-max-skill

# Vercel skills (terminal, project folder; --copy: a copy, not a link)
npx skills add vercel-labs/agent-skills --skill vercel-react-best-practices --agent claude-code --copy -y
npx skills add vercel-labs/agent-skills --skill web-design-guidelines --agent claude-code --copy -y
npx skills add vercel-labs/agent-skills --skill vercel-composition-patterns --agent claude-code --copy -y
npx skills add vercel/components.build --skill building-components --agent claude-code --copy -y
npx skills add vercel-labs/before-and-after --skill before-and-after --agent claude-code --copy -y
npx skills add petergyang/no-ai-slop --skill no-ai-slop --agent claude-code --copy -y

# Context7: .mcp.json in the project root (https://mcp.context7.com/mcp)
```

> **Folder name warning:** the folder stays `HealthPproject`. If it is renamed, project-scoped plugin and skill
> registrations stay at the old path and silently stop loading (this happened in an earlier project); reinstall them.

---

# 7. Data Rules

- Folders: `data/raw` (downloaded, never modified) → `data/interim` → `data/processed` (model inputs).
- Every source's URL, version/release date and checksum (where available) is pinned in `config/sources.toml`.
  There is no "download the latest" logic; changing a version is a deliberate decision.
- Scripts must be **idempotent and resumable**: skip when the output exists (unless `--force`), resume interrupted
  downloads.
- Every processing step must be runnable on a small sample (`--chrom 22`, `--limit N`). It is **always** validated on
  a small sample before a full run.
- Very large raw files (gnomAD) are processed chromosome by chromosome; raw files are deleted only with the user's
  approval and an explicit flag (`--delete-raw`).
- Parquet files are partitioned by chromosome and sorted by position.
- Coordinates are 1-based (VCF style). REF is checked against the reference FASTA; mismatching records are dropped
  and counted in the summary. Indels are stored left-normalized.
- **ClinVar labeling rules:**
  - Germline classification is used.
  - Pathogenic + Likely pathogenic → `1`; Benign + Likely benign → `0`.
  - VUS, conflicting interpretations, risk factor, drug response, other, etc. → excluded from training.
    VUS records are kept in a separate table for the interface demo.
  - Review status of at least **1 star** ("criteria provided" or higher).
  - Records without GRCh38 coordinates are dropped; deduplicate on `(chrom, pos, ref, alt)`.
- Raw data, processed data, model weights and ESM outputs are **never** committed.
- The README lists the license and attribution for every data source.

---

# 8. Data Leakage and Evaluation Rules

In this project the most expensive mistake is not a code bug but **unnoticed data leakage.** This section is not negotiable.

- **Gene-disjoint split:** train / validation / test are **disjoint by gene**; all variants of a gene fall into one
  set. A fixed seed is used. The split file is stored under `data/processed/splits/` and its summary (gene and
  variant counts, file hash) is written to `reports/`. The split is never silently regenerated.
- **Phase 2 types** reuse the gene split fixed in Phase 1; a gene's set never changes with variant type.
- **Temporal test:** variants added after the pinned training ClinVar release form a separate temporal test set;
  no training variant appears in it.
- **The test set is locked:** it is opened only in "Final evaluation" tasks. Hyperparameters, thresholds,
  calibration, feature selection and model selection use **the validation set only**. If the model changes after the
  test set has been opened, the report says so explicitly.
- Every label-derived statistic (e.g. "pathogenic fraction in gene", target encoding) is computed **from the
  training set only**. Scalers, PCA and similar transforms are fit on the training set only.
- The same variant cannot appear in two sets. This and gene disjointness are verified by **automated tests**.
- **Metrics:** ROC-AUC, PR-AUC, MCC, balanced accuracy, Brier score, ECE (calibration); confidence intervals via
  gene-level bootstrap. Results are reported **separately per variant type**; a single headline number mixing types
  is never reported on its own.
- **Mandatory ablations** (in final reports): without allele frequency, without ESM, LightGBM vs dual tower.
- Poor results are reported honestly.
- Every training / evaluation summary records seeds, data versions and library versions.

## Leakage Auditor (Subagent)

The **leakage-auditor** in `.claude/agents/leakage-auditor.md` does not run continuously. It is invoked with the
Agent tool only at these moments, reviews only the relevant code and summaries, reports, and exits:

- When labeling rules or split code are written or changed
- When a new feature or data source is added to the feature pipeline
- In every "Final evaluation" task, **before the test set is opened**

When invoking it, name the files and summary files to review; do not have it scan the whole project.
🔴 Critical findings are fixed before moving on. 🟡 and 🟢 findings are reported briefly to the user.

## Code Review

Run the built-in `/code-review` skill at the end of tasks where a bug is expensive downstream: **1.3** and **8.1**
(our annotator), **5.1** (dual tower model) and **7.2** (API). Review only the task's diff. Fix confirmed findings
before the completion report; mention the rest briefly. Do not run it on other tasks unless the user asks.

---

# 9. Context and Token Efficiency

The user's explicit priority is **saving tokens.**

- Once Graphify is installed, query the graph first to explore the code; keep file reads targeted.
- Do not scan the whole repository by default. Do not read unrelated modules "for general familiarity".
- Do not reread large files already read in this session unless they changed.
- **Never read large data files in full.** Schema, row count and the first few rows (`head`, `schema`, `count`
  with polars) are enough.
- **Long-running jobs** (downloads, ESM computation, training):
  - Run in the background, or the user runs them in their own terminal and reports the result.
  - Logs are written to files under `logs/`; progress is not streamed to the terminal (except rare one-line status).
  - Every script ends by writing a **small summary JSON** (at most ~50 lines): counts, durations, metrics, warnings.
    The assistant reads **only this summary**, not the log. Logs are read only for diagnosis, and only in a targeted
    way (`tail`, `grep`).
- Prefer: `git diff --stat`, `git diff -- <path>`, targeted `pytest -k`, targeted search, direct file reads.
- Avoid: repeated full test runs, dumping large generated files, reading inside `node_modules/`, `.venv/`, `data/`,
  `artifacts/`, `.git/`.
- Do not narrate every tool call.

---

# 10. Task Discipline

Implement only the **current task** in `TASKS.md` unless the user explicitly changes the scope.

Before starting:
1. Restate the task in at most five lines.
2. Name the owning module.
3. List the files likely to change.
4. Name relevant existing tests and patterns.
5. Report any direct conflict with the documented architecture.

- Do not write a long plan for a small task.
- Do not add features from later tasks because they seem useful.
- Do not combine unrelated cleanup, refactoring or feature work with the current task. Do not touch unrelated files.
- Do not rename public functions, Parquet columns, modules or folders unless required.
- When a requirement is slightly ambiguous, follow existing patterns and choose the smallest reversible solution
  consistent with the documentation.
- Ask the user **only** when the ambiguity could materially affect: data sources / the forbidden-score boundary,
  label definition, splitting and leakage, use of the test set, meaning of the model output, the API contract,
  disk usage or deletions.
- If a task turns out larger than expected, split it into subtasks after asking the user.
- Do not move to the next task until the user confirms the task's **Try it** step.

---

# 11. Change Strategy

- Small, reviewable patches; no broad rewrites.
- Edit existing files in place; do not regenerate a whole file when a focused edit is enough.
- Do not change working code purely for style.
- Reuse existing repository patterns before creating new abstractions.
- Do not introduce, without a concrete need: new layers, base classes, design patterns, frameworks, services,
  databases, queues, experiment-tracking platforms.
- CLI commands and API endpoints stay thin; logic lives in the owning module.
- The API uses explicit request and response schemas (pydantic); internal data structures are never exposed directly.
- **All code is English:** identifiers, strings, log messages, commit messages.
- **No comments and no docstrings in code.** Code must explain itself through clear naming, small functions and tests.
- The web interface is English by default and may offer a Turkish language option; this is the only place Turkish
  appears in the repository besides the Turkish section of the README.

---

# 12. Testing Rules

Every behavioral change adds or updates appropriate tests.

Priority tests:
- **Our annotator:** unit tests with known variants (both strands, codon boundaries, first and last exons) and bulk
  agreement against ClinVar HGVS.
- **Leakage invariants:** gene overlap between sets = 0, variant overlap = 0, no fitting on non-training data.
- **Labeling rules:** which ClinVar record maps to which label.
- **Feature lookup:** expected frequency / conservation values for known positions.
- **ONNX parity:** PyTorch and ONNX outputs match within tolerance.
- **API:** schema validation, invalid input, file size limits.
- **Frontend:** `tsc --noEmit`, `vite build`, a few Vitest tests for critical components.

- Tests use only small sample data under `tests/fixtures/`; they never depend on the full datasets.
- Do not write tests that restate implementation details; test externally observable behavior and invariants.
- For bug fixes: first a test that reproduces the bug → confirm it fails for the expected reason → smallest fix →
  run the relevant test scope.
- Do not run the full suite after every trivial edit when a targeted run is sufficient.
- Before declaring completion, run the smallest build and test commands that validate the change.

---

# 13. Security and Privacy

- `.env`, `.env.*`, credential files and private keys are never read, displayed or committed.
  No secrets are expected in this project; if ever needed, add a `.env.example` with empty values.
- FastAPI binds to `127.0.0.1` only; CORS allows only the local Vite origin.
- Validate all external input: variant text format, chromosome/position range, REF/ALT letters, uploaded VCF size
  and row limits.
- Uploaded VCF data is processed in memory only; never persisted to disk, never logged.
- No executables are downloaded from untrusted sources; data is downloaded only from the official URLs in
  `config/sources.toml`.
- Deletions (raw data, old outputs) are never performed without the user's approval.

---

# 14. GitHub

- `.gitignore` is created first: `data/`, `artifacts/`, `logs/`, `graphify-out/`, `.venv/`, `__pycache__/`,
  `web/node_modules/`, `web/dist/`, `.env`, `.env.*` (except `.env.example`), `.claude/settings.local.json`.
- The repository is created on GitHub as **private** under the name `Healthcare`. Going public is the user's decision.
- A **single `README.md`** contains an English section followed by a Turkish section (linked from the top) and is
  kept current as tasks progress: what the project does, results (per-type metrics), architecture, data sources and
  licenses, setup and run steps, screenshots, limitations, disclaimer.
- Someone cloning the repository must be able to download the data and reproduce the model on their own machine by
  following the README.
- Suggest a meaningful English commit message at the end of every task (e.g. `Add MANE-based missense annotator`).
  Never commit or push without the user's approval.

---

# 15. Communication

- Talk to the user in **Turkish**, plainly; be concise and implementation-focused. (Repository content stays English.)
- The user is not a medical professional: explain genetics concepts the first time they come up, simply and with
  analogies (DNA = recipe book, gene = recipe, protein = the machine a recipe builds, variant = a changed letter).
- Do not explain basic Python / ML / React concepts unless asked.
- Give early notice only for: requirement conflicts, forbidden-score risk, data leakage risk, disk / data integrity
  risk, or a test or metric result that changes the approach.
- Do not paste complete files; give file paths, short summaries and focused diffs.

---

# 16. Completion Report

Report only, in Turkish, using these Turkish headings:

## Değişenler
- Files changed
- Data schema / API changes
- Behavior added

## Doğrulama
- Build / test commands run and results
- Key numbers or metrics from the summary JSON, if any

## Kararlar
- Assumptions made
- Durable decisions introduced

## Kalanlar
- Known risks
- Deferred work
- **Dene:** what the user should run and what they should see
- Next task

Then: tick the box in `TASKS.md`, update the "Current task" line, add at most 3 lines to "Completion Notes", run
`graphify update .`, update the README if needed, and suggest a commit message.
Update the "Decisions" section of this file only when a durable architectural or domain decision was made.

---

# 17. Definition of Done

A task is complete only when:

- The task's acceptance criteria are met.
- No forbidden score or source was used.
- The gene split and the test lock are intact.
- Required validation (input, REF agreement, schema) exists.
- Relevant tests pass.
- No unrelated functionality was added.
- Documentation reflects any durable change.
- The report states remaining risks honestly.

---

# 18. Decisions

| Date | Decision | Rationale |
|---|---|---|
| 2026-10-11 | Not competing; the problem is solved as a personal CV project with public data. | The team learned about the competition too late. |
| 2026-10-11 | Precomputed pathogenicity scores are never used (Section 2). | Credibility of the "our own model from scratch" claim. |
| 2026-10-11 | Conservation (phyloP/phastCons), gnomAD frequency and gene constraint are allowed. | Raw statistics, not trained on labels. |
| 2026-10-11 | ESM-2 650M is used, run locally by us. | Masked language model that never saw labels; matches the user's dual-tower experience. |
| 2026-10-11 | Single model, dual tower, tower dropout; Phase 1 missense, Phase 2 other coding types. | Hardest and most valuable part first; one ONNX output. |
| 2026-10-11 | Deep intronic, UTR and regulatory variants are out of scope. | No reliable prediction without the forbidden splice scores. |
| 2026-10-11 | Our own annotator (MANE Select based); no dependency on VEP. | "Fully ours" principle; runs natively on Windows. |
| 2026-10-11 | Gene-disjoint split + temporal test + locked test set; LightGBM baseline. | Honest evaluation; prevents memorization and leakage. |
| 2026-10-11 | Backend Python (FastAPI + ONNX Runtime + DuckDB), frontend React + TypeScript + Vite. | UI decoupled from backend; speed on large data is solved in the query layer. |
| 2026-10-11 | Repository fully English, no code comments/docstrings; single bilingual README; UI may offer Turkish. | User's explicit requirement; international CV presentation. |
| 2026-10-11 | ClinVar: training = October 2025 release (cumulative); temporal test = variants first added between October 2025 and October 2026; reclassified variants reported separately, never mixed into the test set. | Honest "future variants" test; reclassifications are a distinct, interesting finding. |
| 2026-10-11 | `/code-review` at the end of tasks 1.3, 5.1, 7.2, 8.1 only. | Fresh-eyes review where bugs are costly, without per-task token cost. |
| 2026-10-12 | Model features are exactly the `f_`-prefixed columns; IDs, `label` and everything else are never fed to a model. | Whitelist makes ClinVar metadata leakage structurally impossible. |
| 2026-10-12 | gnomAD absence: AF = 0 when the position is covered in exomes, null + `f_gnomad_covered = 0` when not; the no-AF ablation drops every column in `FREQUENCY_FEATURES`. | Separates "not observed" from "not measured" without a hidden frequency back door. |
| 2026-10-12 | Variants with a missense annotation in more than one MANE gene are dropped (no ClinVar-based tie-break); a missense call in one gene plus a non-coding call in an overlapping gene is kept. | Same rule at training and inference; keeps gene-disjoint splits clean. |
| 2026-10-12 | Gene split covers all 19,229 MANE genes; genes with overlapping CDS or a shared labeled position form one group; labeled groups are stratified greedily (variant + pathogenic counts), unlabeled groups get a seeded 70/15/15 draw. Pinned by sha256 of `genes.parquet` and of the training feature table. | Phase 2 never needs new assignments; overlapping genes cannot leak across sets. |
| 2026-10-12 | Test and temporal rows live only in `data/processed/splits/locked/`; models read data through `split.genes.load_split`, which raises `LockedSplitError` for them unless `allow_locked=True` (final evaluation tasks only). | Makes the test lock a code-level guarantee, not a convention. |
| 2026-10-13 | Final missense model = average of LightGBM and dual tower probabilities (validation ROC-AUC 0.9901; best or tied on gnomAD-seen and unseen subsets); frequency-based routing rejected (0.9894). Both models are also reported separately. | Ensemble is never worse and better calibrated than LightGBM alone. |
| 2026-10-13 | Three-way decision on the per-variant calibrated probability: ≥ 0.90 "Likely pathogenic", ≤ 0.10 "Likely benign", otherwise "Uncertain" (ACMG Bayesian posterior). Calibrator (Platt on the ensemble logit) is fit on validation only; method picked by gene-grouped 5-fold out-of-fold Brier. | Set-level 90% precision made the two bands overlap (P ≥ 0.31, B ≤ 0.91), leaving no "Uncertain" band. |
| 2026-10-13 | One ONNX file holds LightGBM (onnxmltools, opset 15) + dual tower (legacy `dynamo=False` exporter, no onnxscript) + mean + Platt; dual tower preprocessing runs inside the graph; input = raw `f_` features in LightGBM column order (NaN = missing) + ESM blocks. ONNX float32 tree comparisons are the model's definition. | Serving needs only onnxruntime; f_delta_* float64 noise splits cannot be reproduced in float32 (max LightGBM diff 0.014, 0 decision changes on validation). |
| 2026-10-16 | Serving reuses the training feature functions (`join_sources` → `esm_features` → `finalize`) on a one-row table. The explanation replaces one feature group at a time with 32 seeded training rows and reports the mean change in calibrated log-odds; the protein embedding is not part of it. | No train/serve feature skew (0 of 100 validation variants differ). Masking with NaN misled, because missingness is itself a learned signal; probability changes vanish near 0 and 1, log-odds do not. |
| 2026-10-16 | Web tasks also use the Vercel skills vercel-composition-patterns, building-components and before-and-after (installed at 7.3). | Fit React + Vite; Next.js, React Native and deploy skills do not. |
| 2026-10-11 | Environment: Windows 11, uv, Docker available / no WSL2, data on A: (990 Pro). | User's hardware. |

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
