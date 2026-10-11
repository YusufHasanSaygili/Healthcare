# Healthcare — Genetic Variant Pathogenicity Prediction

[English](#english) · [Türkçe](#türkçe)

> **Disclaimer:** This project is for research and educational purposes only. It must not be used for clinical
> diagnosis or treatment decisions.

---

## English

A from-scratch AI model that predicts whether germline genetic variants are **pathogenic** or **benign**, with a web
interface that runs locally.

- No precomputed pathogenicity scores (CADD, REVEL, AlphaMissense, SpliceAI, …) are used anywhere. The model learns
  only from raw public data: ClinVar labels, the GRCh38 reference genome, MANE Select transcripts, gnomAD allele
  frequencies, evolutionary conservation, UniProt curated features and ESM-2 protein language model outputs computed
  locally.
- Evaluation uses gene-disjoint splits and a temporal hold-out; results are reported per variant type.

**Status:** early development (Phase 0).

### Requirements

- Windows 11 / Linux, Python 3.12, [uv](https://docs.astral.sh/uv/)
- NVIDIA GPU with a recent driver (developed on an RTX 4070 Ti, 12 GB)

### Setup

```bash
uv sync
uv run healthcare doctor
```

---

## Türkçe

Germline genetik varyantların **patojenik** mi **benign** mi olduğunu tahmin eden, sıfırdan eğitilmiş bir yapay
zeka modeli ve yerelde çalışan bir web arayüzü.

- Hiçbir hazır patojenite skoru (CADD, REVEL, AlphaMissense, SpliceAI, …) kullanılmaz. Model yalnızca ham açık
  verilerden öğrenir: ClinVar etiketleri, GRCh38 referans genomu, MANE Select transkriptleri, gnomAD alel
  frekansları, evrimsel korunmuşluk, UniProt küratörlü özellikleri ve yerelde hesaplanan ESM-2 protein dil modeli
  çıktıları.
- Değerlendirme gen bazında ayrık setler ve zamansal bir test setiyle yapılır; sonuçlar varyant türüne göre ayrı
  raporlanır.

**Durum:** erken geliştirme (Aşama 0).

### Gereksinimler

- Windows 11 / Linux, Python 3.12, [uv](https://docs.astral.sh/uv/)
- Güncel sürücülü NVIDIA GPU (RTX 4070 Ti, 12 GB üzerinde geliştirildi)

### Kurulum

```bash
uv sync
uv run healthcare doctor
```

> **Sorumluluk reddi:** Bu proje yalnızca araştırma ve eğitim amaçlıdır; klinik tanı veya tedavi kararları için
> kullanılamaz.
