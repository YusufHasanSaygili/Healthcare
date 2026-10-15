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

**Status:** missense model v1 finished and evaluated; prediction service and web interface in progress.

### Results (missense SNVs, model v1)

Final model: calibrated average of LightGBM and a dual tower network (ESM-2 protein tower + measurement tower),
exported to a single ONNX file. Thresholds and calibration come from the validation set only. Full report:
[`reports/missense_v1.md`](reports/missense_v1.md).

| Evaluation set | Variants | ROC-AUC [95% gene CI] | PR-AUC | Uncertain share |
|---|---|---|---|---|
| Gene-disjoint test | 27,501 | 0.990 [0.988–0.992] | 0.977 | 12.2% |
| Temporal (new in ClinVar after Oct 2025) | 18,205 | 0.989 [0.987–0.990] | 0.980 | 15.1% |
| Temporal, genes unseen in training | 2,919 | 0.988 [0.984–0.991] | 0.973 | 14.4% |
| Test, without allele frequency | 27,501 | 0.969 [0.966–0.972] | 0.926 | – |
| Test, without ESM-2 | 27,501 | 0.985 [0.983–0.987] | 0.964 | – |

Calibration drifts on the newer temporal variants (ECE 0.033 vs 0.004 on test). Other variant types are not
covered yet.

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

**Durum:** missense modeli v1 tamamlandı ve değerlendirildi; tahmin servisi ve web arayüzü geliştiriliyor.

### Sonuçlar (missense SNV'ler, model v1)

Final model: LightGBM ile çift kuleli ağın (ESM-2 protein kulesi + ölçüm kulesi) kalibre edilmiş ortalaması, tek
bir ONNX dosyası. Eşikler ve kalibrasyon yalnızca doğrulama setinden belirlendi. Tam rapor (İngilizce):
[`reports/missense_v1.md`](reports/missense_v1.md).

| Değerlendirme seti | Varyant | ROC-AUC [%95 gen GA] | PR-AUC | Belirsiz payı |
|---|---|---|---|---|
| Gen bazında ayrık test | 27.501 | 0,990 [0,988–0,992] | 0,977 | %12,2 |
| Zamansal (Ekim 2025 sonrası ClinVar'a eklenen) | 18.205 | 0,989 [0,987–0,990] | 0,980 | %15,1 |
| Zamansal, eğitimde görülmemiş genler | 2.919 | 0,988 [0,984–0,991] | 0,973 | %14,4 |
| Test, alel frekansı olmadan | 27.501 | 0,969 [0,966–0,972] | 0,926 | – |
| Test, ESM-2 olmadan | 27.501 | 0,985 [0,983–0,987] | 0,964 | – |

Yeni (zamansal) varyantlarda kalibrasyon kayıyor (ECE 0,033; testte 0,004). Diğer varyant türleri henüz kapsanmıyor.

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
