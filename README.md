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

Calibration drifts on the newer temporal variants (ECE 0.033 vs 0.004 on test).

### Results (all coding types, model v2, gene-disjoint test)

Joint LightGBM, averaged with a joint dual tower for missense, in-frame indels and synonymous (routes and
per-type Platt calibration chosen on validation only). ROC-AUC per variant type; the no-allele-frequency column
matters most for loss-of-function types, whose ClinVar benign labels largely follow population frequency. Full
report: [`reports/all_types_v2.md`](reports/all_types_v2.md).

| Type | Variants (benign) | ROC-AUC [95% gene CI] | Without allele frequency |
|---|---|---|---|
| Missense | 27,500 (19,497) | 0.990 [0.989–0.992] | 0.961 |
| Stop gained | 8,531 (67) | 0.972 [0.954–0.989] | 0.827 |
| Canonical splice | 6,293 (91) | 0.950 [0.926–0.973] | 0.840 |
| Frameshift | 14,176 (138) | 0.902 [0.841–0.959] | 0.721 |
| In-frame deletion | 490 (191) | 0.990 [0.982–0.997] | 0.948 |
| In-frame insertion | 213 (126) | 0.985 [0.967–0.998] | 0.894 |
| Synonymous | 98,509 (81 pathogenic) | 0.889 [0.836–0.936] | 0.831 |

Start loss, stop loss and in-frame delins have too few of one class in the test set to score.

### Requirements

- Windows 11 / Linux, Python 3.12, [uv](https://docs.astral.sh/uv/)
- NVIDIA GPU with a recent driver (developed on an RTX 4070 Ti, 12 GB)

### Setup

```bash
uv sync
uv run healthcare doctor
```

### Predict a variant

After the data is built and the v2 model is exported (`uv run healthcare export-onnx-v2`), any coding SNV or small
indel can be queried by genomic position (VCF style, left-normalized automatically) or a missense change by protein
notation:

```bash
uv run healthcare predict 17:43063931 G>A
uv run healthcare predict 17:43124027 ACT>A
uv run healthcare predict BRCA1 p.Arg1699Trp
```

The output shows the variant type, the calibrated probability, the three-way decision and the evidence groups that
moved it most. Non-coding variants get an explicit out-of-scope message.

The same predictions are available from a local API (bound to `127.0.0.1` only); its interactive documentation
opens at `http://127.0.0.1:8000/docs`:

```bash
uv run healthcare serve
```

### Web interface

With the API running, start the interface in a second terminal and open `http://127.0.0.1:5173`:

```bash
cd web
npm install
npm run dev
```

After changing an API schema, regenerate the frontend types with `uv run healthcare openapi` and then
`npm run gen:api` in `web/`; `npm run build` type-checks and builds the interface, and `npm test` runs its unit
tests. The interface has light and dark themes and an English/Turkish language switch.

Three pages: **Single variant** (genomic or protein change, with the explanation), **Batch (VCF)** (upload a
plain-text GRCh38 VCF of up to 2 MB and 500 SNVs, see the results table and download a CSV; `web/public/sample.vcf`
is a ready-made example) and **Model card** (thresholds, test and temporal metrics, baseline and ablations, data
versions, limitations). Uploaded files are read in memory by the local API and never saved.

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

Yeni (zamansal) varyantlarda kalibrasyon kayıyor (ECE 0,033; testte 0,004).

### Sonuçlar (tüm kodlayan türler, model v2, gen bazında ayrık test)

Ortak LightGBM; missense, in-frame indel ve synonymous için ortak çift kuleli modelle ortalaması alınır (yollar ve
tür bazında Platt kalibrasyonu yalnızca doğrulama setinden). ROC-AUC tür bazındadır; alel frekansı olmadan sütunu
özellikle işlev kaybı türlerinde önemlidir, çünkü bu türlerdeki ClinVar benign etiketleri büyük ölçüde toplumdaki
sıklığa dayanır. Tam rapor (İngilizce): [`reports/all_types_v2.md`](reports/all_types_v2.md).

| Tür | Varyant (benign) | ROC-AUC [%95 gen GA] | Alel frekansı olmadan |
|---|---|---|---|
| Missense | 27.500 (19.497) | 0,990 [0,989–0,992] | 0,961 |
| Stop kazanımı | 8.531 (67) | 0,972 [0,954–0,989] | 0,827 |
| Kanonik splice | 6.293 (91) | 0,950 [0,926–0,973] | 0,840 |
| Frameshift | 14.176 (138) | 0,902 [0,841–0,959] | 0,721 |
| In-frame silme | 490 (191) | 0,990 [0,982–0,997] | 0,948 |
| In-frame ekleme | 213 (126) | 0,985 [0,967–0,998] | 0,894 |
| Synonymous | 98.509 (81 patojenik) | 0,889 [0,836–0,936] | 0,831 |

Start kaybı, stop kaybı ve in-frame delins'te testte bir sınıftan yeterli örnek yok.

### Gereksinimler

- Windows 11 / Linux, Python 3.12, [uv](https://docs.astral.sh/uv/)
- Güncel sürücülü NVIDIA GPU (RTX 4070 Ti, 12 GB üzerinde geliştirildi)

### Kurulum

```bash
uv sync
uv run healthcare doctor
```

### Varyant sorgulama

Veriler hazırlanıp v2 modeli dışa aktarıldıktan sonra (`uv run healthcare export-onnx-v2`) kodlayan her SNV ya da
küçük indel genomik konumla (VCF biçiminde, otomatik sola hizalanır), missense değişimler ise protein gösterimiyle
sorgulanabilir:

```bash
uv run healthcare predict 17:43063931 G>A
uv run healthcare predict 17:43124027 ACT>A
uv run healthcare predict BRCA1 p.Arg1699Trp
```

Çıktıda varyant türü, kalibre olasılık, üç yollu karar ve sonucu en çok etkileyen kanıt grupları yer alır. Kodlamayan
varyantlar için açık bir "kapsam dışı" mesajı verilir.

Aynı tahminler yerel bir API üzerinden de alınabilir (yalnızca `127.0.0.1`); etkileşimli belgeleri
`http://127.0.0.1:8000/docs` adresinde açılır:

```bash
uv run healthcare serve
```

### Web arayüzü

API çalışırken ikinci bir terminalde arayüzü başlatın ve `http://127.0.0.1:5173` adresini açın:

```bash
cd web
npm install
npm run dev
```

API şeması değiştiğinde frontend tiplerini `uv run healthcare openapi` ve ardından `web/` içinde
`npm run gen:api` ile yeniden üretin; `npm run build` tip kontrolü yapar ve arayüzü derler, `npm test` birim
testlerini çalıştırır. Arayüzde açık/koyu tema ve İngilizce/Türkçe dil seçeneği vardır.

Üç sayfa vardır: **Tek varyant** (genomik ya da protein değişimi, açıklamasıyla), **Toplu (VCF)** (en fazla
2 MB ve 500 SNV içeren düz metin GRCh38 VCF yükleyin, sonuç tablosunu görün ve CSV indirin; `web/public/sample.vcf`
hazır bir örnektir) ve **Model kartı** (eşikler, test ve zamansal metrikler, temel model ve ablasyonlar, veri
sürümleri, sınırlamalar). Yüklenen dosyalar yerel API tarafından bellekte okunur ve asla kaydedilmez.

> **Sorumluluk reddi:** Bu proje yalnızca araştırma ve eğitim amaçlıdır; klinik tanı veya tedavi kararları için
> kullanılamaz.
