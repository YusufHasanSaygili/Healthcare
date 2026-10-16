import type { Decision, EvidenceGroup } from "./api/client";
import type { UploadProblem } from "./batch";
import type { LabStatus } from "./components/VitalsMonitor";
import type { Route } from "./route";

type Messages = {
  title: string;
  tagline: string;
  inputLabel: string;
  inputHint: string;
  submit: string;
  submitting: string;
  examplesLabel: string;
  themeToLight: string;
  themeToDark: string;
  languageLabel: string;
  awaitingTitle: string;
  awaitingHint: string;
  specimen: string;
  probabilityLabel: string;
  decisions: Record<Decision, string>;
  gaugeLabel: string;
  gaugeBenign: string;
  gaugeUncertain: string;
  gaugePathogenic: string;
  whyTitle: string;
  whyNote: string;
  strongestOneWay: string;
  strongestMixed: string;
  towardPathogenic: string;
  towardBenign: string;
  and: string;
  groups: Record<EvidenceGroup, string>;
  measurementsTitle: string;
  frequencyLabel: string;
  frequencyNotObserved: string;
  frequencyNotMeasured: string;
  frequencyOneIn: string;
  conservationLabel: string;
  conservationHint: string;
  esmLabel: string;
  esmHint: string;
  granthamLabel: string;
  granthamHint: string;
  regionsLabel: string;
  regionsNone: string;
  regions: Record<string, string>;
  notAvailable: string;
  members: string;
  outOfScopeTitle: string;
  consequences: Record<string, string>;
  invalidInputTitle: string;
  networkError: string;
  unexpectedError: string;
  disclaimer: string;
  footer: string;
  monitor: Record<LabStatus, string>;
  navLabel: string;
  nav: Record<Route, string>;
  batchTitle: string;
  batchHint: string;
  batchFileLabel: string;
  batchSubmit: string;
  batchSubmitting: string;
  batchSample: string;
  batchPrivacy: string;
  batchRejectedTitle: string;
  uploadProblems: Record<UploadProblem, string>;
  batchSummary: string;
  downloadCsv: string;
  columns: { line: string; variant: string; gene: string; protein: string; probability: string; decision: string };
  outOfScope: string;
  lineErrorsTitle: string;
  lineErrorsHidden: string;
  noResults: string;
  modelTitle: string;
  modelLoading: string;
  modelFile: string;
  variantTypes: string;
  thresholdsTitle: string;
  thresholdsText: string;
  evaluationTitle: string;
  evaluationNote: string;
  metrics: { set: string; variants: string; rocAuc: string; prAuc: string; ece: string; uncertain: string };
  ablationsTitle: string;
  ablationsNote: string;
  ablation: string;
  dataTitle: string;
  limitationsTitle: string;
  sourceLanguageNote: string;
};

const en: Messages = {
  title: "Variant Pathogenicity",
  tagline: "Missense model v1 · runs on this computer",
  inputLabel: "Variant",
  inputHint:
    "Genomic position (17:43063931 G>A) or protein change (BRCA1 p.Arg1699Trp). GRCh38, missense SNVs.",
  submit: "Analyze",
  submitting: "Analyzing…",
  examplesLabel: "Examples",
  themeToLight: "Switch to light theme",
  themeToDark: "Switch to dark theme",
  languageLabel: "Language",
  awaitingTitle: "Awaiting specimen",
  awaitingHint: "Enter a variant or pick an example. Nothing leaves this computer.",
  specimen: "Specimen",
  probabilityLabel: "Probability of pathogenicity",
  decisions: {
    likely_pathogenic: "Likely pathogenic",
    uncertain: "Uncertain",
    likely_benign: "Likely benign",
  },
  gaugeLabel: "Probability {value}, in the {band} band",
  gaugeBenign: "Likely benign",
  gaugeUncertain: "Uncertain",
  gaugePathogenic: "Likely pathogenic",
  whyTitle: "What drove the result",
  whyNote:
    "Each bar shows how far the prediction moves when that evidence is replaced by values from typical variants (log-odds).",
  strongestOneWay: "Strongest evidence {direction}: {list}.",
  strongestMixed: "{first} points {firstDirection}; {second} points {secondDirection}.",
  towardPathogenic: "toward pathogenic",
  towardBenign: "toward benign",
  and: "and",
  groups: {
    allele_frequency: "Population frequency",
    conservation: "Evolutionary conservation",
    esm_scores: "Protein language model (ESM-2)",
    gene_constraint: "Gene constraint",
    protein_regions: "Protein regions",
    amino_acid_change: "Amino acid change",
    protein_position: "Position in the protein",
  },
  measurementsTitle: "Measurements",
  frequencyLabel: "Population frequency (gnomAD)",
  frequencyNotObserved: "Not seen in gnomAD",
  frequencyNotMeasured: "Not measured",
  frequencyOneIn: "1 in {value} alleles",
  conservationLabel: "Conservation (phyloP)",
  conservationHint: "Higher means more conserved across species",
  esmLabel: "ESM-2 score (LLR)",
  esmHint: "More negative means the protein model finds the change less likely",
  granthamLabel: "Amino acid distance (Grantham)",
  granthamHint: "Higher means the two amino acids differ more",
  regionsLabel: "Protein regions (UniProt)",
  regionsNone: "None annotated",
  regions: {
    domain: "Domain",
    active_site: "Active site",
    binding_site: "Binding site",
    motif: "Motif",
    region: "Region",
    coiled_coil: "Coiled coil",
    transmembrane: "Transmembrane",
    signal_or_propeptide: "Signal or propeptide",
    disulfide: "Disulfide bond",
    ptm_site: "Modification site",
    helix: "Helix",
    strand: "Strand",
    turn: "Turn",
  },
  notAvailable: "Not available",
  members: "LightGBM {lightgbm} · dual tower {dualTower} · shown value is their calibrated average",
  outOfScopeTitle: "Not analyzed",
  consequences: {
    intergenic: "No MANE Select transcript covers this position.",
    multiple_genes: "The change is missense in more than one gene, so it is excluded.",
    multi_nucleotide: "This protein change needs more than one DNA letter to change.",
    synonymous: "Synonymous variants are not supported yet.",
    stop_gained: "Stop-gained variants are not supported yet.",
    stop_lost: "Stop-lost variants are not supported yet.",
    start_lost: "Start-lost variants are not supported yet.",
    canonical_splice: "Splice-site variants are not supported yet.",
    intronic: "Intronic variants are out of scope.",
    utr: "UTR variants are out of scope.",
  },
  invalidInputTitle: "Check the variant",
  networkError: "The local API is not reachable. Start it with: uv run healthcare serve",
  unexpectedError: "Something went wrong. Try again.",
  disclaimer:
    "For research and educational use only. Not for clinical diagnosis or treatment decisions.",
  footer: "GRCh38 · MANE Select · ClinVar labels · tested on genes the model never saw",
  monitor: {
    idle: "NO SPECIMEN",
    analyzing: "ANALYZING",
    likely_pathogenic: "PATHOGENIC SIGNAL",
    uncertain: "UNCERTAIN SIGNAL",
    likely_benign: "BENIGN SIGNAL",
    out_of_scope: "NO READING",
    error: "SIGNAL LOST",
  },
  navLabel: "Pages",
  nav: { single: "Single variant", batch: "Batch (VCF)", model: "Model card" },
  batchTitle: "Batch query",
  batchHint:
    "Plain-text VCF (GRCh38), tab-separated CHROM POS ID REF ALT. Up to 2 MB and 500 variants; SNVs only.",
  batchFileLabel: "VCF file",
  batchSubmit: "Analyze file",
  batchSubmitting: "Analyzing…",
  batchSample: "Download a sample VCF",
  batchPrivacy: "The file is read in memory by the local API and never saved.",
  batchRejectedTitle: "The file was rejected",
  uploadProblems: {
    empty: "Choose a VCF file first; the selected file is empty.",
    too_large: "The file is larger than 2 MB. Split it into smaller files.",
    compressed: "Compressed files are not supported. Unzip it and upload the plain-text .vcf.",
    over_limit: "The file is over the limit: at most 2 MB and 500 variants. Split it into smaller files.",
    not_text: "The file is not plain UTF-8 text. Export it as an uncompressed .vcf.",
  },
  batchSummary:
    "Read: {variants} · Scored: {scored} · Out of scope: {skipped} · Line errors: {errors} · {seconds} s",
  downloadCsv: "Download CSV",
  columns: {
    line: "Line",
    variant: "Variant",
    gene: "Gene",
    protein: "Protein change",
    probability: "Probability",
    decision: "Decision",
  },
  outOfScope: "Out of scope",
  lineErrorsTitle: "Lines that could not be read",
  lineErrorsHidden: "{count} more not shown.",
  noResults: "No variant in the file could be scored.",
  modelTitle: "Model card",
  modelLoading: "Loading the model card…",
  modelFile: "Model file",
  variantTypes: "Variant types",
  thresholdsTitle: "Decision thresholds",
  thresholdsText: "Likely pathogenic at {high} or above, likely benign at {low} or below, uncertain in between.",
  evaluationTitle: "Evaluation",
  evaluationNote:
    "Final model: calibrated average of LightGBM and the dual tower. Brackets show 95% intervals from a gene-level bootstrap.",
  metrics: {
    set: "Set",
    variants: "Variants",
    rocAuc: "ROC-AUC",
    prAuc: "PR-AUC",
    ece: "ECE",
    uncertain: "Uncertain",
  },
  ablationsTitle: "Baseline and ablations",
  ablationsNote: "ROC-AUC on the gene-disjoint test set, to compare with the final model above.",
  ablation: "Model",
  dataTitle: "Data versions",
  limitationsTitle: "Limitations",
  sourceLanguageNote: "",
};

const tr: Messages = {
  title: "Varyant Patojenitesi",
  tagline: "Missense modeli v1 · bu bilgisayarda çalışır",
  inputLabel: "Varyant",
  inputHint:
    "Genomik konum (17:43063931 G>A) ya da protein değişimi (BRCA1 p.Arg1699Trp). GRCh38, missense SNV'ler.",
  submit: "Analiz et",
  submitting: "Analiz ediliyor…",
  examplesLabel: "Örnekler",
  themeToLight: "Açık temaya geç",
  themeToDark: "Koyu temaya geç",
  languageLabel: "Dil",
  awaitingTitle: "Numune bekleniyor",
  awaitingHint: "Bir varyant girin ya da bir örnek seçin. Hiçbir veri bu bilgisayardan çıkmaz.",
  specimen: "Numune",
  probabilityLabel: "Patojenik olma olasılığı",
  decisions: {
    likely_pathogenic: "Muhtemelen patojenik",
    uncertain: "Belirsiz",
    likely_benign: "Muhtemelen benign",
  },
  gaugeLabel: "Olasılık {value}, {band} aralığında",
  gaugeBenign: "Muhtemelen benign",
  gaugeUncertain: "Belirsiz",
  gaugePathogenic: "Muhtemelen patojenik",
  whyTitle: "Sonucu ne belirledi",
  whyNote:
    "Her çubuk, o kanıt tipik varyantların değerleriyle değiştirildiğinde tahminin ne kadar değiştiğini gösterir (log-oran).",
  strongestOneWay: "En güçlü kanıtlar {direction}: {list}.",
  strongestMixed: "{first} {firstDirection}, {second} ise {secondDirection} işaret ediyor.",
  towardPathogenic: "patojenik yönde",
  towardBenign: "benign yönde",
  and: "ve",
  groups: {
    allele_frequency: "Popülasyon sıklığı",
    conservation: "Evrimsel korunmuşluk",
    esm_scores: "Protein dil modeli (ESM-2)",
    gene_constraint: "Gen kısıtlaması",
    protein_regions: "Protein bölgeleri",
    amino_acid_change: "Amino asit değişimi",
    protein_position: "Proteindeki konum",
  },
  measurementsTitle: "Ölçümler",
  frequencyLabel: "Popülasyon sıklığı (gnomAD)",
  frequencyNotObserved: "gnomAD'de görülmedi",
  frequencyNotMeasured: "Ölçülmedi",
  frequencyOneIn: "{value} alelde 1",
  conservationLabel: "Korunmuşluk (phyloP)",
  conservationHint: "Yüksek değer, türler arasında daha korunmuş demektir",
  esmLabel: "ESM-2 skoru (LLR)",
  esmHint: "Daha negatif değer, protein modelinin bu değişimi daha az beklediğini gösterir",
  granthamLabel: "Amino asit farkı (Grantham)",
  granthamHint: "Yüksek değer, iki amino asidin birbirinden daha farklı olduğunu gösterir",
  regionsLabel: "Protein bölgeleri (UniProt)",
  regionsNone: "İşaretli bölge yok",
  regions: {
    domain: "Alan (domain)",
    active_site: "Aktif bölge",
    binding_site: "Bağlanma bölgesi",
    motif: "Motif",
    region: "Bölge",
    coiled_coil: "Sarmal sarmal (coiled coil)",
    transmembrane: "Zar geçişi",
    signal_or_propeptide: "Sinyal ya da propeptit",
    disulfide: "Disülfit bağı",
    ptm_site: "Modifikasyon bölgesi",
    helix: "Heliks",
    strand: "Beta zinciri",
    turn: "Dönüş",
  },
  notAvailable: "Yok",
  members:
    "LightGBM {lightgbm} · çift kuleli ağ {dualTower} · gösterilen değer ikisinin kalibre ortalamasıdır",
  outOfScopeTitle: "Analiz edilmedi",
  consequences: {
    intergenic: "Bu konumu kapsayan bir MANE Select transkripti yok.",
    multiple_genes: "Değişim birden fazla gende missense olduğu için kapsam dışı.",
    multi_nucleotide: "Bu protein değişimi için birden fazla DNA harfinin değişmesi gerekir.",
    synonymous: "Sinonim varyantlar henüz desteklenmiyor.",
    stop_gained: "Erken durdurma (stop-gained) varyantları henüz desteklenmiyor.",
    stop_lost: "Durdurma kaybı (stop-lost) varyantları henüz desteklenmiyor.",
    start_lost: "Başlangıç kaybı (start-lost) varyantları henüz desteklenmiyor.",
    canonical_splice: "Kırpılma bölgesi (splice) varyantları henüz desteklenmiyor.",
    intronic: "İntronik varyantlar kapsam dışı.",
    utr: "UTR varyantları kapsam dışı.",
  },
  invalidInputTitle: "Varyantı kontrol edin",
  networkError: "Yerel API'ye ulaşılamıyor. Şu komutla başlatın: uv run healthcare serve",
  unexpectedError: "Bir sorun oluştu. Tekrar deneyin.",
  disclaimer:
    "Yalnızca araştırma ve eğitim amaçlıdır. Klinik tanı veya tedavi kararlarında kullanılamaz.",
  footer: "GRCh38 · MANE Select · ClinVar etiketleri · modelin hiç görmediği genlerde test edildi",
  monitor: {
    idle: "NUMUNE YOK",
    analyzing: "ANALİZ",
    likely_pathogenic: "PATOJENİK SİNYAL",
    uncertain: "BELİRSİZ SİNYAL",
    likely_benign: "BENİGN SİNYAL",
    out_of_scope: "OKUMA YOK",
    error: "SİNYAL KAYBI",
  },
  navLabel: "Sayfalar",
  nav: { single: "Tek varyant", batch: "Toplu (VCF)", model: "Model kartı" },
  batchTitle: "Toplu sorgu",
  batchHint:
    "Düz metin VCF (GRCh38), sekmeyle ayrılmış CHROM POS ID REF ALT. En fazla 2 MB ve 500 varyant; yalnızca SNV.",
  batchFileLabel: "VCF dosyası",
  batchSubmit: "Dosyayı analiz et",
  batchSubmitting: "Analiz ediliyor…",
  batchSample: "Örnek VCF indir",
  batchPrivacy: "Dosya yerel API tarafından bellekte okunur ve asla kaydedilmez.",
  batchRejectedTitle: "Dosya reddedildi",
  uploadProblems: {
    empty: "Önce bir VCF dosyası seçin; seçilen dosya boş.",
    too_large: "Dosya 2 MB'tan büyük. Daha küçük dosyalara bölün.",
    compressed: "Sıkıştırılmış dosyalar desteklenmiyor. Açıp düz metin .vcf olarak yükleyin.",
    over_limit: "Dosya sınırı aşıyor: en fazla 2 MB ve 500 varyant. Daha küçük dosyalara bölün.",
    not_text: "Dosya düz UTF-8 metin değil. Sıkıştırılmamış .vcf olarak dışa aktarın.",
  },
  batchSummary:
    "Okunan: {variants} · Skorlanan: {scored} · Kapsam dışı: {skipped} · Satır hatası: {errors} · {seconds} sn",
  downloadCsv: "CSV indir",
  columns: {
    line: "Satır",
    variant: "Varyant",
    gene: "Gen",
    protein: "Protein değişimi",
    probability: "Olasılık",
    decision: "Karar",
  },
  outOfScope: "Kapsam dışı",
  lineErrorsTitle: "Okunamayan satırlar",
  lineErrorsHidden: "{count} satır daha gösterilmiyor.",
  noResults: "Dosyadaki hiçbir varyant skorlanamadı.",
  modelTitle: "Model kartı",
  modelLoading: "Model kartı yükleniyor…",
  modelFile: "Model dosyası",
  variantTypes: "Varyant türleri",
  thresholdsTitle: "Karar eşikleri",
  thresholdsText: "{high} ve üstü olası patojenik, {low} ve altı olası benign, arası belirsiz.",
  evaluationTitle: "Değerlendirme",
  evaluationNote:
    "Nihai model: LightGBM ve çift kuleli modelin kalibre edilmiş ortalaması. Köşeli parantezler gen düzeyinde bootstrap ile %95 aralıklardır.",
  metrics: {
    set: "Küme",
    variants: "Varyant",
    rocAuc: "ROC-AUC",
    prAuc: "PR-AUC",
    ece: "ECE",
    uncertain: "Belirsiz",
  },
  ablationsTitle: "Temel model ve ablasyonlar",
  ablationsNote: "Gen-ayrık test kümesinde ROC-AUC; yukarıdaki nihai modelle karşılaştırın.",
  ablation: "Model",
  dataTitle: "Veri sürümleri",
  limitationsTitle: "Sınırlamalar",
  sourceLanguageNote: "Küme adları, veri sürümleri ve sınırlamalar modelin raporundan geldiği için İngilizcedir.",
};

export const MESSAGES = { en, tr } as const;

export function fill(template: string, values: Record<string, string>): string {
  return template.replace(/\{(\w+)\}/g, (match, key: string) => values[key] ?? match);
}
