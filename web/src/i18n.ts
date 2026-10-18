import type { Decision, EvidenceGroup, VariantType } from "./api/client";
import type { UploadProblem } from "./batch";
import type { LabStatus } from "./components/VitalsMonitor";
import type { Direction, ReasonKey } from "./explain";
import type { Route } from "./route";

export type Messages = {
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
  metrics: {
    set: string;
    type: string;
    variants: string;
    rocAuc: string;
    prAuc: string;
    ece: string;
    uncertain: string;
  };
  ablationsTitle: string;
  ablationsNote: string;
  ablation: string;
  dataTitle: string;
  limitationsTitle: string;
  sourceLanguageNote: string;
  explainTitle: string;
  explainVariant: string;
  explainProtein: string;
  explainIndel: string;
  indelDeleted: string;
  indelInserted: string;
  indelReplaced: string;
  explainType: Record<Exclude<VariantType, "missense">, string>;
  typeNames: Record<VariantType, string>;
  explainWhy: Record<Decision, string>;
  explainToward: Record<Direction, string>;
  explainReasons: Record<ReasonKey, string>;
  explainNoReasons: string;
  geneFunctionTitle: string;
  geneFunctionSource: string;
  aminoAcids: Record<string, string>;
};

const en: Messages = {
  title: "Variant Pathogenicity",
  tagline: "All coding types, model v2 · runs on this computer",
  inputLabel: "Variant",
  inputHint:
    "Genomic position (17:43063931 G>A) or protein change (BRCA1 p.Arg1699Trp). GRCh38; coding SNVs and small indels, VCF style (7:117559592 CTTT>C).",
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
    truncation: "Truncation and NMD",
    splicing: "Splicing",
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
  members: "LightGBM {lightgbm} · dual tower {dualTower} · shown value is calibrated per variant type (LightGBM alone or the average of both)",
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
    "Plain-text VCF (GRCh38), tab-separated CHROM POS ID REF ALT. Up to 2 MB and 500 variants; coding SNVs and small indels.",
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
    "Final model: per variant type, LightGBM alone or the average of LightGBM and the dual tower (chosen on validation), then calibrated per type. Brackets show 95% intervals from a gene-level bootstrap.",
  metrics: {
    set: "Set",
    type: "Type",
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
  explainTitle: "What does this mean?",
  explainVariant:
    "This is a variant: a one-letter change in the DNA. At position {position} of chromosome {chrom}, the reference genome has {ref}, but here there is {alt}.",
  explainProtein:
    "It changes part {residue} of the protein made from the {gene} gene: {aaRef} should be there, but {aaAlt} is made instead.",
  explainIndel:
    "This is a variant: a change of letters in the DNA. At position {position} of chromosome {chrom}, {change}.",
  indelDeleted: "{count} letter(s) are missing",
  indelInserted: "{count} extra letter(s) are inserted",
  indelReplaced: "{deleted} letter(s) are replaced by {inserted} other letter(s)",
  explainType: {
    synonymous:
      "The letter changes, but the amino acid it codes for stays the same ({hgvs}). The protein itself is not changed; any harm would come from how the recipe is copied or cut.",
    stop_gained:
      "It turns part {residue} of the {gene} protein into a stop signal, so the protein is cut short ({hgvs}).",
    stop_lost:
      "It removes the stop signal at the end of the {gene} protein, so the protein runs on too long.",
    start_lost:
      "It breaks the start signal of the {gene} protein, so the protein may not be made from the right place.",
    canonical_splice:
      "It hits the cut-and-join signal at an exon edge of {gene}. Part of the recipe may be skipped or kept by mistake.",
    frameshift:
      "The recipe is read three letters at a time. Adding or removing letters that are not a multiple of three shifts the reading from part {residue} of the {gene} protein onward, and the protein usually ends early ({hgvs}).",
    inframe_deletion:
      "It removes a few amino acids from the {gene} protein ({hgvs}); the rest of the protein stays intact.",
    inframe_insertion:
      "It adds a few amino acids to the {gene} protein ({hgvs}); the rest of the protein stays intact.",
    inframe_indel:
      "It replaces a few amino acids of the {gene} protein ({hgvs}); the rest of the protein stays intact.",
  },
  typeNames: {
    missense: "Missense",
    synonymous: "Synonymous",
    stop_gained: "Nonsense (stop gained)",
    stop_lost: "Stop lost",
    start_lost: "Start lost",
    canonical_splice: "Canonical splice site",
    frameshift: "Frameshift",
    inframe_deletion: "In-frame deletion",
    inframe_insertion: "In-frame insertion",
    inframe_indel: "In-frame deletion-insertion",
  },
  explainWhy: {
    likely_pathogenic: "Why it may be harmful",
    uncertain: "Why the result is uncertain",
    likely_benign: "Why it is probably harmless",
  },
  explainToward: { pathogenic: "Points toward harmful", benign: "Points toward harmless" },
  explainReasons: {
    conserved:
      "This position has barely changed across species over millions of years. Places that never change are usually important.",
    notConserved: "This position changes freely between species, so it is probably not critical.",
    absent:
      "This change has never been seen in gnomAD, a database of hundreds of thousands of people.",
    common: "This change is seen in healthy populations ({frequency} in gnomAD). Changes that common are rarely harmful.",
    esmUnusual:
      "The protein language model (ESM-2) finds the new amino acid very unusual at this spot.",
    esmTolerated: "The protein language model (ESM-2) finds the new amino acid acceptable at this spot.",
    chemDifferent:
      "The new amino acid is chemically very different from the original: size, charge and affinity for water change.",
    chemSimilar: "The new amino acid is chemically similar to the original.",
    inRegion: "It falls in a known working part of the protein: {regions}.",
    nmdTriggered:
      "The early stop is far from the end of the message, so the cell will probably destroy it (nonsense-mediated decay) and make no protein.",
    nmdEscapes:
      "The early stop is near the end of the message, so it probably escapes destruction and a shortened protein is made.",
    mostlyLost: "More than half of the protein is lost.",
    mostlyKept: "Only the very end of the protein is affected.",
  },
  explainNoReasons:
    "No single measurement stands out; the result comes from all measurements together.",
  geneFunctionTitle: "What the {gene} gene does",
  geneFunctionSource: "Source: UniProtKB/Swiss-Prot.",
  aminoAcids: {
    Ala: "alanine",
    Arg: "arginine",
    Asn: "asparagine",
    Asp: "aspartic acid",
    Cys: "cysteine",
    Gln: "glutamine",
    Glu: "glutamic acid",
    Gly: "glycine",
    His: "histidine",
    Ile: "isoleucine",
    Leu: "leucine",
    Lys: "lysine",
    Met: "methionine",
    Phe: "phenylalanine",
    Pro: "proline",
    Ser: "serine",
    Thr: "threonine",
    Trp: "tryptophan",
    Tyr: "tyrosine",
    Val: "valine",
    Sec: "selenocysteine",
  },
};

const tr: Messages = {
  title: "Varyant Patojenitesi",
  tagline: "Tüm kodlayan türler, model v2 · bu bilgisayarda çalışır",
  inputLabel: "Varyant",
  inputHint:
    "Genomik konum (17:43063931 G>A) ya da protein değişimi (BRCA1 p.Arg1699Trp). GRCh38; kodlayan SNV'ler ve küçük indel'ler, VCF biçiminde (7:117559592 CTTT>C).",
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
    truncation: "Kesilme ve NMD",
    splicing: "Kırpılma (splicing)",
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
    "LightGBM {lightgbm} · çift kuleli ağ {dualTower} · gösterilen değer varyant türüne göre kalibre edilir (yalnız LightGBM ya da ikisinin ortalaması)",
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
    "Düz metin VCF (GRCh38), sekmeyle ayrılmış CHROM POS ID REF ALT. En fazla 2 MB ve 500 varyant; kodlayan SNV'ler ve küçük indel'ler.",
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
    "Nihai model: her varyant türü için doğrulama kümesinde seçilen yalnız LightGBM ya da LightGBM ile çift kuleli modelin ortalaması, ardından türe göre kalibrasyon. Köşeli parantezler gen düzeyinde bootstrap ile %95 aralıklardır.",
  metrics: {
    set: "Küme",
    type: "Tür",
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
  explainTitle: "Bu ne demek?",
  explainVariant:
    "Bu bir varyant: DNA'da tek harflik bir değişiklik. {chrom}. kromozomun {position}. harfinde referans genomda {ref} varken burada {alt} var.",
  explainProtein:
    "{gene} geninden yapılan proteinin {residue}. parçasını değiştiriyor: orada {aaRef} olması gerekirken {aaAlt} yapılıyor.",
  explainIndel:
    "Bu bir varyant: DNA'da harf değişikliği. {chrom}. kromozomun {position}. harfinde {change}.",
  indelDeleted: "{count} harf eksik",
  indelInserted: "{count} fazladan harf eklenmiş",
  indelReplaced: "{deleted} harf, {inserted} başka harfle değiştirilmiş",
  explainType: {
    synonymous:
      "Harf değişiyor ama kodladığı aminoasit aynı kalıyor ({hgvs}). Proteinin kendisi değişmiyor; bir zarar varsa tarifin kopyalanması ya da kesilip birleştirilmesinden gelir.",
    stop_gained:
      "{gene} proteininin {residue}. parçasını bir dur işaretine çeviriyor; protein yarıda kesiliyor ({hgvs}).",
    stop_lost:
      "{gene} proteininin sonundaki dur işaretini siliyor; protein olması gerekenden uzun devam ediyor.",
    start_lost:
      "{gene} proteininin başlangıç işaretini bozuyor; protein doğru yerden üretilemeyebilir.",
    canonical_splice:
      "{gene} geninde bir exon kenarındaki kes-birleştir işaretine denk geliyor. Tarifin bir parçası yanlışlıkla atlanabilir ya da içeride kalabilir.",
    frameshift:
      "Tarif üçer harf okunur. Üçün katı olmayan sayıda harf eklenip silinince {gene} proteininin {residue}. parçasından itibaren okuma kayar ve protein genelde erken biter ({hgvs}).",
    inframe_deletion:
      "{gene} proteininden birkaç aminoasit çıkarıyor ({hgvs}); proteinin geri kalanı sağlam kalıyor.",
    inframe_insertion:
      "{gene} proteinine birkaç aminoasit ekliyor ({hgvs}); proteinin geri kalanı sağlam kalıyor.",
    inframe_indel:
      "{gene} proteininin birkaç aminoasidini değiştiriyor ({hgvs}); proteinin geri kalanı sağlam kalıyor.",
  },
  typeNames: {
    missense: "Missense",
    synonymous: "Synonymous (sessiz)",
    stop_gained: "Nonsense (erken dur)",
    stop_lost: "Dur işareti kaybı",
    start_lost: "Başlangıç kaybı",
    canonical_splice: "Kanonik splice bölgesi",
    frameshift: "Frameshift (çerçeve kayması)",
    inframe_deletion: "Çerçeve içi silme",
    inframe_insertion: "Çerçeve içi ekleme",
    inframe_indel: "Çerçeve içi silme-ekleme",
  },
  explainWhy: {
    likely_pathogenic: "Neden zararlı olabilir?",
    uncertain: "Sonuç neden belirsiz?",
    likely_benign: "Neden muhtemelen zararsız?",
  },
  explainToward: { pathogenic: "Zararlı yönü gösterenler", benign: "Zararsız yönü gösterenler" },
  explainReasons: {
    conserved:
      "Bu konum türler arasında milyonlarca yıldır neredeyse hiç değişmemiş. Hiç değişmeyen yerler genelde önemlidir.",
    notConserved: "Bu konum türler arasında serbestçe değişiyor; muhtemelen kritik bir yer değil.",
    absent: "Bu değişiklik yüz binlerce kişiyi içeren gnomAD veritabanında hiç görülmemiş.",
    common:
      "Bu değişiklik sağlıklı topluluklarda görülüyor (gnomAD'de {frequency}). Bu kadar yaygın değişiklikler nadiren zararlıdır.",
    esmUnusual: "Protein dil modeli (ESM-2) yeni aminoasidi bu noktada çok alışılmadık buluyor.",
    esmTolerated: "Protein dil modeli (ESM-2) yeni aminoasidi bu noktada kabul edilebilir buluyor.",
    chemDifferent:
      "Yeni aminoasit kimyasal olarak orijinalinden çok farklı: boyutu, yükü ve suya yatkınlığı değişiyor.",
    chemSimilar: "Yeni aminoasit kimyasal olarak orijinaline benziyor.",
    inRegion: "Proteinin bilinen çalışan bir bölgesine denk geliyor: {regions}.",
    nmdTriggered:
      "Erken dur işareti mesajın sonundan uzakta; hücre büyük ihtimalle bu mesajı imha eder (NMD) ve protein hiç üretilmez.",
    nmdEscapes:
      "Erken dur işareti mesajın sonuna yakın; mesaj büyük ihtimalle imhadan kaçar ve kısalmış bir protein üretilir.",
    mostlyLost: "Proteinin yarısından fazlası kayboluyor.",
    mostlyKept: "Sadece proteinin en son kısmı etkileniyor.",
  },
  explainNoReasons:
    "Tek başına öne çıkan bir ölçüm yok; sonuç tüm ölçümlerin birleşiminden geliyor.",
  geneFunctionTitle: "{gene} geni ne iş yapar?",
  geneFunctionSource: "Kaynak: UniProtKB/Swiss-Prot (metin İngilizcedir).",
  aminoAcids: {
    Ala: "alanin",
    Arg: "arginin",
    Asn: "asparagin",
    Asp: "aspartik asit",
    Cys: "sistein",
    Gln: "glutamin",
    Glu: "glutamik asit",
    Gly: "glisin",
    His: "histidin",
    Ile: "izolösin",
    Leu: "lösin",
    Lys: "lizin",
    Met: "metiyonin",
    Phe: "fenilalanin",
    Pro: "prolin",
    Ser: "serin",
    Thr: "treonin",
    Trp: "triptofan",
    Tyr: "tirozin",
    Val: "valin",
    Sec: "selenosistein",
  },
};

export const MESSAGES = { en, tr } as const;

export function fill(template: string, values: Record<string, string>): string {
  return template.replace(/\{(\w+)\}/g, (match, key: string) => values[key] ?? match);
}
