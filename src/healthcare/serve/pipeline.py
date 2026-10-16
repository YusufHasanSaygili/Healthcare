import json
import statistics
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import onnxruntime as ort
import polars as pl

from healthcare.annotate.snv import (
    SNV_PATTERN,
    Annotation,
    TranscriptIndex,
    parse_protein_change,
    parse_snv,
    protein_change_snvs,
)
from healthcare.data.annotations import load_transcripts
from healthcare.esm import store, variant_inputs
from healthcare.features import missense
from healthcare.features.esm_features import ESM_FEATURES, esm_features

FEATURE_GROUPS = {
    "allele_frequency": tuple(missense.FREQUENCY_FEATURES),
    "conservation": ("f_phylop", "f_phastcons"),
    "esm_scores": tuple(ESM_FEATURES),
    "gene_constraint": tuple(missense.CONSTRAINT_FEATURES.values()),
    "protein_regions": (*missense.UNIPROT_GROUPS, "f_uniprot_mapped"),
    "amino_acid_change": (
        "f_aa_ref",
        "f_aa_alt",
        "f_grantham",
        "f_blosum62",
        "f_delta_charge",
        "f_delta_hydropathy",
        "f_delta_polarity",
        "f_delta_volume",
    ),
    "protein_position": ("f_relative_position", "f_protein_length"),
}
PHASE_TWO = {"synonymous", "stop_gained", "stop_lost", "start_lost", "canonical_splice"}
NON_CODING = {"intronic", "utr"}
OUTPUTS = ["probability", "lightgbm", "dual_tower"]
BACKGROUND_SEED = 20261016
Variant = tuple[str, int, str, str]
MAX_BATCH_VARIANTS = 500
LOG_ODDS_EPSILON = 1e-6
EXAMPLES = (
    "BRCA1 p.Arg1699Trp",
    "BRCA1 p.Ser1613Gly",
    "TP53 p.Arg175His",
    "CFTR p.Gly551Asp",
    "BRCA1 p.Pro871Leu",
)


@dataclass(frozen=True)
class Prediction:
    variant: str
    gene: str
    transcript_id: str
    hgvs_p: str
    probability: float
    decision: str
    lightgbm: float
    dual_tower: float
    contributions: list[tuple[str, float]]
    features: dict[str, float | str | None] = field(repr=False)
    gene_function: str | None = None


@dataclass(frozen=True)
class InvalidVariant:
    variant: str
    reason: str


@dataclass(frozen=True)
class OutOfScope:
    variant: str
    gene: str | None
    consequence: str
    reason: str


@dataclass(frozen=True)
class ParsedVcf:
    variants: list[tuple[int, Variant]]
    errors: list[tuple[int, str]]


def parse_vcf(text: str, max_variants: int = MAX_BATCH_VARIANTS) -> ParsedVcf:
    variants, errors = [], []
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("	")
        if len(fields) < 5:
            errors.append((line_number, "expected tab-separated CHROM POS ID REF ALT columns"))
            continue
        chrom, pos, _, ref, alts = (field.strip() for field in fields[:5])
        for alt in alts.split(","):
            if len(ref) != 1 or len(alt) != 1:
                errors.append((line_number, f"{ref}>{alt} is not an SNV (indels are planned)"))
                continue
            try:
                variants.append((line_number, parse_snv(f"{chrom}:{pos} {ref}>{alt}")))
            except ValueError as error:
                errors.append((line_number, str(error)))
            if len(variants) > max_variants:
                raise ValueError(f"the file has more than {max_variants} variants")
    return ParsedVcf(variants, errors)


def variant_text(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos} {ref}>{alt}"


def scope_reason(consequence: str) -> str:
    if consequence in PHASE_TWO:
        return f"{consequence} variants are not supported yet (planned for phase 2)"
    if consequence in NON_CODING:
        return f"{consequence} variants are permanently out of scope (non-coding)"
    return f"{consequence} variants are out of scope"


def classify(name: str, annotations: list[Annotation]) -> Annotation | OutOfScope:
    if not annotations:
        return OutOfScope(name, None, "intergenic", "no MANE Select transcript overlaps it")
    if any(a.consequence == "ref_mismatch" for a in annotations):
        raise ValueError(f"REF does not match the reference genome at {name}")
    missense = [a for a in annotations if a.consequence == "missense"]
    if len(missense) > 1:
        genes = ", ".join(a.gene_name for a in missense)
        return OutOfScope(
            name, genes, "multiple_genes", "it is missense in more than one MANE gene (excluded)"
        )
    if missense:
        return missense[0]
    annotation = min(annotations, key=lambda a: a.consequence not in PHASE_TWO)
    return OutOfScope(
        name, annotation.gene_name, annotation.consequence, scope_reason(annotation.consequence)
    )


def decide(probability: float, thresholds: dict) -> str:
    if probability >= thresholds["likely_pathogenic_min"]:
        return "likely_pathogenic"
    if probability <= thresholds["likely_benign_max"]:
        return "likely_benign"
    return "uncertain"


def explanation_rows(features: np.ndarray, background: np.ndarray, names: list[str]) -> np.ndarray:
    rows = [features[None, :]]
    for columns in FEATURE_GROUPS.values():
        index = [names.index(c) for c in columns]
        replaced = np.repeat(features[None, :], len(background), axis=0)
        replaced[:, index] = background[:, index]
        rows.append(replaced)
    return np.concatenate(rows)


def log_odds(probabilities: np.ndarray) -> np.ndarray:
    clipped = np.clip(probabilities.astype(np.float64), LOG_ODDS_EPSILON, 1 - LOG_ODDS_EPSILON)
    return np.log(clipped / (1 - clipped))


def group_contributions(probabilities: np.ndarray, background_size: int) -> list[tuple[str, float]]:
    logits = log_odds(probabilities)
    replaced = logits[1:].reshape(len(FEATURE_GROUPS), background_size).mean(axis=1)
    deltas = logits[0] - replaced
    pairs = [(group, round(float(d), 3)) for group, d in zip(FEATURE_GROUPS, deltas, strict=True)]
    return sorted(pairs, key=lambda item: -abs(item[1]))


def background_sample(train: pl.DataFrame, names: list[str], size: int = 32) -> np.ndarray:
    sample = train.sample(min(size, train.height), seed=BACKGROUND_SEED)
    return missense.to_matrix(sample, names).astype(np.float32)


class Predictor:
    def __init__(
        self,
        inputs: missense.Inputs,
        model_path: Path,
        train: pl.DataFrame,
        functions: dict[str, str] | None = None,
    ):
        transcripts = load_transcripts(inputs.transcripts)
        self.functions = functions or {}
        self.inputs = inputs
        self.index = TranscriptIndex(transcripts)
        self.by_gene = {t.gene_name.upper(): t for t in transcripts}
        self.session = ort.InferenceSession(str(model_path), providers=["CPUExecutionProvider"])
        metadata = {
            key: json.loads(value)
            for key, value in self.session.get_modelmeta().custom_metadata_map.items()
        }
        self.feature_names: list[str] = metadata["feature_names"]
        self.thresholds: dict = metadata["thresholds"]
        self.blocks = [variant_inputs.BLOCKS.index(b) for b in metadata["esm_blocks"]]
        grouped = sorted(c for columns in FEATURE_GROUPS.values() for c in columns)
        if grouped != sorted(self.feature_names):
            raise ValueError("explanation groups do not partition the model's features")
        self.background = background_sample(train, self.feature_names)

    def predict(self, text: str) -> list[Prediction | OutOfScope]:
        if SNV_PATTERN.match(text.strip()):
            return [self.predict_variant(*parse_snv(text))]
        gene, aa_ref, position, aa_alt = parse_protein_change(text)
        transcript = self.by_gene.get(gene)
        if transcript is None:
            raise ValueError(f"{gene} is not a MANE Select gene")
        variants = protein_change_snvs(transcript, aa_ref, position, aa_alt)
        if not variants:
            return [
                OutOfScope(
                    text.strip(),
                    transcript.gene_name,
                    "multi_nucleotide",
                    "this change needs more than one nucleotide change (not an SNV)",
                )
            ]
        return [self.predict_variant(*variant) for variant in variants]

    def predict_variant(self, chrom: str, pos: int, ref: str, alt: str) -> Prediction | OutOfScope:
        result = self.predict_many([(chrom, pos, ref, alt)])[0]
        if isinstance(result, InvalidVariant):
            raise ValueError(result.reason)
        return result

    def predict_many(
        self, variants: list[Variant]
    ) -> list[Prediction | OutOfScope | InvalidVariant]:
        outcomes: list[Annotation | OutOfScope | InvalidVariant] = []
        for chrom, pos, ref, alt in variants:
            name = variant_text(chrom, pos, ref, alt)
            try:
                outcomes.append(classify(name, self.index.annotate(chrom, pos, ref, alt)))
            except ValueError as error:
                outcomes.append(InvalidVariant(name, str(error)))
        pending = {
            variant: outcome
            for variant, outcome in zip(variants, outcomes, strict=True)
            if isinstance(outcome, Annotation)
        }
        predictions = self.score(pending) if pending else {}
        return [
            predictions[variant] if isinstance(outcome, Annotation) else outcome
            for variant, outcome in zip(variants, outcomes, strict=True)
        ]

    def score(self, pending: dict[Variant, Annotation]) -> dict[Variant, Prediction]:
        table = self.features(pending)
        matrix = missense.to_matrix(table, self.feature_names).astype(np.float32)
        features = table.select(self.feature_names).rows(named=True)
        embeddings: dict[str, np.ndarray | None] = {}
        predictions = {}
        for row, key in enumerate(table.select(missense.KEY).iter_rows()):
            annotation = pending[key]
            if annotation.transcript_id not in embeddings:
                embeddings[annotation.transcript_id] = self.embeddings(annotation.transcript_id)
            rows = explanation_rows(matrix[row], self.background, self.feature_names)
            protein = np.repeat(
                self.protein(embeddings[annotation.transcript_id], annotation.protein_position),
                len(rows),
                0,
            )
            outputs = self.session.run(OUTPUTS, {"features": rows, "protein": protein})
            probability = float(outputs[0][0])
            predictions[key] = Prediction(
                variant=variant_text(*key),
                gene=annotation.gene_name,
                transcript_id=annotation.transcript_id,
                hgvs_p=annotation.hgvs_p,
                probability=probability,
                decision=decide(probability, self.thresholds),
                lightgbm=float(outputs[1][0]),
                dual_tower=float(outputs[2][0]),
                contributions=group_contributions(outputs[0], len(self.background)),
                features=features[row],
                gene_function=self.functions.get(annotation.transcript_id),
            )
        return predictions

    def features(self, pending: dict[Variant, Annotation]) -> pl.DataFrame:
        rows = pl.DataFrame(
            [
                (*variant, a.gene_name, a.transcript_id, a.protein_position, a.aa_ref, a.aa_alt)
                for variant, a in pending.items()
            ],
            schema={
                "chrom": pl.String,
                "pos": pl.Int64,
                "ref": pl.String,
                "alt": pl.String,
                "gene_name": pl.String,
                "transcript_id": pl.String,
                "protein_position": pl.Int64,
                "aa_ref": pl.String,
                "aa_alt": pl.String,
            },
            orient="row",
        )
        joined = missense.join_sources(rows, self.inputs)
        return missense.finalize(esm_features(joined, self.inputs.esm_store))

    def embeddings(self, transcript_id: str) -> np.ndarray | None:
        if not store.is_done(self.inputs.esm_store, transcript_id):
            return None
        return np.asarray(store.load(self.inputs.esm_store, transcript_id)[1])

    def protein(self, embeddings: np.ndarray | None, position: int) -> np.ndarray:
        width = len(self.blocks) * 1280
        if embeddings is None:
            return np.zeros((1, width), dtype=np.float32)
        blocks = variant_inputs.variant_blocks(embeddings, position)[self.blocks]
        return blocks.reshape(1, width).astype(np.float32)


def parity_check(predictor: Predictor, table: pl.DataFrame, blocks_path: Path) -> dict:
    expected = predictor.session.run(
        OUTPUTS,
        {
            "features": missense.to_matrix(table, predictor.feature_names).astype(np.float32),
            "protein": variant_inputs.load(blocks_path, table)[:, predictor.blocks]
            .reshape(table.height, -1)
            .astype(np.float32),
        },
    )[0]
    feature_mismatches, probability_diffs, timings = 0, [], []
    for i, row in enumerate(table.iter_rows(named=True)):
        started = time.perf_counter()
        result = predictor.predict_variant(row["chrom"], row["pos"], row["ref"], row["alt"])
        if not isinstance(result, Prediction):
            feature_mismatches += len(predictor.feature_names)
            continue
        timings.append(round((time.perf_counter() - started) * 1000, 1))
        probability_diffs.append(abs(result.probability - float(expected[i])))
        for name in predictor.feature_names:
            served, stored = result.features[name], row[name]
            same = served == stored or (
                isinstance(served, float)
                and isinstance(stored, float)
                and (np.isclose(served, stored) or (np.isnan(served) and np.isnan(stored)))
            )
            feature_mismatches += not same
    timings.sort()
    return {
        "variants": table.height,
        "predicted": len(timings),
        "feature_mismatches": feature_mismatches,
        "max_probability_diff": round(max(probability_diffs, default=0.0), 7),
        "latency_ms": {
            "median": round(statistics.median(timings), 1) if timings else None,
            "p95": timings[int(0.95 * len(timings)) - 1] if timings else None,
            "max": timings[-1] if timings else None,
        },
    }
