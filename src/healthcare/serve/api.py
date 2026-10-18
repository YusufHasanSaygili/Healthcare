import time
from typing import Annotated, Literal, Protocol

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field
from starlette.concurrency import run_in_threadpool

from healthcare.features.missense import UNIPROT_GROUPS
from healthcare.serve.pipeline import (
    MAX_BATCH_VARIANTS,
    InvalidVariant,
    OutOfScope,
    Prediction,
    Variant,
    parse_vcf,
)

ALLOWED_ORIGINS = ["http://localhost:5173", "http://127.0.0.1:5173"]
MAX_UPLOAD_BYTES = 2_000_000
MAX_JSON_BYTES = 10_000
MAX_REPORTED_ERRORS = 100
DISCLAIMER = (
    "For research and educational use only; not for clinical diagnosis or treatment decisions."
)

Decision = Literal["likely_pathogenic", "uncertain", "likely_benign"]
EvidenceGroup = Literal[
    "allele_frequency",
    "conservation",
    "esm_scores",
    "gene_constraint",
    "protein_regions",
    "amino_acid_change",
    "protein_position",
    "truncation",
    "splicing",
]
VariantType = Literal[
    "missense",
    "synonymous",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "canonical_splice",
    "frameshift",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_indel",
]
ONNX_FILE = "all_types_v2.onnx"


class PredictionService(Protocol):
    def predict(self, text: str) -> list[Prediction | OutOfScope]: ...

    def predict_many(
        self, variants: list[Variant]
    ) -> list[Prediction | OutOfScope | InvalidVariant]: ...


class PredictRequest(BaseModel):
    variant: str = Field(
        min_length=3,
        max_length=100,
        description="Genomic SNV (17:43063931 G>A) or protein change (BRCA1 p.Arg1699Trp)",
        examples=["BRCA1 p.Arg1699Trp"],
    )


class Contribution(BaseModel):
    group: EvidenceGroup
    log_odds: float


class Evidence(BaseModel):
    gnomad_af: float | None
    phylop: float | None
    esm_llr: float | None
    grantham: float | None
    protein_regions: list[str]
    nmd_escape: float | None = None
    protein_retained: float | None = None


class PredictedVariant(BaseModel):
    status: Literal["predicted"] = "predicted"
    variant: str
    gene: str
    transcript_id: str
    hgvs_p: str | None
    consequence: VariantType
    probability: float
    decision: Decision
    lightgbm: float
    dual_tower: float
    contributions: list[Contribution]
    evidence: Evidence
    gene_function: str | None = Field(
        default=None, description="UniProtKB/Swiss-Prot function summary; display only"
    )


class OutOfScopeVariant(BaseModel):
    status: Literal["out_of_scope"] = "out_of_scope"
    variant: str
    gene: str | None
    consequence: str
    reason: str


VariantResult = Annotated[PredictedVariant | OutOfScopeVariant, Field(discriminator="status")]


class PredictResponse(BaseModel):
    query: str
    results: list[VariantResult]
    milliseconds: float
    disclaimer: str = DISCLAIMER


class BatchItem(BaseModel):
    line: int
    result: VariantResult


class BatchError(BaseModel):
    line: int
    message: str


class BatchResponse(BaseModel):
    variants: int
    results: list[BatchItem]
    errors: list[BatchError]
    error_count: int
    milliseconds: float
    disclaimer: str = DISCLAIMER


class Health(BaseModel):
    status: Literal["ok"] = "ok"
    model: str
    variant_types: list[str]


class Thresholds(BaseModel):
    likely_benign_max: float
    likely_pathogenic_min: float


class EvaluationRow(BaseModel):
    name: str
    variant_type: str
    variants: int
    roc_auc: float
    roc_auc_ci: tuple[float, float]
    pr_auc: float
    ece: float
    uncertain_share: float


class Ablation(BaseModel):
    name: str
    variant_type: str
    roc_auc: float


class ModelCard(BaseModel):
    model: str
    variant_types: list[str]
    thresholds: Thresholds
    data_versions: dict[str, str]
    evaluation: list[EvaluationRow]
    ablations: list[Ablation]
    limitations: list[str]
    disclaimer: str = DISCLAIMER


EVALUATION_SETS = {
    "test": "Gene-disjoint test",
    "temporal_all": "Temporal (added to ClinVar after October 2025)",
    "temporal_test_gene": "Temporal, test genes only",
}
COMPARISONS = {
    "lightgbm": "LightGBM only",
    "dual_tower": "Dual tower only",
    "no_esm": "LightGBM without ESM-2",
    "no_frequency": "LightGBM without allele frequency",
}


def model_card(report: dict, data_versions: dict, thresholds: dict) -> ModelCard:
    results = report["results"]
    rows = [
        EvaluationRow(
            name=name,
            variant_type=variant_type,
            variants=metrics["n"],
            roc_auc=metrics["roc_auc"],
            roc_auc_ci=tuple(metrics["gene_bootstrap_95ci"]["roc_auc"]),
            pr_auc=metrics["pr_auc"],
            ece=metrics["ece"],
            uncertain_share=metrics["decisions"]["uncertain"],
        )
        for key, name in EVALUATION_SETS.items()
        for variant_type, metrics in results[key]["final"].items()
        if metrics["eligible"]
    ]
    ablations = [
        Ablation(name=name, variant_type=variant_type, roc_auc=value)
        for key, name in COMPARISONS.items()
        for variant_type, value in results["test"]["compare_roc_auc"][key].items()
        if value is not None
    ]
    test = results["test"]["final"]
    loss = [t for t in ("frameshift", "stop_gained", "canonical_splice") if test[t]["eligible"]]
    no_frequency = results["test"]["compare_roc_auc"]["no_frequency"]
    return ModelCard(
        model=ONNX_FILE,
        variant_types=sorted(report["routes"]),
        thresholds=Thresholds(**thresholds),
        data_versions=data_versions,
        evaluation=rows,
        ablations=ablations,
        limitations=[
            "Only coding variants in MANE Select transcripts (GRCh38) are scored: missense, "
            "synonymous, stop gained/lost, start lost, canonical splice (±1/2), frameshift and "
            "in-frame indels. Deep intronic, UTR, regulatory and structural variants are out of "
            "scope.",
            "For loss-of-function types ClinVar benign labels largely follow population frequency; "
            "without allele frequency test ROC-AUC drops to "
            + ", ".join(f"{no_frequency[t]:.2f} ({t})" for t in loss)
            + ".",
            "Start loss, stop loss and in-frame delins have too few labelled examples of one class "
            "to evaluate; their predictions are the least certain.",
            "Labels come from ClinVar and inherit its biases toward well-studied genes.",
            "The model is not clinically validated.",
        ],
    )


def evidence(prediction: Prediction) -> Evidence:
    features = prediction.features
    return Evidence(
        gnomad_af=features["f_gnomad_af"],
        phylop=features["f_phylop"],
        esm_llr=features["f_esm_llr"],
        grantham=features["f_grantham"],
        protein_regions=[
            name.removeprefix("f_uniprot_") for name in UNIPROT_GROUPS if features[name] == 1.0
        ],
        nmd_escape=features.get("f_nmd_escape"),
        protein_retained=features.get("f_protein_retained"),
    )


def to_result(result: Prediction | OutOfScope) -> PredictedVariant | OutOfScopeVariant:
    if isinstance(result, OutOfScope):
        return OutOfScopeVariant(
            variant=result.variant,
            gene=result.gene,
            consequence=result.consequence,
            reason=result.reason,
        )
    return PredictedVariant(
        variant=result.variant,
        gene=result.gene,
        transcript_id=result.transcript_id,
        hgvs_p=result.hgvs_p,
        consequence=result.consequence,
        probability=result.probability,
        decision=result.decision,
        lightgbm=result.lightgbm,
        dual_tower=result.dual_tower,
        contributions=[
            Contribution(group=group, log_odds=value) for group, value in result.contributions
        ],
        evidence=evidence(result),
        gene_function=result.gene_function,
    )


async def read_limited(request: Request, limit: int) -> bytes:
    declared = request.headers.get("content-length")
    if declared is not None and declared.isdigit() and int(declared) > limit:
        raise HTTPException(413, f"upload is larger than {limit} bytes")
    body = bytearray()
    async for chunk in request.stream():
        body.extend(chunk)
        if len(body) > limit:
            raise HTTPException(413, f"upload is larger than {limit} bytes")
    return bytes(body)


def decode_vcf(body: bytes) -> str:
    if body[:2] == b"\x1f\x8b":
        raise HTTPException(415, "compressed VCF is not supported; upload plain text")
    try:
        return body.decode("utf-8")
    except UnicodeDecodeError as error:
        raise HTTPException(422, "the VCF is not UTF-8 text") from error


def batch_response(
    parsed_errors: list[tuple[int, str]],
    lines: list[int],
    outcomes: list[Prediction | OutOfScope | InvalidVariant],
    started: float,
) -> BatchResponse:
    errors = list(parsed_errors)
    items = []
    for line, outcome in zip(lines, outcomes, strict=True):
        if isinstance(outcome, InvalidVariant):
            errors.append((line, outcome.reason))
        else:
            items.append(BatchItem(line=line, result=to_result(outcome)))
    errors.sort()
    return BatchResponse(
        variants=len(lines),
        results=items,
        errors=[BatchError(line=line, message=message) for line, message in errors][
            :MAX_REPORTED_ERRORS
        ],
        error_count=len(errors),
        milliseconds=round((time.perf_counter() - started) * 1000, 1),
    )


def create_app(
    predictor: PredictionService | None = None, card: ModelCard | None = None
) -> FastAPI:
    app = FastAPI(
        title="Healthcare variant pathogenicity API",
        description=DISCLAIMER,
        version="1",
    )
    app.add_middleware(
        CORSMiddleware,
        allow_origins=ALLOWED_ORIGINS,
        allow_methods=["GET", "POST"],
        allow_headers=["Content-Type"],
    )

    @app.middleware("http")
    async def limit_json_body(request: Request, call_next):
        declared = request.headers.get("content-length", "")
        too_large = declared.isdigit() and int(declared) > MAX_JSON_BYTES
        if request.url.path == "/api/predict" and too_large:
            return JSONResponse({"detail": "request body is too large"}, status_code=413)
        return await call_next(request)

    @app.get("/api/health")
    def health() -> Health:
        return Health(model=card.model, variant_types=card.variant_types)

    @app.post("/api/predict")
    def predict(request: PredictRequest) -> PredictResponse:
        started = time.perf_counter()
        try:
            results = predictor.predict(request.variant)
        except ValueError as error:
            raise HTTPException(422, str(error)) from error
        return PredictResponse(
            query=request.variant,
            results=[to_result(result) for result in results],
            milliseconds=round((time.perf_counter() - started) * 1000, 1),
        )

    @app.post(
        "/api/predict/batch",
        openapi_extra={
            "requestBody": {
                "required": True,
                "description": f"Plain-text VCF, at most {MAX_UPLOAD_BYTES} bytes and "
                f"{MAX_BATCH_VARIANTS} variants; never stored",
                "content": {"text/plain": {"schema": {"type": "string"}}},
            }
        },
    )
    async def predict_batch(request: Request) -> BatchResponse:
        started = time.perf_counter()
        text = decode_vcf(await read_limited(request, MAX_UPLOAD_BYTES))
        try:
            parsed = parse_vcf(text)
        except ValueError as error:
            raise HTTPException(413, str(error)) from error
        lines = [line for line, _ in parsed.variants]
        variants = [variant for _, variant in parsed.variants]
        outcomes = await run_in_threadpool(predictor.predict_many, variants) if variants else []
        return batch_response(parsed.errors, lines, outcomes, started)

    @app.get("/api/model-card")
    def get_model_card() -> ModelCard:
        return card

    return app


def openapi_schema() -> dict:
    return create_app().openapi()
