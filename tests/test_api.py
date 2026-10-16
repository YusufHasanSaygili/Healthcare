import gzip

import pytest
from fastapi.testclient import TestClient

from healthcare.features.missense import UNIPROT_GROUPS
from healthcare.serve import api
from healthcare.serve.pipeline import InvalidVariant, OutOfScope, Prediction

THRESHOLDS = {"likely_benign_max": 0.1, "likely_pathogenic_min": 0.9}
FEATURES = {name: 0.0 for name in UNIPROT_GROUPS} | {
    "f_uniprot_domain": 1.0,
    "f_gnomad_af": 0.0,
    "f_phylop": 7.5,
    "f_esm_llr": -12.0,
    "f_grantham": 101.0,
}


def prediction(variant: str) -> Prediction:
    return Prediction(
        variant=variant,
        gene="BRCA1",
        transcript_id="ENST00000357654.9",
        hgvs_p="p.Arg1699Trp",
        probability=0.95,
        decision="likely_pathogenic",
        lightgbm=0.9,
        dual_tower=0.92,
        contributions=[("esm_scores", 3.5), ("allele_frequency", 0.4)],
        features=FEATURES,
    )


class FakePredictor:
    def __init__(self):
        self.batches = []

    def predict(self, text):
        if text == "bad input":
            raise ValueError("expected an SNV like 17:43045712 G>A")
        if text.startswith("1:"):
            return [OutOfScope("1:1000 A>G", None, "intergenic", "no MANE Select transcript")]
        return [prediction("17:43063931 G>A")]

    def predict_many(self, variants):
        self.batches.append(variants)
        return [
            InvalidVariant(f"{c}:{p} {r}>{a}", "REF does not match")
            if p == 999
            else prediction(f"{c}:{p} {r}>{a}")
            for c, p, r, a in variants
        ]


def report() -> dict:
    def metrics(roc):
        return {
            "roc_auc": roc,
            "pr_auc": 0.97,
            "ece": 0.004,
            "gene_bootstrap_95ci": {"roc_auc": [roc - 0.002, roc + 0.002]},
        }

    def result(roc, ece=0.004):
        return {
            "variants": 100,
            "models": {
                "ensemble_calibrated": metrics(roc) | {"ece": ece},
                "lightgbm": metrics(0.989),
                "dual_tower": metrics(0.988),
                "no_esm/ensemble_uncalibrated": metrics(0.985),
                "no_frequency/ensemble_uncalibrated": metrics(0.969),
            },
            "decisions": {"uncertain": {"share": 0.12}},
        }

    return {
        "onnx": "missense_v1.onnx",
        "variant_type": "missense",
        "data_versions": {"clinvar": "2025-10-06"},
        "results": {
            "test": result(0.990),
            "temporal_all": result(0.989, ece=0.033),
            "temporal_test_gene": result(0.988),
        },
    }


@pytest.fixture
def predictor():
    return FakePredictor()


@pytest.fixture
def client(predictor):
    return TestClient(api.create_app(predictor, api.model_card(report(), THRESHOLDS)))


def test_health_reports_model_and_scope(client):
    body = client.get("/api/health").json()
    assert body == {"status": "ok", "model": "missense_v1.onnx", "variant_types": ["missense"]}


def test_predict_returns_decision_contributions_and_evidence(client):
    response = client.post("/api/predict", json={"variant": "BRCA1 p.Arg1699Trp"})
    assert response.status_code == 200
    body = response.json()
    result = body["results"][0]
    assert result["status"] == "predicted"
    assert result["decision"] == "likely_pathogenic"
    assert result["contributions"][0] == {"group": "esm_scores", "log_odds": 3.5}
    assert result["evidence"]["protein_regions"] == ["domain"]
    assert "not for clinical" in body["disclaimer"]


def test_predict_reports_out_of_scope_variants(client):
    result = client.post("/api/predict", json={"variant": "1:1000 A>G"}).json()["results"][0]
    assert result["status"] == "out_of_scope" and result["consequence"] == "intergenic"


@pytest.mark.parametrize(
    ("payload", "status"),
    [
        ({"variant": "bad input"}, 422),
        ({"variant": "x"}, 422),
        ({"variant": "A" * 101}, 422),
        ({}, 422),
    ],
)
def test_predict_rejects_invalid_input(client, payload, status):
    assert client.post("/api/predict", json=payload).status_code == status


def test_predict_rejects_oversized_json_bodies(client):
    response = client.post("/api/predict", json={"variant": "BRCA1 p.Arg1699Trp", "x": "A" * 20000})
    assert response.status_code == 413


def test_predict_error_message_reaches_the_client(client):
    response = client.post("/api/predict", json={"variant": "bad input"})
    assert "expected an SNV" in response.json()["detail"]


VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\n"
    "chr17\t43063931\t.\tG\tA,T\n"
    "17\t43063932\t.\tGA\tG\n"
    "17\t999\t.\tC\tT\n"
    "MT\t10\t.\tA\tG\n"
    "17\t5\n"
)


def test_batch_predicts_snvs_and_reports_line_errors(client, predictor):
    response = client.post("/api/predict/batch", content=VCF.encode())
    assert response.status_code == 200
    body = response.json()
    assert predictor.batches == [
        [("17", 43063931, "G", "A"), ("17", 43063931, "G", "T"), ("17", 999, "C", "T")]
    ]
    assert body["variants"] == 3
    assert [item["line"] for item in body["results"]] == [3, 3]
    assert [error["line"] for error in body["errors"]] == [4, 5, 6, 7]
    assert "not an SNV" in body["errors"][0]["message"]
    assert "REF does not match" in body["errors"][1]["message"]
    assert body["error_count"] == 4


def test_batch_rejects_oversized_uploads(client, monkeypatch):
    monkeypatch.setattr(api, "MAX_UPLOAD_BYTES", 100)
    response = client.post("/api/predict/batch", content=VCF.encode() * 10)
    assert response.status_code == 413


def test_batch_rejects_too_many_variants(client):
    rows = "".join(f"17\t{1000 + i}\t.\tA\tG\n" for i in range(api.MAX_BATCH_VARIANTS + 1))
    response = client.post("/api/predict/batch", content=rows.encode())
    assert response.status_code == 413
    assert "more than" in response.json()["detail"]


def test_batch_rejects_compressed_and_binary_uploads(client):
    assert client.post("/api/predict/batch", content=gzip.compress(VCF.encode())).status_code == 415
    assert client.post("/api/predict/batch", content=b"\xff\xfe\x00bad").status_code == 422


def test_batch_with_no_variants_returns_an_empty_result(client, predictor):
    body = client.post("/api/predict/batch", content=b"#CHROM\tPOS\n").json()
    assert body["variants"] == 0 and body["results"] == [] and predictor.batches == []


def test_cors_allows_only_the_local_vite_origin(client):
    allowed = client.get("/api/health", headers={"Origin": "http://localhost:5173"})
    other = client.get("/api/health", headers={"Origin": "http://evil.example"})
    assert allowed.headers["access-control-allow-origin"] == "http://localhost:5173"
    assert "access-control-allow-origin" not in other.headers


def test_model_card_lists_metrics_ablations_and_limitations(client):
    card = client.get("/api/model-card").json()
    assert card["thresholds"] == THRESHOLDS
    assert [row["roc_auc"] for row in card["evaluation"]] == [0.990, 0.989, 0.988]
    assert {a["name"] for a in card["ablations"]} >= {"Without allele frequency", "Without ESM-2"}
    assert any("0.033" in line for line in card["limitations"])
    assert any("0.969" in line for line in card["limitations"])


def test_openapi_documents_the_plain_text_batch_body(client):
    schema = client.get("/openapi.json").json()
    body = schema["paths"]["/api/predict/batch"]["post"]["requestBody"]
    assert "text/plain" in body["content"]
