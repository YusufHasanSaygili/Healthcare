import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
import torch
from sklearn.metrics import log_loss
from torch import nn

from healthcare.esm import variant_inputs
from healthcare.evaluate.metrics import classification_metrics
from healthcare.features.esm_features import ESM_FEATURES
from healthcare.features.missense import (
    AMINO_ACID_CODES,
    CATEGORICAL_FEATURES,
    FREQUENCY_FEATURES,
)
from healthcare.models.baseline import early_stopping_genes
from healthcare.split.genes import load_split

ESM_FORMS = {
    "none": [],
    "position": [0],
    "window3": [0, 1],
    "window7": [0, 2],
    "protein": [0, 3],
    "all": [0, 1, 2, 3],
}
LOG10_COLUMNS = {"f_gnomad_af", "f_gnomad_af_grpmax"}
LOG1P_COLUMNS = {"f_gnomad_ac", "f_gnomad_nhomalt", "f_gnomad_exome_an", "f_protein_length"}
AA_VOCAB = len(AMINO_ACID_CODES) + 1
MAX_EPOCHS = 60
ESM_SCALARS = tuple(ESM_FEATURES)
NO_FREQUENCY = tuple(FREQUENCY_FEATURES)
PATIENCE = 4
VARIANT_TYPE = "f_variant_type"
VARIANT_TYPE_COUNT = 10


def transform(column: str, values: np.ndarray) -> np.ndarray:
    if column in LOG10_COLUMNS:
        return np.log10(values + 1e-6)
    if column in LOG1P_COLUMNS:
        return np.log1p(values)
    return values


@dataclass
class TabularEncoder:
    columns: list[str]
    medians: np.ndarray
    means: np.ndarray
    stds: np.ndarray
    flag_columns: list[str]
    variant_types: int = 0

    @classmethod
    def fit(cls, table: pl.DataFrame, dropped: tuple[str, ...] = ()) -> "TabularEncoder":
        columns = [
            c
            for c in table.columns
            if c.startswith("f_")
            and c not in CATEGORICAL_FEATURES
            and c != VARIANT_TYPE
            and not any(c.startswith(d) for d in dropped)
        ]
        raw = cls.raw(table, columns)
        medians = np.nanmedian(raw, axis=0)
        filled = np.where(np.isnan(raw), medians, raw)
        flags = [c for c, has in zip(columns, np.isnan(raw).any(axis=0), strict=True) if has]
        types = VARIANT_TYPE_COUNT if VARIANT_TYPE in table.columns else 0
        return cls(columns, medians, filled.mean(axis=0), filled.std(axis=0) + 1e-6, flags, types)

    @staticmethod
    def raw(table: pl.DataFrame, columns: list[str]) -> np.ndarray:
        matrix = table.select([pl.col(c).cast(pl.Float64) for c in columns]).to_numpy()
        matrix = matrix.astype(np.float64, copy=True)
        for i, column in enumerate(columns):
            matrix[:, i] = transform(column, matrix[:, i])
        return matrix

    def encode(self, table: pl.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        raw = self.raw(table, self.columns)
        missing = np.isnan(raw)
        scaled = (np.where(missing, self.medians, raw) - self.means) / self.stds
        flags = missing[:, [self.columns.index(c) for c in self.flag_columns]]
        parts = [scaled, flags]
        if self.variant_types:
            codes = table[VARIANT_TYPE].cast(pl.Int64).to_numpy()
            parts.append(np.eye(self.variant_types)[codes])
        numeric = np.hstack(parts).astype(np.float32)
        amino_acids = np.stack(
            [
                table[c]
                .fill_null("")
                .replace_strict(AMINO_ACID_CODES, default=AA_VOCAB - 1)
                .to_numpy()
                for c in CATEGORICAL_FEATURES
            ],
            axis=1,
        ).astype(np.int64)
        return numeric, amino_acids

    @property
    def width(self) -> int:
        return len(self.columns) + len(self.flag_columns) + self.variant_types


class DualTower(nn.Module):
    def __init__(self, numeric_width: int, protein_blocks: int, config: dict):
        super().__init__()
        d, dropout = config["hidden"], config["dropout"]
        self.amino_acids = nn.Embedding(AA_VOCAB, config["aa_dim"])
        layers, width = [], numeric_width + 2 * config["aa_dim"]
        for _ in range(config["tabular_depth"]):
            layers += [nn.Linear(width, d), nn.GELU(), nn.Dropout(dropout)]
            width = d
        self.tabular = nn.Sequential(*layers)
        self.protein = (
            nn.Sequential(
                nn.Linear(protein_blocks * 1280, d),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(d, d),
                nn.GELU(),
            )
            if protein_blocks
            else None
        )
        fusion, width = [], 2 * d
        for _ in range(config["fusion_depth"]):
            fusion += [nn.Linear(width, d), nn.GELU(), nn.Dropout(dropout)]
            width = d
        self.fusion = nn.Sequential(*fusion, nn.Linear(width, 1))
        self.hidden = d

    def forward(self, numeric, amino_acids, protein, protein_mask):
        tabular = self.tabular(torch.cat([numeric, self.amino_acids(amino_acids).flatten(1)], 1))
        if self.protein is None:
            protein_out = torch.zeros_like(tabular)
        else:
            protein_out = self.protein(protein) * protein_mask[:, None]
        return self.fusion(torch.cat([tabular, protein_out], 1)).squeeze(1)


@dataclass
class Tensors:
    numeric: torch.Tensor
    amino_acids: torch.Tensor
    protein: torch.Tensor
    labels: torch.Tensor
    available: torch.Tensor


class Data:
    def __init__(
        self,
        splits_dir: Path,
        esm_blocks_path: Path,
        device: str = "cuda",
        loader: Callable[[Path, str], pl.DataFrame] = load_split,
        esm_columns: list[int] | None = None,
    ):
        train = loader(splits_dir, "train")
        stop = early_stopping_genes(train)
        self.tables = {
            "fit": train.filter(~pl.col("gene_name").is_in(stop)),
            "stop": train.filter(pl.col("gene_name").is_in(stop)),
            "validation": loader(splits_dir, "validation"),
        }
        self.esm_columns = esm_columns or list(range(len(variant_inputs.BLOCKS)))
        self.blocks = {
            name: variant_inputs.load(esm_blocks_path, table, esm_columns)
            for name, table in self.tables.items()
        }
        self.device = device
        self.cache: dict = {}

    def protein_blocks(self, name: str) -> torch.Tensor:
        key = ("protein", name)
        if key not in self.cache:
            self.cache[key] = torch.from_numpy(np.ascontiguousarray(self.blocks[name])).to(
                self.device
            )
        return self.cache[key]

    def encoded(self, dropped: tuple[str, ...]) -> tuple[TabularEncoder, dict]:
        key = ("tabular", dropped)
        if key not in self.cache:
            for stale in [k for k in self.cache if k[0] == "tabular"]:
                del self.cache[stale]
            encoder = TabularEncoder.fit(self.tables["fit"], dropped)
            arrays = {}
            for name, table in self.tables.items():
                numeric, amino_acids = encoder.encode(table)
                available = table["protein_position"].is_not_null().to_numpy()
                arrays[name] = (
                    torch.from_numpy(numeric).to(self.device),
                    torch.from_numpy(amino_acids).to(self.device),
                    torch.from_numpy(table["label"].to_numpy().astype(np.float32)).to(self.device),
                    torch.from_numpy(available.astype(np.float32)).to(self.device),
                )
            self.cache[key] = (encoder, arrays)
        return self.cache[key]

    def tensors(self, dropped: tuple[str, ...], esm_form: str) -> tuple[TabularEncoder, dict]:
        encoder, arrays = self.encoded(dropped)
        blocks = [self.esm_columns.index(b) for b in ESM_FORMS[esm_form]]
        tensors = {}
        for name, (numeric, amino_acids, labels, available) in arrays.items():
            protein = (
                self.protein_blocks(name)[:, blocks].flatten(1)
                if blocks
                else torch.zeros((labels.shape[0], 0), dtype=torch.float16, device=self.device)
            )
            tensors[name] = Tensors(numeric, amino_acids, protein, labels, available)
        return encoder, tensors


@torch.inference_mode()
def predict(model: DualTower, tensors: Tensors, batch: int = 8192) -> np.ndarray:
    model.eval()
    outputs = []
    for start in range(0, tensors.labels.shape[0], batch):
        end = min(start + batch, tensors.labels.shape[0])
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = model(
                tensors.numeric[start:end],
                tensors.amino_acids[start:end],
                tensors.protein[start:end].float(),
                tensors.available[start:end],
            )
        outputs.append(torch.sigmoid(logits.float()).cpu().numpy())
    return np.concatenate(outputs)


def train(
    data: Data, config: dict, seed: int, dropped: tuple[str, ...] = (), history: bool = False
) -> dict:
    started = time.monotonic()
    torch.manual_seed(seed)
    np.random.seed(seed)
    encoder, tensors = data.tensors(dropped, config["esm_form"])
    fit, stop = tensors["fit"], tensors["stop"]
    model = DualTower(encoder.width, len(ESM_FORMS[config["esm_form"]]), config).to(data.device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config["lr"], weight_decay=config["weight_decay"]
    )
    loss_fn = nn.BCEWithLogitsLoss()
    generator = torch.Generator(device=data.device).manual_seed(seed)
    best_loss, best_state, best_epoch, stale = float("inf"), None, 0, 0
    epochs = []
    n = fit.labels.shape[0]
    for epoch in range(MAX_EPOCHS):
        model.train()
        order = torch.randperm(n, device=data.device, generator=generator)
        for start in range(0, n, config["batch_size"]):
            index = order[start : start + config["batch_size"]]
            mask = (
                torch.rand(index.shape[0], device=data.device, generator=generator)
                >= config["tower_dropout"]
            ).float() * fit.available[index]
            with torch.autocast("cuda", dtype=torch.bfloat16):
                logits = model(
                    fit.numeric[index],
                    fit.amino_acids[index],
                    fit.protein[index].float(),
                    mask,
                )
                loss = loss_fn(logits.float(), fit.labels[index])
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        probabilities = predict(model, stop)
        if not np.isfinite(probabilities).all():
            raise FloatingPointError("non-finite predictions; training diverged")
        stop_loss = float(log_loss(stop.labels.cpu().numpy(), probabilities, labels=[0, 1]))
        if history:
            validation = tensors["validation"]
            scores = classification_metrics(
                validation.labels.cpu().numpy(), predict(model, validation)
            )
            epochs.append(
                {
                    "epoch": epoch + 1,
                    "stop_loss": round(stop_loss, 5),
                    "roc_auc": round(scores["roc_auc"], 4),
                    "pr_auc": round(scores["pr_auc"], 4),
                }
            )
        if stop_loss < best_loss - 1e-5:
            best_loss, best_epoch, stale = stop_loss, epoch, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= PATIENCE:
                break
    if best_state is None:
        raise FloatingPointError("early-stopping loss never improved")
    model.load_state_dict(best_state)
    validation = classification_metrics(
        tensors["validation"].labels.cpu().numpy(), predict(model, tensors["validation"])
    )
    return {
        "model": model,
        "encoder": encoder,
        "validation": validation,
        "best_epoch": best_epoch,
        "history": epochs,
        "stop_loss": best_loss,
        "parameters": sum(p.numel() for p in model.parameters()),
        "seconds": round(time.monotonic() - started, 1),
    }
