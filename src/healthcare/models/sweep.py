import json
import math
import random
import statistics
import time
from pathlib import Path

import numpy as np
import optuna
import torch
from optuna.storages import JournalStorage
from optuna.storages.journal import JournalFileBackend

from healthcare.esm.model import REVISION as ESM_REVISION
from healthcare.evaluate.metrics import rounded
from healthcare.models import dual_tower

SEED = 20261012
STUDY = "dual_tower_missense_v1"
RANDOM_FRACTION = 0.4
CONFIRM_TOP = 10
CONFIRM_SEEDS = (0, 1, 2, 3, 4)
SPACE = {
    "hidden": [128, 192, 256, 384, 512],
    "tabular_depth": [1, 2, 3],
    "fusion_depth": [1, 2],
    "aa_dim": [4, 8, 16],
    "batch_size": [256, 512, 1024, 2048],
    "esm_form": list(dual_tower.ESM_FORMS)[1:],
}
RANGES = {
    "dropout": (0.0, 0.4),
    "tower_dropout": (0.0, 0.5),
    "lr": (1e-4, 3e-3),
    "weight_decay": (1e-6, 1e-2),
}
LOG_SCALE = {"lr", "weight_decay"}


def suggest(trial: optuna.Trial) -> dict:
    config = {name: trial.suggest_categorical(name, choices) for name, choices in SPACE.items()}
    for name, (low, high) in RANGES.items():
        config[name] = trial.suggest_float(name, low, high, log=name in LOG_SCALE)
    return config


def random_config(rng: random.Random) -> dict:
    config = {name: rng.choice(choices) for name, choices in SPACE.items()}
    for name, (low, high) in RANGES.items():
        if name in LOG_SCALE:
            config[name] = 10 ** rng.uniform(math.log10(low), math.log10(high))
        else:
            config[name] = rng.uniform(low, high)
    return config


def is_enqueued(study: optuna.Study, trial: optuna.Trial) -> bool:
    return "fixed_params" in study.get_trials(deepcopy=False)[trial.number].system_attrs


def open_study(journal: Path) -> optuna.Study:
    journal.parent.mkdir(parents=True, exist_ok=True)
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    return optuna.create_study(
        study_name=STUDY,
        storage=JournalStorage(JournalFileBackend(str(journal))),
        sampler=optuna.samplers.TPESampler(seed=SEED, multivariate=True, n_startup_trials=20),
        direction="maximize",
        load_if_exists=True,
    )


def run_sweep(
    data: dual_tower.Data, journal: Path, results: Path, hours: float, max_trials: int | None
) -> dict:
    started = time.monotonic()
    study = open_study(journal)
    rng = random.Random(SEED + len(study.trials))
    ran, stopped = 0, None
    with results.open("a", encoding="utf-8") as handle:
        while time.monotonic() - started < hours * 3600 and (
            max_trials is None or ran < max_trials
        ):
            waiting = study.get_trials(deepcopy=False, states=(optuna.trial.TrialState.WAITING,))
            if not waiting and rng.random() < RANDOM_FRACTION:
                study.enqueue_trial(random_config(rng))
            trial = study.ask()
            is_random = is_enqueued(study, trial)
            trial.set_user_attr("random", is_random)
            config = suggest(trial)
            try:
                outcome = dual_tower.train(data, config, seed=SEED)
            except (torch.OutOfMemoryError, FloatingPointError):
                torch.cuda.empty_cache()
                study.tell(trial, state=optuna.trial.TrialState.FAIL)
                continue
            except torch.AcceleratorError as error:
                study.tell(trial, state=optuna.trial.TrialState.FAIL)
                stopped = f"CUDA error at trial {trial.number}: {str(error).splitlines()[0]}"
                break
            metrics = rounded(outcome["validation"])
            trial.set_user_attr("validation", metrics)
            trial.set_user_attr("best_epoch", outcome["best_epoch"])
            trial.set_user_attr("seconds", outcome["seconds"])
            study.tell(trial, metrics["roc_auc"])
            handle.write(
                json.dumps(
                    {
                        "trial": trial.number,
                        "random": is_random,
                        "config": config,
                        "validation": metrics,
                        "best_epoch": outcome["best_epoch"],
                        "seconds": outcome["seconds"],
                    }
                )
                + "\n"
            )
            handle.flush()
            ran += 1
            del outcome
    return sweep_summary(study, ran, time.monotonic() - started) | {"stopped_by_error": stopped}


def completed(study: optuna.Study) -> list:
    return [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]


def sweep_summary(study: optuna.Study, ran: int, seconds: float) -> dict:
    trials = sorted(completed(study), key=lambda t: -t.value)
    return {
        "study": STUDY,
        "trials_completed_total": len(trials),
        "trials_this_run": ran,
        "random_trials": sum(bool(t.user_attrs.get("random")) for t in trials),
        "hours_this_run": round(seconds / 3600, 2),
        "best": [
            {
                "trial": t.number,
                "roc_auc": t.value,
                "pr_auc": t.user_attrs["validation"]["pr_auc"],
                "params": t.params,
            }
            for t in trials[:5]
        ],
    }


def seeds_summary(runs: list[dict]) -> dict:
    values = {m: [r["validation"][m] for r in runs] for m in ("roc_auc", "pr_auc", "mcc", "ece")}
    return {
        metric: {"mean": round(statistics.mean(v), 4), "std": round(statistics.pstdev(v), 4)}
        for metric, v in values.items()
    }


def run_confirmation(
    data: dual_tower.Data, journal: Path, artifacts_dir: Path, reports_dir: Path
) -> dict:
    started = time.monotonic()
    study = open_study(journal)
    top = sorted(completed(study), key=lambda t: -t.value)[:CONFIRM_TOP]
    candidates = []
    for trial in top:
        runs = [dual_tower.train(data, trial.params, seed) for seed in CONFIRM_SEEDS]
        candidates.append(
            {
                "trial": trial.number,
                "params": trial.params,
                "sweep_roc_auc": trial.value,
                "seeds": seeds_summary(runs),
                "runs": runs,
            }
        )
    for candidate in candidates:
        roc = candidate["seeds"]["roc_auc"]
        candidate["score"] = round(roc["mean"] - roc["std"], 4)
    winner = max(candidates, key=lambda c: c["score"])

    ablations = {}
    for name, dropped, form in (
        ("no_frequency", dual_tower.NO_FREQUENCY, winner["params"]["esm_form"]),
        ("no_protein_tower", (), "none"),
        ("no_esm", dual_tower.ESM_SCALARS, "none"),
        ("no_frequency_no_esm", dual_tower.NO_FREQUENCY + dual_tower.ESM_SCALARS, "none"),
    ):
        config = winner["params"] | {"esm_form": form}
        runs = [dual_tower.train(data, config, seed, dropped) for seed in CONFIRM_SEEDS]
        ablations[name] = seeds_summary(runs)

    best_run = winner["runs"][0]
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dict": best_run["model"].state_dict(),
            "config": winner["params"],
            "encoder": {
                key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in vars(best_run["encoder"]).items()
            },
            "seed": CONFIRM_SEEDS[0],
            "esm_revision": ESM_REVISION,
        },
        artifacts_dir / "dual_tower.pt",
    )
    baseline = json.loads((reports_dir / "baseline_lightgbm.json").read_text())
    report = {
        "model": "Dual tower (missense, phase 1)",
        "test_set_opened": False,
        "study": STUDY,
        "trials_completed": len(completed(study)),
        "confirmation_seeds": list(CONFIRM_SEEDS),
        "selection": "highest validation ROC-AUC mean minus std over 5 seeds",
        "note": "validation was also used to rank the top configs, so these numbers are mildly "
        "optimistic; the locked test set (task 6.1) gives the unbiased estimate",
        "winner": {
            "trial": winner["trial"],
            "params": winner["params"],
            "seeds": winner["seeds"],
            "parameters": best_run["parameters"],
        },
        "top_candidates": [
            {
                "trial": c["trial"],
                "sweep_roc_auc": round(c["sweep_roc_auc"], 4),
                "score": c["score"],
                "roc_auc": c["seeds"]["roc_auc"],
                "pr_auc": c["seeds"]["pr_auc"],
            }
            for c in sorted(candidates, key=lambda c: -c["score"])
        ],
        "ablations": ablations,
        "lightgbm_validation": {
            "full": baseline["ablations"]["full"]["validation"],
            "no_frequency": baseline["ablations"]["no_frequency"]["validation"],
        },
        "esm_revision": ESM_REVISION,
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "dual_tower.json").write_text(json.dumps(report, indent=2))
    return report
