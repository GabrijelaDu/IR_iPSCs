#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterable, List, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "resubmission" / "configs" / "hl_revision_runs.json"
ACTIVE_CONFIG_PATH = CONFIG_PATH
CHECKPOINT_RE = re.compile(r"val_auroc=(?P<score>[0-9.]+)\.ckpt$")
DEFAULT_ATTRIBUTION_MODE = "final_logit_linearized"
VALID_ATTRIBUTION_MODES = (DEFAULT_ATTRIBUTION_MODE, "branch_signed")


def normalize_attribution_mode(mode: str | None) -> str:
    normalized = DEFAULT_ATTRIBUTION_MODE if mode is None else str(mode).strip()
    if normalized not in VALID_ATTRIBUTION_MODES:
        raise SystemExit(
            f"Unsupported attribution mode '{normalized}'. Available: {', '.join(VALID_ATTRIBUTION_MODES)}"
        )
    return normalized



def parse_attribution_modes(raw: str | Sequence[str] | None) -> List[str]:
    if raw is None:
        return [DEFAULT_ATTRIBUTION_MODE]
    if isinstance(raw, str):
        chunks = raw.split(",")
    else:
        chunks = []
        for item in raw:
            chunks.extend(str(item).split(","))
    modes: List[str] = []
    for chunk in chunks:
        chunk = chunk.strip()
        if not chunk:
            continue
        mode = normalize_attribution_mode(chunk)
        if mode not in modes:
            modes.append(mode)
    return modes or [DEFAULT_ATTRIBUTION_MODE]



def attribution_mode_prefix(base_prefix: str, cam_mode: str | None = None) -> str:
    mode = normalize_attribution_mode(cam_mode)
    if mode == DEFAULT_ATTRIBUTION_MODE:
        return base_prefix
    return f"{base_prefix}_{mode}"



def attribution_mode_suffix(cam_mode: str | None = None) -> str:
    mode = normalize_attribution_mode(cam_mode)
    if mode == DEFAULT_ATTRIBUTION_MODE:
        return ""
    return f"_{mode}"



def attribution_mode_label(cam_mode: str | None = None) -> str:
    mode = normalize_attribution_mode(cam_mode)
    return {
        DEFAULT_ATTRIBUTION_MODE: "final-logit oriented",
        "branch_signed": "branch-signed",
    }[mode]



def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> Dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def resolve_config_path(config_path: str | Path | None = None) -> Path:
    if config_path is None:
        return ACTIVE_CONFIG_PATH
    path = Path(config_path)
    return path if path.is_absolute() else REPO_ROOT / path



def set_active_config_path(config_path: str | Path) -> Path:
    global ACTIVE_CONFIG_PATH
    ACTIVE_CONFIG_PATH = resolve_config_path(config_path)
    return ACTIVE_CONFIG_PATH



def get_active_config_path() -> Path:
    return ACTIVE_CONFIG_PATH



def load_config(config_path: str | Path | None = None) -> Dict[str, Any]:
    resolved = resolve_config_path(config_path)
    if config_path is not None:
        set_active_config_path(resolved)
    with resolved.open() as handle:
        return json.load(handle)


def resolve_repo_path(path_value: str) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def resolve_run(run_name: str, config_path: str | Path | None = None) -> Dict[str, Any]:
    config = load_config(config_path)
    runs = config["runs"]
    if run_name not in runs:
        raise SystemExit(f"Unknown run '{run_name}'. Available: {', '.join(sorted(runs))}")
    spec = deep_merge(config["defaults"], runs[run_name])
    spec["name"] = run_name
    spec["metadata_csv"] = resolve_repo_path(spec["metadata_csv"])
    spec["fasta_path"] = resolve_repo_path(spec["fasta_path"])
    spec["parnet_weights"] = resolve_repo_path(spec["parnet_weights"])
    spec["result_dir"] = resolve_repo_path(spec["result_dir"])
    spec["training_env"] = resolve_repo_path(spec["training_env"])
    spec["modisco_env"] = resolve_repo_path(spec["modisco_env"])
    spec["interpretation_dir"] = spec["result_dir"] / spec["interpretation_dirname"]
    spec["modisco_executable"] = spec["modisco_env"] / "bin" / "modisco"
    spec["modisco_python"] = spec["modisco_env"] / "bin" / "python"
    spec["training_python"] = spec["training_env"] / "bin" / "python"
    spec["umap_outputs"] = [resolve_repo_path(path) for path in spec.get("umap_outputs", [])]
    return spec


def parse_fold_list(raw_value: str | None) -> List[int] | None:
    if raw_value is None:
        return None
    folds: List[int] = []
    for chunk in raw_value.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        fold = int(chunk)
        if fold < 1:
            raise SystemExit(f"Fold numbers are 1-based; got '{chunk}'.")
        folds.append(fold)
    return sorted(set(folds))


def checkpoint_score(path: Path) -> float:
    match = CHECKPOINT_RE.search(path.name)
    if not match:
        return float("-inf")
    return float(match.group("score"))


def checkpoint_candidates(result_dir: Path) -> List[Path]:
    return sorted(result_dir.glob("beststate_fold*/best_model_epoch=*_val_auroc=*.ckpt"))


def choose_checkpoint(spec: Mapping[str, Any], fold: int | None = None, explicit: str | None = None) -> Path:
    if explicit:
        path = resolve_repo_path(explicit)
        if not path.exists():
            raise SystemExit(f"Checkpoint not found: {path}")
        return path
    candidates = checkpoint_candidates(spec["result_dir"])
    if fold is not None:
        candidates = [path for path in candidates if path.parent.name == f"beststate_fold{fold}"]
    if not candidates:
        raise SystemExit(f"No checkpoints found under {spec['result_dir']}")
    return max(candidates, key=checkpoint_score)


def fold_output_dir(spec: Mapping[str, Any], fold: int) -> Path:
    return spec["result_dir"] / f"beststate_fold{fold}"


def shell_env_lines(mapping: Mapping[str, Path | str | int | float]) -> str:
    lines = []
    for key, value in mapping.items():
        lines.append(f"{key}={shlex.quote(str(value))}")
    return "\n".join(lines)


def print_run_listing() -> None:
    config = load_config()
    for run_name in sorted(config["runs"]):
        spec = resolve_run(run_name)
        print(f"{run_name}\t{spec['result_dir']}\t{spec['description']}")


def list_runs() -> List[Dict[str, Any]]:
    config = load_config()
    runs: List[Dict[str, Any]] = []
    for run_name in sorted(config["runs"]):
        spec = resolve_run(run_name)
        runs.append({
            "run": run_name,
            "description": spec["description"],
            "result_dir": str(spec["result_dir"]),
        })
    return runs


def command_resolve_run(args: argparse.Namespace) -> None:
    spec = resolve_run(args.run)
    if args.format == "json":
        payload = dict(spec)
        for key, value in list(payload.items()):
            if isinstance(value, Path):
                payload[key] = str(value)
            elif isinstance(value, list):
                payload[key] = [str(item) if isinstance(item, Path) else item for item in value]
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    env_mapping = {
        "RUN_NAME": spec["name"],
        "RESULT_DIR": spec["result_dir"],
        "INTERPRETATION_DIR": spec["interpretation_dir"],
        "MODISCO_PREFIX": spec["modisco_prefix"],
        "MODISCO_BIN": spec["modisco_executable"],
        "MODISCO_PYTHON": spec["modisco_python"],
        "TRAINING_PYTHON": spec["training_python"],
    }
    print(shell_env_lines(env_mapping))


def print_check(ok: bool, label: str, path: Path | None = None) -> None:
    status = "OK" if ok else "MISSING"
    suffix = f" :: {path}" if path is not None else ""
    print(f"[{status}] {label}{suffix}")


def validate_run(
    run: str,
    motif_db: str | None = None,
    require_trained: bool = False,
    require_umap: bool = False,
    require_modisco_inputs: bool = False,
    cam_modes: str | Sequence[str] | None = None,
    verbose: bool = True,
) -> List[str]:
    spec = resolve_run(run)
    missing: List[str] = []

    checks: List[tuple[str, Path]] = [
        ("metadata", spec["metadata_csv"]),
        ("fasta", spec["fasta_path"]),
        ("parnet weights", spec["parnet_weights"]),
        ("training env", spec["training_env"]),
        ("training python", spec["training_python"]),
        ("modisco env", spec["modisco_env"]),
        ("modisco python", spec["modisco_python"]),
        ("modisco executable", spec["modisco_executable"]),
    ]
    if motif_db:
        checks.append(("motif database", resolve_repo_path(motif_db)))

    for label, path in checks:
        ok = path.exists()
        if verbose:
            print_check(ok, label, path)
        if not ok:
            missing.append(label)

    result_dir = spec["result_dir"]
    if verbose:
        print_check(result_dir.exists(), "result directory", result_dir)
    if result_dir.exists():
        for fold in range(1, spec["n_splits"] + 1):
            fold_dir = fold_output_dir(spec, fold)
            if verbose:
                print_check(fold_dir.exists(), f"fold {fold} directory", fold_dir)
        candidates = checkpoint_candidates(result_dir)
        if candidates:
            best_ckpt = max(candidates, key=checkpoint_score)
            if verbose:
                print_check(True, "best checkpoint", best_ckpt)
        else:
            if verbose:
                print_check(False, "best checkpoint", result_dir)
            if require_trained:
                missing.append("best checkpoint")
        for umap_path in spec["umap_outputs"]:
            ok = umap_path.exists()
            if verbose:
                print_check(ok, "UMAP export", umap_path)
            if require_umap and not ok:
                missing.append(f"UMAP:{umap_path}")
        if require_modisco_inputs:
            for cam_mode in parse_attribution_modes(cam_modes):
                prefix = attribution_mode_prefix(spec["modisco_prefix"], cam_mode)
                for suffix in ("seq_left", "cam_left", "seq_right", "cam_right"):
                    path = spec["interpretation_dir"] / f"{prefix}_{suffix}.npz"
                    ok = path.exists()
                    if verbose:
                        print_check(ok, f"modisco input [{cam_mode}]", path)
                    if not ok:
                        missing.append(f"modisco:{path}")

    return missing


def command_validate(args: argparse.Namespace) -> None:
    missing = validate_run(
        run=args.run,
        motif_db=args.motif_db,
        require_trained=args.require_trained,
        require_umap=args.require_umap,
        require_modisco_inputs=args.require_modisco_inputs,
        cam_modes=args.cam_modes,
        verbose=True,
    )
    if missing:
        raise SystemExit(1)


def import_torch_with_retry() -> Any:
    import importlib

    try:
        torch = importlib.import_module("torch")
        if not hasattr(torch, "version") or not hasattr(torch, "__version__"):
            raise AttributeError("partially initialized torch module")
        return torch
    except AttributeError as exc:
        partial_torch = sys.modules.get("torch")
        if partial_torch is None:
            raise
        stale_modules = [
            name for name in list(sys.modules)
            if name == "torch" or name.startswith("torch.")
        ]
        for name in stale_modules:
            sys.modules.pop(name, None)
        torch = importlib.import_module("torch")
        if not hasattr(torch, "version") or not hasattr(torch, "__version__"):
            raise RuntimeError(
                "PyTorch import did not complete cleanly. Restart the notebook kernel and try again."
            ) from exc
        return torch



def filter_metadata_by_pir_cutoff(tab_met, spec: Mapping[str, Any], pd_mod=None):
    if pd_mod is None:
        import pandas as pd_mod  # type: ignore
    out = tab_met.copy()
    if "PIR_Nuc_baseline" not in out.columns:
        raise SystemExit("metadata_csv is missing PIR_Nuc_baseline, which is required for config-driven metadata filtering.")
    out["PIR_Nuc_baseline"] = pd_mod.to_numeric(out["PIR_Nuc_baseline"], errors="coerce")
    cutoff = spec.get("pir_nuc_baseline_min")
    if cutoff is None:
        return out
    return out.loc[out["PIR_Nuc_baseline"].notna() & (out["PIR_Nuc_baseline"] >= float(cutoff))].copy()


def ensure_training_imports() -> Dict[str, Any]:
    try:
        import json as _json  # noqa: F401
        import numpy as np
        import pandas as pd
        torch = import_torch_with_retry()
        import umap
        from sklearn.model_selection import StratifiedKFold
        from torch.utils.data import DataLoader
    except ModuleNotFoundError as exc:
        env = resolve_run(next(iter(load_config()["runs"])))["training_env"]
        raise SystemExit(
            "Training dependencies are missing in the current Python.\n"
            f"Use {env / 'bin' / 'python'} to run this script."
        ) from exc

    sys.path.insert(0, str(REPO_ROOT))
    import ir_toolkit.ir_toolkit as irt

    return {
        "np": np,
        "pd": pd,
        "torch": torch,
        "umap": umap,
        "StratifiedKFold": StratifiedKFold,
        "DataLoader": DataLoader,
        "irt": irt,
    }


def load_training_data(spec: Mapping[str, Any], deps: Mapping[str, Any]) -> Dict[str, Any]:
    import pandas as pd
    irt = deps["irt"]
    tab_met = pd.read_csv(spec["metadata_csv"], index_col=0)
    tab_met = filter_metadata_by_pir_cutoff(tab_met, spec, pd_mod=pd)
    source_data = irt.prepare_dataset(
        tab_met,
        spec["class_column"],
        tuple(spec["class_limits"]),
    )
    source_data = source_data.sample(frac=1, random_state=spec["random_seed"])
    source_seqs = irt.load_fasta_to_dict(
        str(spec["fasta_path"]),
        ids_keep=set(source_data.index),
        id_func=lambda value: value.split("_")[0],
    )
    all_ids = source_data.index.tolist()
    y = source_data["class"]
    return {
        "source_data": source_data,
        "source_seqs": source_seqs,
        "all_ids": all_ids,
        "y": y,
    }


def model_kwargs_from_spec(spec: Mapping[str, Any]) -> Dict[str, Any]:
    kwargs = dict(spec["model_kwargs"])
    kwargs["attn_hidden"] = kwargs.pop("attn_hidden")
    return kwargs


def training_kwargs_from_spec(spec: Mapping[str, Any], fold: int) -> Dict[str, Any]:
    kwargs = dict(spec["training_kwargs"])
    kwargs["feat_dim"] = spec["model_kwargs"]["feat_dim"]
    kwargs["attn_hidden_dim"] = spec["model_kwargs"]["attn_hidden"]
    kwargs["use_attention_pooling"] = spec["model_kwargs"]["use_attention_pooling"]
    kwargs["use_joint_classifier"] = spec["model_kwargs"]["use_joint_classifier"]
    kwargs["use_middle"] = spec["model_kwargs"]["use_middle"]
    kwargs["use_simple_fusion"] = spec["model_kwargs"]["use_simple_fusion"]
    kwargs["experiment_name"] = f"{spec['name']}_fold{fold}"
    return kwargs


def load_backbone(spec: Mapping[str, Any], deps: Mapping[str, Any]):
    torch = deps["torch"]
    irt = deps["irt"]
    device = "cuda" if torch.cuda.is_available() else "cpu"
    parnet_core = torch.load(str(spec["parnet_weights"]), map_location=device, weights_only=False)
    return irt.ParnetBackbone(parnet_core)


def load_best_model(
    spec: Mapping[str, Any],
    deps: Mapping[str, Any],
    checkpoint_path: Path,
):
    torch = deps["torch"]
    irt = deps["irt"]
    device = "cuda" if torch.cuda.is_available() else "cpu"
    backbone = load_backbone(spec, deps)
    model = irt.IntronEndLightning(
        model=irt.IntronEndCAMModel(backbone=backbone, **model_kwargs_from_spec(spec)),
        pos_weight=1.0,
    )
    state = torch.load(str(checkpoint_path), map_location=device)
    model.load_state_dict(state["state_dict"], strict=True)
    model.to(device)
    model.eval()
    return model


def build_loader(
    ids: Sequence[str],
    source_seqs,
    source_data,
    spec: Mapping[str, Any],
    deps: Mapping[str, Any],
    shuffle: bool,
    drop_last: bool,
    num_workers: int | None = None,
):
    irt = deps["irt"]
    DataLoader = deps["DataLoader"]
    dataset = irt.IntronsEndsDataset(
        list(ids),
        source_seqs,
        source_data,
        L_left=spec["flank_length"],
        L_right=spec["flank_length"],
        use_middle=False,
    )
    loader = DataLoader(
        dataset,
        batch_size=spec["batch_size"],
        shuffle=shuffle,
        num_workers=spec["num_workers"] if num_workers is None else num_workers,
        collate_fn=irt.introns_collate_fn,
        drop_last=drop_last,
    )
    return dataset, loader


def save_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def command_train(args: argparse.Namespace) -> None:
    spec = resolve_run(args.run)
    deps = ensure_training_imports()
    pd = deps["pd"]
    shutil.which("nvidia-smi")
    training = load_training_data(spec, deps)
    all_ids = training["all_ids"]
    y = training["y"]
    skf = deps["StratifiedKFold"](n_splits=spec["n_splits"], shuffle=True, random_state=spec["random_seed"])
    selected_folds = parse_fold_list(args.folds)

    dataset_export = spec["result_dir"] / "training_dataset.csv"
    spec["result_dir"].mkdir(parents=True, exist_ok=True)
    training["source_data"].to_csv(dataset_export)

    for fold_index, (train_idx, val_idx) in enumerate(skf.split(all_ids, y), start=1):
        if selected_folds and fold_index not in selected_folds:
            continue
        model_dir = fold_output_dir(spec, fold_index)
        existing_ckpts = list(model_dir.glob("best_model_epoch=*_val_auroc=*.ckpt")) if model_dir.exists() else []
        if existing_ckpts and not args.overwrite:
            print(f"Skipping fold {fold_index}: checkpoint already exists in {model_dir}")
            continue

        train_ids = [all_ids[i] for i in train_idx]
        val_ids = [all_ids[i] for i in val_idx]
        _, train_loader = build_loader(train_ids, training["source_seqs"], training["source_data"], spec, deps, shuffle=True, drop_last=True)
        _, val_loader = build_loader(val_ids, training["source_seqs"], training["source_data"], spec, deps, shuffle=False, drop_last=True)

        model_dir.mkdir(parents=True, exist_ok=True)
        split_df = pd.DataFrame(
            {
                "id": train_ids + val_ids,
                "split": ["train"] * len(train_ids) + ["val"] * len(val_ids),
                "fold": fold_index,
            }
        )
        split_df.to_csv(model_dir / "fold_split.csv", index=False)

        backbone = load_backbone(spec, deps)
        backbone.train()
        train_kwargs = training_kwargs_from_spec(spec, fold_index)
        save_json(model_dir / "hyperparams.json", train_kwargs)
        trainer, model, log_collector, checkpoint_callback = deps["irt"].setup_training(
            train_loader,
            backbone_model=backbone,
            **train_kwargs,
        )

        deps["irt"].run_training(trainer, model, train_loader, val_loader)
        best_model_path = Path(checkpoint_callback.best_model_path)
        if not best_model_path.exists():
            raise SystemExit(f"Training did not produce a best checkpoint for fold {fold_index}")

        reloaded_backbone = load_backbone(spec, deps)
        best_model = deps["irt"].IntronEndLightning.load_from_checkpoint(
            str(best_model_path),
            model=deps["irt"].IntronEndCAMModel(backbone=reloaded_backbone, **model_kwargs_from_spec(spec)),
        )
        shutil.copy2(best_model_path, model_dir / best_model_path.name)

        results = deps["irt"].evaluate_model(best_model, val_loader)
        pd.DataFrame(
            {
                "id": results["ids"],
                "prediction": results["predictions"],
                "label": results["labels"],
                "binary_prediction": results["binary_predictions"],
            }
        ).to_csv(model_dir / "evaluation_results.csv", index=False)
        deps["irt"].plot_training_history(log_collector, save_path=model_dir / "training_history.pdf", show=False)
        deps["irt"].visualize_cam_examples(results, training["source_seqs"], num_examples=5, save_path=model_dir / "cam_examples.pdf", show=False)
        deps["irt"].plot_roc_prc(results["labels"], results["predictions"], save_path=model_dir / "roc_prc_curves.pdf", show=False)


def train_run(run: str, folds: str | None = None, overwrite: bool = False) -> None:
    command_train(SimpleNamespace(run=run, folds=folds, overwrite=overwrite))


def extract_side_means(loader, model_backbone, torch_mod, pd_mod, device: str, concat_both: bool = True):
    ids_all: List[str] = []
    chunks = []
    with torch_mod.no_grad():
        for batch in loader:
            ids, left, left_mask, right, right_mask, _middle, _labels = batch
            left = left.to(device)
            right = right.to(device)
            left_mask = left_mask.to(device)
            right_mask = right_mask.to(device)

            left_feats = model_backbone(left)
            right_feats = model_backbone(right)

            left_den = left_mask.sum(dim=1, keepdim=True).clamp_min(1e-8)
            right_den = right_mask.sum(dim=1, keepdim=True).clamp_min(1e-8)
            left_vec = (left_feats * left_mask.unsqueeze(1).to(left_feats.dtype)).sum(dim=2) / left_den
            right_vec = (right_feats * right_mask.unsqueeze(1).to(right_feats.dtype)).sum(dim=2) / right_den

            embedding = torch_mod.cat([left_vec, right_vec], dim=1) if concat_both else left_vec
            chunks.append(embedding.cpu())
            ids_all.extend(ids)

    embeds = torch_mod.cat(chunks, dim=0).detach().cpu().numpy()
    return pd_mod.DataFrame(embeds, index=ids_all)


def command_export_umap(args: argparse.Namespace) -> None:
    spec = resolve_run(args.run)
    deps = ensure_training_imports()
    checkpoint_path = choose_checkpoint(spec, fold=args.fold, explicit=args.checkpoint)
    training = load_training_data(spec, deps)
    _dataset, loader = build_loader(
        training["all_ids"],
        training["source_seqs"],
        training["source_data"],
        spec,
        deps,
        shuffle=False,
        drop_last=False,
        num_workers=0,
    )
    model = load_best_model(spec, deps, checkpoint_path)
    device = "cuda" if deps["torch"].cuda.is_available() else "cpu"
    embeddings = extract_side_means(loader, model.model.backbone.eval(), deps["torch"], deps["pd"], device=device, concat_both=True)

    reducer = deps["umap"].UMAP(**spec["umap_kwargs"])
    transformed = reducer.fit_transform(embeddings.select_dtypes(include=["number"]))
    umap_df = deps["pd"].DataFrame(transformed, index=embeddings.index, columns=spec["umap_columns"])

    for output_path in spec["umap_outputs"]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        umap_df.to_csv(output_path, index=True)
        print(f"Saved UMAP export to {output_path}")


def export_umap_for_run(run: str, fold: int | None = None, checkpoint: str | None = None) -> None:
    command_export_umap(SimpleNamespace(run=run, fold=fold, checkpoint=checkpoint))


def prediction_table(loader, model, torch_mod, device: str, threshold: float):
    rows = []
    with torch_mod.no_grad():
        for batch in loader:
            ids, left, left_mask, right, right_mask, middle, _labels = batch
            left = left.to(device)
            left_mask = left_mask.to(device)
            right = right.to(device)
            right_mask = right_mask.to(device)
            middle = middle.to(device) if middle is not None else None
            outputs = model(left, left_mask, right, right_mask, middle)
            logits = outputs["seq_logit"]
            probs = torch_mod.sigmoid(logits.float()).view(-1)
            preds = (probs >= threshold).to(torch_mod.int32)
            rows.extend(
                {"id": seq_id, "prob": float(prob), "pred": int(pred)}
                for seq_id, prob, pred in zip(ids, probs.cpu(), preds.cpu())
            )
    return rows


def batch_item(batch: Any, key: int):
    return batch[key]


def fusion_branch_alphas(
    wrapped_model,
    torch_mod,
    logit_left,
    logit_right,
    logit_middle=None,
):
    components = [logit_left.detach(), logit_right.detach()]
    if getattr(wrapped_model, "use_middle", False) and logit_middle is not None:
        components.append(logit_middle.detach())

    joint = torch_mod.stack(components, dim=1).requires_grad_(True)
    with torch_mod.enable_grad():
        if getattr(wrapped_model, "use_joint_classifier", False):
            seq_logit = wrapped_model.joint_classifier(joint).squeeze(-1)
        elif getattr(wrapped_model, "use_simple_fusion", False):
            seq_logit = wrapped_model.simple_fusion(joint).squeeze(-1)
        else:
            seq_logit = joint.mean(dim=1)
        grads = torch_mod.autograd.grad(seq_logit.sum(), joint, create_graph=False, retain_graph=False)[0]

    return seq_logit.detach(), grads.detach()


def export_cams_npz_split(
    model,
    dataloader,
    out_dir: Path,
    prefix: str,
    torch_mod,
    np_mod,
    device: str,
    sign_multiplier: float,
    save_attention: bool,
    cam_mode: str = DEFAULT_ATTRIBUTION_MODE,
) -> Dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    model.eval()
    wrapped = getattr(model, "model", model)
    cam_mode = normalize_attribution_mode(cam_mode)

    seq_left_chunks = []
    seq_right_chunks = []
    cam_left_chunks = []
    cam_right_chunks = []
    attn_left_chunks = []
    attn_right_chunks = []

    for batch in dataloader:
        with torch_mod.no_grad():
            seq_left = batch_item(batch, 1).to(device, non_blocking=True)
            left_mask = batch_item(batch, 2).to(device, non_blocking=True).to(seq_left.dtype)
            seq_right = batch_item(batch, 3).to(device, non_blocking=True)
            right_mask = batch_item(batch, 4).to(device, non_blocking=True).to(seq_right.dtype)

            left_feats = wrapped.backbone(seq_left)
            right_feats = wrapped.backbone(seq_right)

            left_head = wrapped.left_cam_head
            left_attn_logits = left_head.attention_pooling(left_feats).squeeze(1)
            left_attn = torch_mod.softmax(left_attn_logits.masked_fill(left_mask == 0, -1e9), dim=-1)
            left_attn = (left_attn * left_mask) / left_attn.sum(dim=-1, keepdim=True).clamp_min(1e-8)
            left_w = left_head.gap_classifier.weight.squeeze(0)
            left_b = left_head.gap_classifier.bias.squeeze(0)
            left_cam = torch_mod.einsum("c,bcl->bl", left_w, left_feats)
            left_cam = (left_cam + left_b) * left_mask
            left_cam = left_attn * left_cam
            logit_left = left_cam.sum(dim=-1)

            right_head = wrapped.right_cam_head
            right_attn_logits = right_head.attention_pooling(right_feats).squeeze(1)
            right_attn = torch_mod.softmax(right_attn_logits.masked_fill(right_mask == 0, -1e9), dim=-1)
            right_attn = (right_attn * right_mask) / right_attn.sum(dim=-1, keepdim=True).clamp_min(1e-8)
            right_w = right_head.gap_classifier.weight.squeeze(0)
            right_b = right_head.gap_classifier.bias.squeeze(0)
            right_cam = torch_mod.einsum("c,bcl->bl", right_w, right_feats)
            right_cam = (right_cam + right_b) * right_mask
            right_cam = right_attn * right_cam
            logit_right = right_cam.sum(dim=-1)

        branch_left_cam = left_cam * sign_multiplier
        branch_right_cam = right_cam * sign_multiplier
        if cam_mode == DEFAULT_ATTRIBUTION_MODE:
            _seq_logit, fusion_grads = fusion_branch_alphas(
                wrapped,
                torch_mod,
                logit_left,
                logit_right,
                logit_middle=None,
            )
            left_alpha = fusion_grads[:, 0].to(left_cam.dtype)
            right_alpha = fusion_grads[:, 1].to(right_cam.dtype)
            export_left_cam = branch_left_cam * left_alpha.unsqueeze(1)
            export_right_cam = branch_right_cam * right_alpha.unsqueeze(1)
        else:
            export_left_cam = branch_left_cam
            export_right_cam = branch_right_cam
        left_cam_4c = export_left_cam.unsqueeze(1) * seq_left
        right_cam_4c = export_right_cam.unsqueeze(1) * seq_right

        seq_left_chunks.append(seq_left.detach().cpu().numpy())
        seq_right_chunks.append(seq_right.detach().cpu().numpy())
        cam_left_chunks.append(left_cam_4c.detach().cpu().numpy())
        cam_right_chunks.append(right_cam_4c.detach().cpu().numpy())

        if save_attention:
            attn_left_chunks.append(left_attn.unsqueeze(1).detach().cpu().numpy())
            attn_right_chunks.append(right_attn.unsqueeze(1).detach().cpu().numpy())

    outputs = {
        "seq_left": out_dir / f"{prefix}_seq_left.npz",
        "cam_left": out_dir / f"{prefix}_cam_left.npz",
        "seq_right": out_dir / f"{prefix}_seq_right.npz",
        "cam_right": out_dir / f"{prefix}_cam_right.npz",
    }

    np_mod.savez_compressed(outputs["seq_left"], arr_0=np_mod.concatenate(seq_left_chunks, axis=0).astype(np_mod.float32))
    np_mod.savez_compressed(outputs["cam_left"], arr_0=np_mod.concatenate(cam_left_chunks, axis=0).astype(np_mod.float32))
    np_mod.savez_compressed(outputs["seq_right"], arr_0=np_mod.concatenate(seq_right_chunks, axis=0).astype(np_mod.float32))
    np_mod.savez_compressed(outputs["cam_right"], arr_0=np_mod.concatenate(cam_right_chunks, axis=0).astype(np_mod.float32))

    if save_attention:
        attn_left = out_dir / f"{prefix}_attn_left.npz"
        attn_right = out_dir / f"{prefix}_attn_right.npz"
        np_mod.savez_compressed(attn_left, arr_0=np_mod.concatenate(attn_left_chunks, axis=0).astype(np_mod.float32))
        np_mod.savez_compressed(attn_right, arr_0=np_mod.concatenate(attn_right_chunks, axis=0).astype(np_mod.float32))
        outputs["attn_left"] = attn_left
        outputs["attn_right"] = attn_right

    return {key: str(path) for key, path in outputs.items()}


def command_export_modisco_inputs(args: argparse.Namespace) -> None:
    spec = resolve_run(args.run)
    deps = ensure_training_imports()
    checkpoint_path = choose_checkpoint(spec, fold=args.fold, explicit=args.checkpoint)
    training = load_training_data(spec, deps)
    _dataset, full_loader = build_loader(
        training["all_ids"],
        training["source_seqs"],
        training["source_data"],
        spec,
        deps,
        shuffle=False,
        drop_last=False,
        num_workers=0,
    )
    best_model = load_best_model(spec, deps, checkpoint_path)
    device = "cuda" if deps["torch"].cuda.is_available() else "cpu"
    rows = prediction_table(full_loader, best_model, deps["torch"], device=device, threshold=args.threshold)
    pred_table = deps["pd"].DataFrame(rows, columns=["id", "prob", "pred"])
    pred_table["pred"] = pred_table["pred"].astype(str)
    source_with_id = training["source_data"].copy()
    source_with_id["id"] = source_with_id.index

    if args.correct_only:
        merged = pred_table.merge(source_with_id, on="id")
        selected_ids = merged.loc[merged["pred"] == merged["class"], "id"].tolist()
    else:
        selected_ids = pred_table["id"].tolist()

    _correct_ds, correct_loader = build_loader(
        selected_ids,
        training["source_seqs"],
        training["source_data"],
        spec,
        deps,
        shuffle=False,
        drop_last=False,
        num_workers=0,
    )
    outputs_by_mode: Dict[str, Dict[str, str]] = {}
    for cam_mode in parse_attribution_modes(args.cam_modes):
        outputs_by_mode[cam_mode] = export_cams_npz_split(
            model=best_model,
            dataloader=correct_loader,
            out_dir=spec["interpretation_dir"],
            prefix=attribution_mode_prefix(spec["modisco_prefix"], cam_mode),
            torch_mod=deps["torch"],
            np_mod=deps["np"],
            device=device,
            sign_multiplier=spec["modisco_sign_multiplier"],
            save_attention=True,
            cam_mode=cam_mode,
        )

    pred_table.to_csv(spec["interpretation_dir"] / "prediction_table.csv", index=False)
    deps["pd"].DataFrame({"id": selected_ids}).to_csv(spec["interpretation_dir"] / "selected_ids.csv", index=False)
    save_json(
        spec["interpretation_dir"] / "export_manifest.json",
        {
            "run": spec["name"],
            "checkpoint": str(checkpoint_path),
            "correct_only": args.correct_only,
            "selected_ids": len(selected_ids),
            "outputs_by_mode": outputs_by_mode,
            "cam_modes": list(outputs_by_mode),
        },
    )
    for cam_mode, outputs in outputs_by_mode.items():
        print(f"CAM mode: {cam_mode}")
        for path in outputs.values():
            print(f"Saved {path}")


def export_modisco_inputs_for_run(
    run: str,
    fold: int | None = None,
    checkpoint: str | None = None,
    threshold: float = 0.5,
    correct_only: bool = True,
    cam_modes: str | Sequence[str] | None = None,
) -> None:
    command_export_modisco_inputs(
        SimpleNamespace(
            run=run,
            fold=fold,
            checkpoint=checkpoint,
            threshold=threshold,
            correct_only=correct_only,
            cam_modes=cam_modes,
        )
    )


def run_modisco_for_run(
    run: str,
    motif_db: str | Path,
    sides: str = "left,right",
    window: int = 256,
    size: int = 31,
    max_seqlets: int = 1000,
    skip_report: bool = False,
    check: bool = True,
    cam_mode: str = DEFAULT_ATTRIBUTION_MODE,
):
    spec = resolve_run(run)
    script_path = REPO_ROOT / "resubmission" / "scripts" / "run_modisco_report.sh"
    cmd = [
        str(script_path),
        "--run",
        run,
        "--motif-db",
        str(resolve_repo_path(motif_db)),
        "--side",
        sides,
        "--window",
        str(window),
        "--size",
        str(size),
        "--max-seqlets",
        str(max_seqlets),
        "--cam-mode",
        normalize_attribution_mode(cam_mode),
    ]
    if skip_report:
        cmd.append("--skip-report")
    env = os.environ.copy()
    env["PYTHON_BIN"] = str(spec["training_python"])
    return subprocess.run(cmd, cwd=str(REPO_ROOT), env=env, check=check)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean training and interpretation workflow for the half-life revision models.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("list-runs", help="List the configured revision runs.")

    resolve_parser = subparsers.add_parser("resolve-run", help="Resolve one run to absolute paths.")
    resolve_parser.add_argument("--run", required=True)
    resolve_parser.add_argument("--format", choices=("json", "shell"), default="json")

    validate_parser = subparsers.add_parser("validate", help="Validate paths and optional outputs for one run.")
    validate_parser.add_argument("--run", required=True)
    validate_parser.add_argument("--motif-db")
    validate_parser.add_argument("--require-trained", action="store_true")
    validate_parser.add_argument("--require-umap", action="store_true")
    validate_parser.add_argument("--require-modisco-inputs", action="store_true")
    validate_parser.add_argument("--cam-modes", default=DEFAULT_ATTRIBUTION_MODE)

    train_parser = subparsers.add_parser("train", help="Train one configured revision run.")
    train_parser.add_argument("--run", required=True)
    train_parser.add_argument("--folds", help="Comma-separated list of 1-based folds to train.")
    train_parser.add_argument("--overwrite", action="store_true")

    umap_parser = subparsers.add_parser("export-umap", help="Export UMAP embeddings for the best checkpoint of one run.")
    umap_parser.add_argument("--run", required=True)
    umap_parser.add_argument("--fold", type=int)
    umap_parser.add_argument("--checkpoint")

    modisco_parser = subparsers.add_parser("export-modisco-inputs", help="Export sequence/CAM arrays for tf-modisco.")
    modisco_parser.add_argument("--run", required=True)
    modisco_parser.add_argument("--fold", type=int)
    modisco_parser.add_argument("--checkpoint")
    modisco_parser.add_argument("--threshold", type=float, default=0.5)
    modisco_parser.add_argument("--correct-only", action="store_true", default=True)
    modisco_parser.add_argument("--all-examples", action="store_false", dest="correct_only")
    modisco_parser.add_argument("--cam-modes", default=DEFAULT_ATTRIBUTION_MODE)

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "list-runs":
        print_run_listing()
        return
    if args.command == "resolve-run":
        command_resolve_run(args)
        return
    if args.command == "validate":
        command_validate(args)
        return
    if args.command == "train":
        command_train(args)
        return
    if args.command == "export-umap":
        command_export_umap(args)
        return
    if args.command == "export-modisco-inputs":
        command_export_modisco_inputs(args)
        return
    raise SystemExit(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
