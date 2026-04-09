#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
from tqdm.auto import tqdm

from lets_plot import *

REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE_PATH = REPO_ROOT / "resubmission" / "scripts" / "hl_revision_pipeline.py"

spec = importlib.util.spec_from_file_location("hl_revision_pipeline", PIPELINE_PATH)
hlrev = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(hlrev)

LetsPlot.setup_html()
base_size = 14
theme_settings = theme(
    axis_text=element_text(size=base_size),
    axis_title=element_text(size=base_size * 1.2),
    legend_title=element_text(size=base_size * 1.2),
    legend_text=element_text(size=base_size * 0.9),
    plot_title=element_text(size=base_size * 1.4),
    axis_ticks_y=element_line(color='black', size=0.5),
    axis_line_y=element_line(color='black', size=0.5),
)

NEGATIVE_COLOR = "#a76794"
POSITIVE_COLOR = "#2a1f3f"
CLASS_COLORS = None
CAM_CLIP = 0.1


def cam_values_for_mode(pack, cam_mode):
    cam_mode = hlrev.normalize_attribution_mode(cam_mode)
    if cam_mode == hlrev.DEFAULT_ATTRIBUTION_MODE:
        return pack["cam_final_contrib"]
    return pack["cam_branch_contrib"]



def assign_plot_classes(df, value_col, class_limits, negative_label, positive_label):
    low, high = class_limits
    out = df.copy()
    out["plot_class"] = "else"
    out.loc[out[value_col] < low, "plot_class"] = negative_label
    out.loc[out[value_col] >= high, "plot_class"] = positive_label
    return out


def plot_2d_density_contours_continuous(
    x, y, cmap_name="viridis", kde_smooth=1, levels=11, linewidth=2.3,
    base_alpha=1, min_alpha=0.5, alpha_decay="linear", decay_rate=2.0
):
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)
    kde.set_bandwidth(bw_method=kde.factor * kde_smooth)
    x_margin = (x.max() - x.min()) * 0.1
    y_margin = (y.max() - y.min()) * 0.1
    xmin, xmax = x.min() - x_margin, x.max() + x_margin
    ymin, ymax = y.min() - y_margin, y.max() + y_margin
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    min_level = np.percentile(zz, 50)
    levels_vals = np.linspace(min_level, zz.max(), levels)
    if alpha_decay == "linear":
        alphas = np.linspace(base_alpha, min_alpha, len(levels_vals))
    else:
        exp_range = np.linspace(0, 1, len(levels_vals))
        alphas = base_alpha * np.exp(-decay_rate * exp_range)
    cmap = plt.get_cmap(cmap_name)
    for i, lev in enumerate(levels_vals):
        plt.contour(xx, yy, zz, levels=[lev], colors=[cmap(i / max(levels - 1, 1))], linewidths=linewidth, alpha=alphas[::-1][i])


def plot_2dd_scatter(
    tab_data, contour_col, x_name, y_name, contour_col_classes, scat_alpha=0.3,
    dlevels=11, linewidth=2.3, kde_smooth=1.5, fig_size=(6.5, 5.25),
    show_plot=True, file_name=None
):
    plt.figure(figsize=fig_size)
    sns.scatterplot(data=tab_data, x=x_name, y=y_name, edgecolor=None, color="gray", alpha=scat_alpha, s=8, legend=False, rasterized=True)
    tab_data = tab_data[tab_data[contour_col].isin(list(contour_col_classes.keys()))].copy()
    for cls, cls_col in contour_col_classes.items():
        subset = tab_data[tab_data[contour_col] == cls]
        if subset.shape[0] < 3:
            continue
        cls_cmap = LinearSegmentedColormap.from_list(cls, ["white", cls_col], N=dlevels)
        plot_2d_density_contours_continuous(subset[x_name], subset[y_name], cmap_name=cls_cmap, kde_smooth=kde_smooth, levels=dlevels, linewidth=linewidth, alpha_decay="exponent")
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    if file_name:
        plt.savefig(file_name, bbox_inches="tight", dpi=300)
    if show_plot:
        plt.show()
    else:
        plt.close()


def _head_bundle(deps, head, feats, mask):
    torch_mod = deps["torch"]
    B, C, L = feats.shape
    if mask is None:
        mask = torch_mod.ones(B, L, device=feats.device, dtype=feats.dtype)
    mask = mask.to(feats.dtype)
    has_attn = hasattr(head, "attention_pooling")
    if has_attn:
        attn_logits = head.attention_pooling(feats).squeeze(1)
        attn = torch_mod.softmax(attn_logits.masked_fill(mask == 0, -1e9), dim=-1)
        attn = (attn * mask) / attn.sum(dim=-1, keepdim=True).clamp_min(1e-8)
        pooled = torch_mod.einsum("bcl,bl->bc", feats, attn)
        denom = None
    else:
        attn = None
        denom = mask.sum(dim=-1, keepdim=True).clamp_min(1e-8)
        pooled = (feats * mask.unsqueeze(1)).sum(-1) / denom
    logits = head.gap_classifier(pooled).squeeze(-1)
    w = head.gap_classifier.weight.squeeze(0)
    b = head.gap_classifier.bias.squeeze(0)
    cam_old = torch_mod.einsum("c,bcl->bl", w, feats) * mask
    cam_signed = (cam_old + b) * mask
    if has_attn:
        cam_contrib = attn * cam_signed
    else:
        cam_contrib = (cam_signed / denom) * mask
    recon = cam_contrib.sum(dim=-1)
    return dict(cam_old=cam_old, cam_signed=cam_signed, attn=attn, cam_contrib=cam_contrib, seq_logit=logits, recon=recon)


def _fusion_alphas(deps, model_core, logit_left, logit_right):
    torch_mod = deps["torch"]
    joint = torch_mod.stack([logit_left.detach(), logit_right.detach()], dim=1).requires_grad_(True)
    with torch_mod.enable_grad():
        if getattr(model_core, "use_joint_classifier", False):
            joint_logit = model_core.joint_classifier(joint).squeeze(-1)
        elif getattr(model_core, "use_simple_fusion", False):
            joint_logit = model_core.simple_fusion(joint).squeeze(-1)
        else:
            joint_logit = joint.mean(dim=1)
        grads = torch_mod.autograd.grad(joint_logit.sum(), joint, create_graph=False, retain_graph=False)[0]
    return joint_logit.detach(), grads.detach()


def cam_bundle_for_sample(deps, lightning_model, left_onehot, left_mask, right_onehot, right_mask, device="cuda"):
    torch_mod = deps["torch"]
    was_train = lightning_model.training
    lightning_model.eval()
    model_core = lightning_model.model
    with torch_mod.no_grad():
        left_onehot = left_onehot.to(device)
        right_onehot = right_onehot.to(device)
        left_mask = left_mask.to(device) if left_mask is not None else None
        right_mask = right_mask.to(device) if right_mask is not None else None
        feats_left = model_core.backbone(left_onehot)
        feats_right = model_core.backbone(right_onehot)
        pack_left = _head_bundle(deps, model_core.left_cam_head, feats_left, left_mask)
        pack_right = _head_bundle(deps, model_core.right_cam_head, feats_right, right_mask)
    joint_logit, grads = _fusion_alphas(deps, model_core, pack_left["seq_logit"], pack_right["seq_logit"])
    pack_left["fusion_alpha"] = grads[:, 0]
    pack_right["fusion_alpha"] = grads[:, 1]
    pack_left["cam_branch_contrib"] = pack_left["cam_contrib"]
    pack_right["cam_branch_contrib"] = pack_right["cam_contrib"]
    pack_left["cam_final_contrib"] = pack_left["cam_contrib"] * pack_left["fusion_alpha"].unsqueeze(1)
    pack_right["cam_final_contrib"] = pack_right["cam_contrib"] * pack_right["fusion_alpha"].unsqueeze(1)
    if was_train:
        lightning_model.train()
    return {"left": pack_left, "right": pack_right, "joint_logit": joint_logit}


def compute_roc_prc(y_true, y_pred):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    return {
        "roc": {"fpr": fpr, "tpr": tpr},
        "prc": {"precision": precision, "recall": recall},
    }


def export_plots_for_run(run, fold=None, checkpoint=None, top_n=50, cam_min_length=300, cam_mode=hlrev.DEFAULT_ATTRIBUTION_MODE):
    run_spec = hlrev.resolve_run(run)
    cam_mode = hlrev.normalize_attribution_mode(cam_mode)
    cam_suffix = hlrev.attribution_mode_suffix(cam_mode)
    modisco_prefix = hlrev.attribution_mode_prefix(run_spec["modisco_prefix"], cam_mode)
    plot_dir = run_spec["result_dir"] / "plot_exports"
    plot_dir.mkdir(parents=True, exist_ok=True)

    if run_spec["class_column"] == "stability_ohe":
        negative_label = "Unstable"
        positive_label = "Stable"
    else:
        negative_label = "SL-RIs"
        positive_label = "LL-RIs"

    class_colors = {negative_label: NEGATIVE_COLOR, positive_label: POSITIVE_COLOR}

    deps = hlrev.ensure_training_imports()
    checkpoint_path = hlrev.choose_checkpoint(run_spec, fold=fold, explicit=checkpoint)
    training = hlrev.load_training_data(run_spec, deps)
    source_data = training["source_data"]
    source_seqs = training["source_seqs"]

    umap_path = next(path for path in run_spec["umap_outputs"] if path.exists())
    tab_embeds = pd.read_csv(umap_path, index_col=0)
    tab_embeds.index.name = "EVENT"
    tab_meta_plot = pd.read_csv(run_spec["metadata_csv"], index_col=0)
    tab_meta_plot = hlrev.filter_metadata_by_pir_cutoff(tab_meta_plot, run_spec, pd_mod=pd)
    tab_umap = tab_meta_plot.merge(tab_embeds, left_index=True, right_index=True, how="inner")
    tab_umap = assign_plot_classes(tab_umap, run_spec["class_column"], run_spec["class_limits"], negative_label, positive_label)
    x_name, y_name = run_spec["umap_columns"]
    plot_2dd_scatter(tab_umap, "plot_class", x_name=x_name, y_name=y_name, contour_col_classes=class_colors, file_name=plot_dir / f"{run_spec['name']}_umap_density.png")

    best_model = hlrev.load_best_model(run_spec, deps, checkpoint_path)
    device = "cuda" if deps["torch"].cuda.is_available() else "cpu"
    eval_path = checkpoint_path.parent / "evaluation_results.csv"
    res_val = pd.read_csv(eval_path)
    res_val = res_val[res_val["label"] == res_val["binary_prediction"]].copy()
    length_df = pd.DataFrame({"id": list(source_seqs.keys()), "length": [len(seq) for seq in source_seqs.values()]})
    res_val = res_val.merge(length_df, on="id", how="left")
    res_val = res_val[res_val["length"] >= cam_min_length].copy()
    top_pos = res_val[res_val.label == 1].sort_values(by="prediction", ascending=False).head(top_n).copy()
    top_neg = res_val[res_val.label == 0].sort_values(by="prediction", ascending=True).head(top_n).copy()
    binder_grps = []
    class_to_df = {positive_label: top_pos, negative_label: top_neg}
    for grp_name, grp_df in class_to_df.items():
        grp_ids = grp_df.id.tolist()
        ds_grp = deps["irt"].IntronsEndsDataset(grp_ids, source_seqs, source_data, L_left=run_spec["flank_length"], L_right=run_spec["flank_length"], use_middle=False)
        binder_seqs = []
        for seq_id in tqdm(grp_ids, total=len(grp_ids)):
            packs = cam_bundle_for_sample(deps, best_model, ds_grp[seq_id]["left_onehot"].unsqueeze(0), ds_grp[seq_id]["left_mask"].unsqueeze(0), ds_grp[seq_id]["right_onehot"].unsqueeze(0), ds_grp[seq_id]["right_mask"].unsqueeze(0), device=device)
            binder_cams = []
            for side, cam_key in [("left", "cam_left"), ("right", "cam_right")]:
                cam_contrib = cam_values_for_mode(packs[side], cam_mode).squeeze(0).detach().cpu().numpy()
                binder_cams.append(pd.DataFrame({"id": seq_id, "class": grp_name, "cam_head": f"{cam_key}_contrib", "position": np.arange(1, cam_contrib.size + 1), "cam_scores": cam_contrib}))
            binder_seqs.append(pd.concat(binder_cams, axis=0))
        binder_grps.append(pd.concat(binder_seqs, axis=0))
    tab_topcam = pd.concat(binder_grps, axis=0)
    tab_topcam["cam_scores"] = tab_topcam["cam_scores"].clip(-CAM_CLIP, CAM_CLIP)
    cam_plot = (
        ggplot(data=tab_topcam)
        + geom_tile(aes(x="position", y="id", fill="cam_scores"), tooltips="none", width=1)
        + scale_fill_gradient2(low=NEGATIVE_COLOR, mid="white", high=POSITIVE_COLOR, midpoint=0)
        + facet_grid(x="cam_head", y="class", scales="free")
        + labs(x="Position at intron end", y="Introns", fill="CAM score")
        + theme_settings
        + theme(axis_text_y=element_blank())
        + ggsize(1120, 400)
    )
    ggsave(cam_plot, path=str(plot_dir), filename=f"{run_spec['name']}_CAMprofiles{cam_suffix}.svg")

    speck_candidates = [
        REPO_ROOT / "resubmission/data/external/ENCORE/HepG2-HeLa_FeaturesMatrix.csv",
        REPO_ROOT / "data/speckle_rbps/HepG2-HeLa_FeaturesMatrix.csv",
    ]
    speck_path = next((path for path in speck_candidates if path.exists()), None)
    speck_rbps = []
    if speck_path is not None:
        speck_tab = pd.read_csv(speck_path)
        if "Speckles" in speck_tab.columns and "RBP" in speck_tab.columns:
            speck_rbps = speck_tab.dropna(subset=["Speckles"])["RBP"].tolist()

    modisco_tables = []
    for side in ["left", "right"]:
        motif_html = run_spec["interpretation_dir"] / f"{modisco_prefix}_{side}_modisco" / "report" / "motifs.html"
        dfs = pd.read_html(motif_html, flavor="lxml")
        binder = []
        for i in range(0, 3):
            tmp = dfs[0][["pattern", "num_seqlets", f"match{i}", f"qval{i}"]].copy()
            tmp = tmp.rename(columns={f"match{i}": "match", f"qval{i}": "qval"})
            tmp["match_n"] = str(i + 1)
            binder.append(tmp)
        side_df = pd.concat(binder, axis=0)
        side_df["qval_log"] = -np.log10(side_df["qval"])
        side_df["match_name"] = side_df["match"].astype(str).str.split("_").str[0]
        side_df["cam_side"] = side
        modisco_tables.append(side_df)
    dfs_bpes = pd.concat(modisco_tables, axis=0)
    dfs_bpes["significant"] = dfs_bpes["qval"] < 0.05
    dfs_bpes["label"] = np.where(dfs_bpes["significant"], dfs_bpes["match_name"], np.nan)
    dfs_bpes["significance"] = np.where(dfs_bpes["significant"], "Significant", "NS")
    dfs_bpes["pattern_name"] = dfs_bpes["pattern"].astype(str).str.split("_").str[0]
    dfs_bpes["in_speckle"] = dfs_bpes["match_name"].isin(speck_rbps).map({True: "yes", False: "no"})
    hits_plot = (
        ggplot(data=dfs_bpes)
        + geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="red")
        + geom_point(aes(x="num_seqlets", y="qval_log", fill="significance", alpha="significance"), shape=21)
        + scale_fill_manual(values={"NS": "grey", "Significant": "red"})
        + scale_alpha_manual(values={"NS": 0.1, "Significant": 1})
        + geom_text_repel(aes(x="num_seqlets", y="qval_log", label="label", color="in_speckle"), size=6, na_text="")
        + scale_color_manual(values={"yes": "red", "no": "black"})
        + coord_cartesian(xlim=[0, 100], ylim=[-0.2, 4])
        + facet_grid(y="pattern_name")
        + theme_settings
        + ggsize(550, 700)
    )
    ggsave(hits_plot, path=str(plot_dir), filename=f"{run_spec['name']}_modisco_hits{cam_suffix}.svg")

    val_metrics = []
    val_curves_prc = []
    val_curves_roc = []
    model_label = run_spec["name"]
    for fld_dir in sorted(run_spec["result_dir"].glob("beststate_fold*")):
        eval_file = fld_dir / "evaluation_results.csv"
        if not eval_file.exists():
            continue
        res_val_fold = pd.read_csv(eval_file)
        val_auc = roc_auc_score(res_val_fold["label"], res_val_fold["prediction"])
        val_ap = average_precision_score(res_val_fold["label"], res_val_fold["prediction"])
        fold_name = fld_dir.name.replace("beststate_fold", "")
        val_metrics.append({"model": model_label, "fold": fold_name, "ROC": val_auc, "PRC": val_ap})
        curves = compute_roc_prc(res_val_fold["label"], res_val_fold["prediction"])
        val_curves_roc.append(pd.DataFrame({"model": model_label, "fold": fold_name, "tpr": curves["roc"]["tpr"], "fpr": curves["roc"]["fpr"]}))
        val_curves_prc.append(pd.DataFrame({"model": model_label, "fold": fold_name, "precision": curves["prc"]["precision"], "recall": curves["prc"]["recall"]}))

    tab_metrics = pd.DataFrame(val_metrics).melt(value_vars=["ROC", "PRC"], value_name="AU", var_name="metrics", id_vars=["model", "fold"])
    tab_curves_roc = pd.concat(val_curves_roc, axis=0)
    tab_curves_prc = pd.concat(val_curves_prc, axis=0)
    tab_curves_prc_plt = tab_curves_prc.query("precision > 0").copy()
    pos_frac = tab_curves_prc_plt["precision"].min()
    best_prc_row = max(val_metrics, key=lambda row: row["PRC"])
    best_roc_row = max(val_metrics, key=lambda row: row["ROC"])
    prc_annot = pd.DataFrame({
        "x": [0.98],
        "y": [0.98],
        "label": [f"Best fold {best_prc_row['fold']}: AUPRC={best_prc_row['PRC']:.3f}"],
    })
    roc_annot = pd.DataFrame({
        "x": [0.98],
        "y": [0.04],
        "label": [f"Best fold {best_roc_row['fold']}: AUROC={best_roc_row['ROC']:.3f}"],
    })
    perf_prc = (
        ggplot(tab_curves_prc_plt)
        + geom_line(aes(group="fold", y="precision", x="recall"), color=POSITIVE_COLOR, size=0.8, alpha=0.7, tooltips="none")
        + geom_hline(yintercept=pos_frac, linetype="dashed", color=POSITIVE_COLOR, alpha=0.7, size=1)
        + geom_text(aes(x="x", y="y", label="label"), data=prc_annot, color=POSITIVE_COLOR, hjust=1, vjust=1)
        + theme_settings
        + labs(x="Recall", y="Precision", title=f"{run_spec['name']} PRC across folds")
        + coord_cartesian(xlim=[0, 1], ylim=[0, 1])
        + theme_settings
        + ggsize(550, 450)
    )
    perf_roc = (
        ggplot(tab_curves_roc)
        + geom_line(aes(group="fold", y="tpr", x="fpr"), color=POSITIVE_COLOR, size=0.8, alpha=0.7)
        + geom_abline(intercept=0, slope=1, linetype="dashed", size=1)
        + geom_text(aes(x="x", y="y", label="label"), data=roc_annot, color=POSITIVE_COLOR, hjust=1, vjust=0)
        + theme_settings
        + labs(y="True Positive Rate", x="False Positive Rate", title=f"{run_spec['name']} ROC across folds")
        + coord_cartesian(xlim=[0, 1], ylim=[0, 1])
        + theme_settings
        + ggsize(550, 450)
    )
    ggsave(perf_prc, filename=f"{run_spec['name']}_auprc_folds.svg", path=str(plot_dir))
    ggsave(perf_roc, filename=f"{run_spec['name']}_auroc_folds.svg", path=str(plot_dir))
    tab_metrics.to_csv(plot_dir / f"{run_spec['name']}_performance_metrics.csv", index=False)
    tab_curves_prc.to_csv(plot_dir / f"{run_spec['name']}_prc_curves.csv", index=False)
    tab_curves_roc.to_csv(plot_dir / f"{run_spec['name']}_roc_curves.csv", index=False)

    print(f"Saved plot exports to {plot_dir} (CAM mode: {cam_mode})")


def build_parser():
    parser = argparse.ArgumentParser(description="Export downstream plots for an HL revision run.")
    parser.add_argument("--run", required=True)
    parser.add_argument("--fold", type=int)
    parser.add_argument("--checkpoint")
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--cam-min-length", type=int, default=300)
    parser.add_argument("--cam-mode", default=hlrev.DEFAULT_ATTRIBUTION_MODE)
    return parser


def main():
    args = build_parser().parse_args()
    export_plots_for_run(
        args.run,
        fold=args.fold,
        checkpoint=args.checkpoint,
        top_n=args.top_n,
        cam_min_length=args.cam_min_length,
        cam_mode=args.cam_mode,
    )


if __name__ == "__main__":
    main()
