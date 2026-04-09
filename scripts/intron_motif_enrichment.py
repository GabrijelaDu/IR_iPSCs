from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import re
import traceback
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import Markdown, display
from numpy.lib.stride_tricks import sliding_window_view


DNA_ALPHABET = "ACGT"
BASE_TO_INDEX = {base: i for i, base in enumerate(DNA_ALPHABET)}
INDEX_TO_BASE = {i: base for base, i in BASE_TO_INDEX.items()}
COMPLEMENT_INDEX = np.array([3, 2, 1, 0], dtype=np.int8)


@dataclass(frozen=True)
class MemeMotif:
    name: str
    matrix: np.ndarray
    background: np.ndarray

    @property
    def width(self) -> int:
        return int(self.matrix.shape[0])


@dataclass
class AnalysisBundle:
    selected_introns: pd.DataFrame
    thresholds: pd.DataFrame
    results: pd.DataFrame
    pairwise: pd.DataFrame
    hits: pd.DataFrame
    null_distributions: dict[tuple[str, str], pd.DataFrame]


@dataclass
class GroupComparisonBundle:
    summary: pd.DataFrame
    resamples: pd.DataFrame


@dataclass
class ObservedScanBundle:
    selected_introns: pd.DataFrame
    thresholds: pd.DataFrame
    results: pd.DataFrame
    hits: pd.DataFrame


def load_fasta_to_dict(
    fasta_path: str | Path,
    ids_keep: set[str] | None = None,
    id_func=lambda x: x,
) -> dict[str, str]:
    seqs: dict[str, str] = {}
    with Path(fasta_path).open("r") as handle:
        current_id: str | None = None
        current_seq: list[str] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(current_seq).upper()
                    if ids_keep is None or current_id in ids_keep:
                        seqs[current_id] = seq
                current_id = id_func(line[1:].split()[0])
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            seq = "".join(current_seq).upper()
            if ids_keep is None or current_id in ids_keep:
                seqs[current_id] = seq
    return seqs


def clean_sequence(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    return "".join(base for base in seq if base in BASE_TO_INDEX)


def encode_sequence(seq: str) -> np.ndarray:
    out = np.full(len(seq), -1, dtype=np.int8)
    for i, base in enumerate(seq.upper()):
        out[i] = BASE_TO_INDEX.get(base, -1)
    return out


def gc_fraction_from_encoded(encoded: np.ndarray) -> float:
    valid = encoded >= 0
    if not valid.any():
        return np.nan
    gc = np.isin(encoded[valid], [1, 2]).sum()
    return float(gc / valid.sum())


def base_frequencies_from_encoded(encoded: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    valid = encoded[encoded >= 0]
    counts = np.bincount(valid, minlength=4).astype(float) + pseudocount
    return counts / counts.sum()


def parse_meme_minimal(path: str | Path) -> dict[str, MemeMotif]:
    lines = Path(path).read_text().splitlines()
    background = np.full(4, 0.25, dtype=float)
    bg_pattern = re.compile(r"([ACGT])\s+([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)")

    for i, line in enumerate(lines):
        if line.startswith("Background letter frequencies"):
            for follow in lines[i + 1 : i + 4]:
                hits = bg_pattern.findall(follow)
                if hits:
                    background = np.array(
                        [float(dict(hits).get(base, 0.25)) for base in DNA_ALPHABET],
                        dtype=float,
                    )
                    background = background / background.sum()
                    break
            break

    motifs: dict[str, MemeMotif] = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("MOTIF "):
            i += 1
            continue

        parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"Could not parse motif name from line: {line}")
        name = parts[1]

        i += 1
        while i < len(lines) and "letter-probability matrix:" not in lines[i]:
            i += 1
        if i >= len(lines):
            raise ValueError(f"Missing matrix block for motif {name}")

        matrix_line = lines[i]
        width_match = re.search(r"\bw=\s*(\d+)", matrix_line)
        if width_match is None:
            raise ValueError(f"Could not parse width from line: {matrix_line}")
        width = int(width_match.group(1))

        matrix_rows: list[list[float]] = []
        for _ in range(width):
            i += 1
            if i >= len(lines):
                raise ValueError(f"Unexpected end of file inside motif {name}")
            row = lines[i].strip()
            values = [float(value) for value in row.split()]
            if len(values) != 4:
                raise ValueError(f"Expected 4 columns in motif {name}, got: {row}")
            matrix_rows.append(values)

        matrix = np.asarray(matrix_rows, dtype=float)
        row_sums = matrix.sum(axis=1, keepdims=True)
        matrix = np.divide(matrix, row_sums, where=row_sums != 0)
        motifs[name] = MemeMotif(name=name, matrix=matrix, background=background.copy())
        i += 1

    if not motifs:
        raise ValueError(f"No motifs were parsed from {path}")
    return motifs


def motif_to_log_odds(
    motif: MemeMotif,
    background: np.ndarray | None = None,
    pseudocount: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    bg = motif.background if background is None else np.asarray(background, dtype=float)
    bg = bg / bg.sum()
    pwm = np.clip(motif.matrix, pseudocount, 1.0)
    fwd = np.log2(pwm / bg[np.newaxis, :])
    rev = fwd[::-1][:, COMPLEMENT_INDEX]
    return fwd, rev


def scan_encoded_sequence(
    encoded: np.ndarray,
    fwd_log_odds: np.ndarray,
    rev_log_odds: np.ndarray,
    threshold: float | None = None,
) -> dict[str, Any]:
    width = int(fwd_log_odds.shape[0])
    if encoded.size < width:
        return {
            "window_count": 0,
            "max_score": np.nan,
            "hit_count": 0,
            "hit_density": np.nan,
            "hit_positions": np.array([], dtype=int),
            "hit_scores": np.array([], dtype=float),
            "hit_strands": np.array([], dtype=object),
            "best_scores": np.array([], dtype=float),
        }

    windows = sliding_window_view(encoded, width)
    valid = (windows >= 0).all(axis=1)
    n_windows = int(windows.shape[0])

    fwd_scores = np.full(n_windows, -np.inf, dtype=float)
    rev_scores = np.full(n_windows, -np.inf, dtype=float)
    if valid.any():
        valid_windows = windows[valid]
        pos = np.arange(width)[:, np.newaxis]
        fwd_scores[valid] = fwd_log_odds[pos, valid_windows.T].sum(axis=0)
        rev_scores[valid] = rev_log_odds[pos, valid_windows.T].sum(axis=0)

    best_scores = np.maximum(fwd_scores, rev_scores)
    strands = np.where(fwd_scores >= rev_scores, "+", "-")
    finite_mask = np.isfinite(best_scores)
    if threshold is None:
        hit_mask = finite_mask
    else:
        hit_mask = finite_mask & (best_scores >= threshold)

    hit_positions = np.flatnonzero(hit_mask)
    hit_scores = best_scores[hit_mask]
    hit_strands = strands[hit_mask]

    return {
        "window_count": n_windows,
        "max_score": float(np.nanmax(best_scores[finite_mask])) if finite_mask.any() else np.nan,
        "hit_count": int(hit_mask.sum()),
        "hit_density": float(hit_mask.sum() / n_windows) if n_windows else np.nan,
        "hit_positions": hit_positions.astype(int),
        "hit_scores": hit_scores.astype(float),
        "hit_strands": hit_strands.astype(object),
        "best_scores": best_scores.astype(float),
    }


def estimate_hit_threshold(
    motif: MemeMotif,
    background: np.ndarray,
    rng: np.random.Generator,
    n_samples: int = 50_000,
    window_fpr: float = 0.01,
) -> float:
    fwd_log_odds, rev_log_odds = motif_to_log_odds(motif, background=background)
    windows = rng.choice(4, size=(n_samples, motif.width), p=background, shuffle=False)
    pos = np.arange(motif.width)[:, np.newaxis]
    fwd_scores = fwd_log_odds[pos, windows.T].sum(axis=0)
    rev_scores = rev_log_odds[pos, windows.T].sum(axis=0)
    best_scores = np.maximum(fwd_scores, rev_scores)
    return float(np.quantile(best_scores, 1.0 - window_fpr))


def simulate_null_distribution(
    encoded: np.ndarray,
    motif: MemeMotif,
    threshold: float,
    rng: np.random.Generator,
    n_null: int = 500,
    null_mode: str = "shuffle",
) -> pd.DataFrame:
    width = motif.width
    if encoded.size < width:
        return pd.DataFrame(columns=["iteration", "max_score", "hit_count", "hit_density"])

    seq_freqs = base_frequencies_from_encoded(encoded)
    fwd_log_odds, rev_log_odds = motif_to_log_odds(motif, background=motif.background)
    rows: list[dict[str, float | int]] = []

    for iteration in range(n_null):
        if null_mode == "shuffle":
            null_encoded = rng.permutation(encoded)
        elif null_mode == "iid":
            null_encoded = rng.choice(4, size=encoded.size, p=seq_freqs).astype(np.int8)
        else:
            raise ValueError(f"Unsupported null_mode: {null_mode}")

        scored = scan_encoded_sequence(null_encoded, fwd_log_odds, rev_log_odds, threshold=threshold)
        rows.append(
            {
                "iteration": iteration,
                "max_score": scored["max_score"],
                "hit_count": scored["hit_count"],
                "hit_density": scored["hit_density"],
            }
        )

    return pd.DataFrame(rows)


def empirical_pvalue(observed: float, null_values: pd.Series | np.ndarray) -> float:
    null = np.asarray(pd.Series(null_values).dropna(), dtype=float)
    if null.size == 0 or not np.isfinite(observed):
        return np.nan
    return float((np.sum(null >= observed) + 1) / (null.size + 1))


def z_from_null(observed: float, null_values: pd.Series | np.ndarray) -> float:
    null = np.asarray(pd.Series(null_values).dropna(), dtype=float)
    if null.size == 0 or not np.isfinite(observed):
        return np.nan
    sd = null.std(ddof=1)
    if sd == 0:
        return np.nan
    return float((observed - null.mean()) / sd)


def benjamini_hochberg(pvalues: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    ranks = np.arange(1, ranked.size + 1)
    adjusted = ranked * ranked.size / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    q_valid = np.empty_like(ranked)
    q_valid[order] = adjusted
    q[valid] = q_valid
    return q


def load_metadata_table(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "EVENT" not in df.columns and df.index.name == "EVENT":
        df = df.reset_index()
    return df


def build_sequence_table(
    fasta_path: str | Path,
    metadata_path: str | Path | None = None,
    seq_id_col: str = "EVENT",
    fasta_id_func=lambda x: x,
) -> pd.DataFrame:
    seq_dict = load_fasta_to_dict(fasta_path, id_func=fasta_id_func)
    rows: list[dict[str, Any]] = []
    for seq_id, raw_seq in seq_dict.items():
        seq = clean_sequence(raw_seq)
        encoded = encode_sequence(seq)
        rows.append(
            {
                seq_id_col: seq_id,
                "sequence": seq,
                "sequence_length": int(encoded.size),
                "gc_content": gc_fraction_from_encoded(encoded),
            }
        )

    seq_df = pd.DataFrame(rows)
    if metadata_path is None or str(metadata_path).strip() == "":
        return seq_df.sort_values(seq_id_col).reset_index(drop=True)

    meta = load_metadata_table(metadata_path)
    if seq_id_col not in meta.columns:
        raise ValueError(f"{seq_id_col} is not present in metadata columns: {sorted(meta.columns)}")

    merged = meta.merge(seq_df, on=seq_id_col, how="inner")
    if "Intron_Position" in merged.columns:
        merged = merged.sort_values(["Intron_Position", seq_id_col], na_position="last")
    else:
        merged = merged.sort_values(seq_id_col)
    return merged.reset_index(drop=True)


def subset_sequence_table(
    table: pd.DataFrame,
    group_col: str | None = None,
    group_value: str | None = None,
) -> pd.DataFrame:
    if group_col is None or group_col == "<all>" or not str(group_value or "").strip():
        return table.copy()

    mask = table[group_col].astype(str) == str(group_value)
    out = table.loc[mask].copy()
    if "Intron_Position" in out.columns:
        out = out.sort_values(["Intron_Position", table.columns[0]], na_position="last")
    return out.reset_index(drop=True)


def scan_sequences_observed_only(
    sequence_table: pd.DataFrame,
    motifs: dict[str, MemeMotif],
    motif_names: list[str],
    seq_id_col: str = "EVENT",
    window_fpr: float = 0.01,
    threshold_samples: int = 50_000,
    random_seed: int = 7,
) -> ObservedScanBundle:
    if sequence_table.empty:
        raise ValueError("No sequences were selected for analysis.")
    if not motif_names:
        raise ValueError("Please select at least one motif.")

    required = {seq_id_col, "sequence"}
    missing = required.difference(sequence_table.columns)
    if missing:
        raise ValueError(f"Missing required sequence table columns: {sorted(missing)}")

    encoded_map = {
        str(row[seq_id_col]): encode_sequence(str(row["sequence"]))
        for _, row in sequence_table.iterrows()
    }
    pooled_encoded = np.concatenate([encoded for encoded in encoded_map.values() if encoded.size > 0])
    pooled_background = base_frequencies_from_encoded(pooled_encoded)

    rng = np.random.default_rng(random_seed)
    threshold_rows: list[dict[str, Any]] = []
    thresholds: dict[str, float] = {}
    for motif_name in motif_names:
        motif = motifs[motif_name]
        threshold = estimate_hit_threshold(
            motif=motif,
            background=pooled_background,
            rng=np.random.default_rng(rng.integers(0, 2**32 - 1)),
            n_samples=threshold_samples,
            window_fpr=window_fpr,
        )
        thresholds[motif_name] = threshold
        threshold_rows.append(
            {
                "motif": motif_name,
                "width": motif.width,
                "window_fpr": window_fpr,
                "score_threshold": threshold,
            }
        )

    result_rows: list[dict[str, Any]] = []
    hit_rows: list[dict[str, Any]] = []
    for _, row in sequence_table.iterrows():
        seq_id = str(row[seq_id_col])
        encoded = encoded_map[seq_id]
        if encoded.size == 0:
            continue

        for motif_name in motif_names:
            motif = motifs[motif_name]
            fwd_log_odds, rev_log_odds = motif_to_log_odds(motif, background=motif.background)
            observed = scan_encoded_sequence(
                encoded=encoded,
                fwd_log_odds=fwd_log_odds,
                rev_log_odds=rev_log_odds,
                threshold=thresholds[motif_name],
            )

            result_row = {
                "motif": motif_name,
                seq_id_col: seq_id,
                "sequence_length": int(encoded.size),
                "window_count": observed["window_count"],
                "gc_content": row.get("gc_content", np.nan),
                "hit_count": observed["hit_count"],
                "hit_density": observed["hit_density"],
                "max_score": observed["max_score"],
                "score_threshold": thresholds[motif_name],
            }
            for col in sequence_table.columns:
                if col in result_row or col == "sequence":
                    continue
                result_row[col] = row[col]
            result_rows.append(result_row)

            scaled_len = max(int(encoded.size - motif.width + 1), 1)
            for pos, score, strand in zip(
                observed["hit_positions"],
                observed["hit_scores"],
                observed["hit_strands"],
            ):
                hit_rows.append(
                    {
                        "motif": motif_name,
                        seq_id_col: seq_id,
                        "position": int(pos),
                        "relative_position": float(pos / scaled_len),
                        "score": float(score),
                        "strand": strand,
                    }
                )

    results = pd.DataFrame(result_rows)
    if results.empty:
        raise ValueError("The selected sequences were shorter than all selected motifs.")

    sort_cols = [col for col in ["motif", "Intron_Position", seq_id_col] if col in results.columns]
    if sort_cols:
        results = results.sort_values(sort_cols).reset_index(drop=True)

    return ObservedScanBundle(
        selected_introns=sequence_table.copy(),
        thresholds=pd.DataFrame(threshold_rows),
        results=results,
        hits=pd.DataFrame(hit_rows),
    )


def analyze_sequences(
    sequence_table: pd.DataFrame,
    motifs: dict[str, MemeMotif],
    motif_names: list[str],
    seq_id_col: str = "EVENT",
    n_null: int = 500,
    null_mode: str = "shuffle",
    window_fpr: float = 0.01,
    threshold_samples: int = 50_000,
    random_seed: int = 7,
    pairwise_stat: str = "hit_count",
) -> AnalysisBundle:
    if sequence_table.empty:
        raise ValueError("No sequences were selected for analysis.")
    if not motif_names:
        raise ValueError("Please select at least one motif.")

    required = {seq_id_col, "sequence"}
    missing = required.difference(sequence_table.columns)
    if missing:
        raise ValueError(f"Missing required sequence table columns: {sorted(missing)}")

    encoded_map = {
        str(row[seq_id_col]): encode_sequence(str(row["sequence"]))
        for _, row in sequence_table.iterrows()
    }
    pooled_encoded = np.concatenate([encoded for encoded in encoded_map.values() if encoded.size > 0])
    pooled_background = base_frequencies_from_encoded(pooled_encoded)

    rng = np.random.default_rng(random_seed)
    threshold_rows: list[dict[str, Any]] = []
    thresholds: dict[str, float] = {}
    for motif_name in motif_names:
        motif = motifs[motif_name]
        threshold = estimate_hit_threshold(
            motif=motif,
            background=pooled_background,
            rng=np.random.default_rng(rng.integers(0, 2**32 - 1)),
            n_samples=threshold_samples,
            window_fpr=window_fpr,
        )
        thresholds[motif_name] = threshold
        threshold_rows.append(
            {
                "motif": motif_name,
                "width": motif.width,
                "window_fpr": window_fpr,
                "score_threshold": threshold,
            }
        )

    result_rows: list[dict[str, Any]] = []
    hit_rows: list[dict[str, Any]] = []
    null_distributions: dict[tuple[str, str], pd.DataFrame] = {}

    for _, row in sequence_table.iterrows():
        seq_id = str(row[seq_id_col])
        encoded = encoded_map[seq_id]
        if encoded.size == 0:
            continue

        for motif_name in motif_names:
            motif = motifs[motif_name]
            threshold = thresholds[motif_name]
            fwd_log_odds, rev_log_odds = motif_to_log_odds(motif, background=motif.background)
            observed = scan_encoded_sequence(encoded, fwd_log_odds, rev_log_odds, threshold=threshold)
            null_df = simulate_null_distribution(
                encoded=encoded,
                motif=motif,
                threshold=threshold,
                rng=np.random.default_rng(rng.integers(0, 2**32 - 1)),
                n_null=n_null,
                null_mode=null_mode,
            )
            null_distributions[(motif_name, seq_id)] = null_df

            count_mean = null_df["hit_count"].mean() if not null_df.empty else np.nan
            max_mean = null_df["max_score"].mean() if not null_df.empty else np.nan

            result_row = {
                "motif": motif_name,
                seq_id_col: seq_id,
                "sequence_length": int(encoded.size),
                "window_count": observed["window_count"],
                "gc_content": row.get("gc_content", np.nan),
                "hit_count": observed["hit_count"],
                "hit_density": observed["hit_density"],
                "max_score": observed["max_score"],
                "null_mean_hit_count": count_mean,
                "null_mean_max_score": max_mean,
                "hit_count_excess": observed["hit_count"] - count_mean if np.isfinite(count_mean) else np.nan,
                "max_score_excess": observed["max_score"] - max_mean if np.isfinite(max_mean) else np.nan,
                "hit_count_pvalue": empirical_pvalue(observed["hit_count"], null_df["hit_count"]),
                "max_score_pvalue": empirical_pvalue(observed["max_score"], null_df["max_score"]),
                "hit_density_pvalue": empirical_pvalue(observed["hit_density"], null_df["hit_density"]),
                "hit_count_z": z_from_null(observed["hit_count"], null_df["hit_count"]),
                "max_score_z": z_from_null(observed["max_score"], null_df["max_score"]),
                "hit_density_z": z_from_null(observed["hit_density"], null_df["hit_density"]),
                "score_threshold": threshold,
                "null_mode": null_mode,
                "n_null": n_null,
            }

            for col in sequence_table.columns:
                if col in result_row or col == "sequence":
                    continue
                result_row[col] = row[col]
            result_rows.append(result_row)

            seq_len = max(int(encoded.size - motif.width + 1), 1)
            for pos, score, strand in zip(
                observed["hit_positions"],
                observed["hit_scores"],
                observed["hit_strands"],
            ):
                hit_rows.append(
                    {
                        "motif": motif_name,
                        seq_id_col: seq_id,
                        "position": int(pos),
                        "relative_position": float(pos / seq_len),
                        "score": float(score),
                        "strand": strand,
                    }
                )

    results = pd.DataFrame(result_rows)
    if results.empty:
        raise ValueError("The selected sequences were shorter than all selected motifs.")

    results["hit_count_qvalue"] = benjamini_hochberg(results["hit_count_pvalue"])
    results["max_score_qvalue"] = benjamini_hochberg(results["max_score_pvalue"])
    results["hit_density_qvalue"] = benjamini_hochberg(results["hit_density_pvalue"])

    sort_cols = [col for col in ["motif", "Intron_Position", seq_id_col] if col in results.columns]
    if sort_cols:
        results = results.sort_values(sort_cols).reset_index(drop=True)

    pairwise = pairwise_intron_comparisons(
        results=results,
        null_distributions=null_distributions,
        seq_id_col=seq_id_col,
        stat=pairwise_stat,
    )

    return AnalysisBundle(
        selected_introns=sequence_table.copy(),
        thresholds=pd.DataFrame(threshold_rows),
        results=results,
        pairwise=pairwise,
        hits=pd.DataFrame(hit_rows),
        null_distributions=null_distributions,
    )


def pairwise_intron_comparisons(
    results: pd.DataFrame,
    null_distributions: dict[tuple[str, str], pd.DataFrame],
    seq_id_col: str = "EVENT",
    stat: str = "hit_count",
) -> pd.DataFrame:
    if results.empty or results[seq_id_col].nunique() < 2:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for motif_name, motif_df in results.groupby("motif", sort=False):
        motif_df = motif_df.reset_index(drop=True)
        for i, j in combinations(range(len(motif_df)), 2):
            left = motif_df.iloc[i]
            right = motif_df.iloc[j]
            left_id = str(left[seq_id_col])
            right_id = str(right[seq_id_col])
            left_null = null_distributions[(motif_name, left_id)][stat].dropna().to_numpy()
            right_null = null_distributions[(motif_name, right_id)][stat].dropna().to_numpy()
            n = min(left_null.size, right_null.size)
            if n == 0:
                continue

            diff_null = left_null[:n] - right_null[:n]
            diff_obs = float(left[stat] - right[stat])
            diff_sd = diff_null.std(ddof=1)
            z_score = np.nan if diff_sd == 0 else float((diff_obs - diff_null.mean()) / diff_sd)
            p_two_sided = float((np.sum(np.abs(diff_null) >= abs(diff_obs)) + 1) / (n + 1))
            higher_id = left_id if diff_obs >= 0 else right_id

            rows.append(
                {
                    "motif": motif_name,
                    "statistic": stat,
                    "left_intron": left_id,
                    "right_intron": right_id,
                    "observed_difference": diff_obs,
                    "null_mean_difference": float(diff_null.mean()),
                    "difference_z": z_score,
                    "pairwise_pvalue": p_two_sided,
                    "higher_intron": higher_id,
                }
            )

    pairwise = pd.DataFrame(rows)
    if pairwise.empty:
        return pairwise
    pairwise["pairwise_qvalue"] = benjamini_hochberg(pairwise["pairwise_pvalue"])
    return pairwise.sort_values(["motif", "pairwise_qvalue", "pairwise_pvalue"]).reset_index(drop=True)


def assign_length_bins(
    lengths: pd.Series | np.ndarray,
    n_bins: int = 5,
) -> pd.Series:
    length_series = pd.Series(lengths).astype(float)
    out = pd.Series(pd.NA, index=length_series.index, dtype="object")
    valid = length_series.dropna()
    if valid.empty:
        return out

    if n_bins <= 1 or valid.nunique() <= 1:
        out.loc[valid.index] = "bin_1"
        return out

    bins = int(min(n_bins, valid.nunique(), valid.shape[0]))
    ranked = valid.rank(method="first")
    codes = pd.qcut(ranked, q=bins, labels=False, duplicates="drop")
    out.loc[valid.index] = codes.map(lambda x: f"bin_{int(x) + 1}")
    return out


def length_matched_group_comparison(
    results: pd.DataFrame,
    group_col: str,
    group_a: str,
    group_b: str,
    value_col: str = "hit_density_z",
    seq_id_col: str = "EVENT",
    motif_col: str = "motif",
    length_col: str = "sequence_length",
    n_length_bins: int = 5,
    n_resamples: int = 1000,
    n_label_permutations: int = 1,
    random_seed: int = 7,
) -> GroupComparisonBundle:
    required = {group_col, value_col, motif_col, length_col, seq_id_col}
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"Missing required columns for group comparison: {sorted(missing)}")

    group_a = str(group_a)
    group_b = str(group_b)
    comp = results.copy()
    comp["_comparison_group"] = comp[group_col].astype(str)
    comp = comp[comp["_comparison_group"].isin([group_a, group_b])].copy()
    if comp.empty:
        raise ValueError(f"No rows found for {group_col} in [{group_a!r}, {group_b!r}].")

    rng = np.random.default_rng(random_seed)
    summary_rows: list[dict[str, Any]] = []
    resample_rows: list[dict[str, Any]] = []

    for motif_name, sub in comp.groupby(motif_col, sort=False):
        sub = sub.dropna(subset=[value_col, length_col]).copy()
        if sub.empty or sub["_comparison_group"].nunique() < 2:
            continue

        sub["_length_bin"] = assign_length_bins(sub[length_col], n_bins=n_length_bins)
        sub = sub.dropna(subset=["_length_bin"]).copy()
        if sub.empty:
            continue

        counts = sub.groupby(["_length_bin", "_comparison_group"]).size().unstack(fill_value=0)
        if group_a not in counts.columns or group_b not in counts.columns:
            continue

        usable_bins = [
            bin_name
            for bin_name in counts.index.tolist()
            if counts.loc[bin_name, group_a] > 0 and counts.loc[bin_name, group_b] > 0
        ]
        if not usable_bins:
            continue

        matched_per_bin = {
            bin_name: int(min(counts.loc[bin_name, group_a], counts.loc[bin_name, group_b]))
            for bin_name in usable_bins
        }
        matched_per_group = int(sum(matched_per_bin.values()))
        if matched_per_group == 0:
            continue

        obs_diffs: list[float] = []
        perm_diffs: list[float] = []

        for resample_idx in range(n_resamples):
            sampled_parts = []
            for bin_name in usable_bins:
                n_take = matched_per_bin[bin_name]
                for group_name in [group_a, group_b]:
                    pool = sub[(sub["_length_bin"] == bin_name) & (sub["_comparison_group"] == group_name)]
                    sampled = pool.sample(
                        n=n_take,
                        replace=False,
                        random_state=int(rng.integers(0, 2**32 - 1)),
                    )
                    sampled_parts.append(sampled)

            matched = pd.concat(sampled_parts, ignore_index=True)
            mean_a = matched.loc[matched["_comparison_group"] == group_a, value_col].mean()
            mean_b = matched.loc[matched["_comparison_group"] == group_b, value_col].mean()
            obs_diff = float(mean_a - mean_b)
            obs_diffs.append(obs_diff)
            resample_rows.append(
                {
                    motif_col: motif_name,
                    "draw_type": "observed",
                    "resample_iteration": resample_idx,
                    "permutation_iteration": pd.NA,
                    "difference": obs_diff,
                    "group_a_mean": float(mean_a),
                    "group_b_mean": float(mean_b),
                    "matched_introns_per_group": matched_per_group,
                }
            )

            matched_bins = matched["_length_bin"].to_numpy(dtype=object)
            for perm_idx in range(n_label_permutations):
                perm_labels = matched["_comparison_group"].to_numpy(dtype=object).copy()
                for bin_name in usable_bins:
                    mask = matched_bins == bin_name
                    perm_labels[mask] = rng.permutation(perm_labels[mask])

                perm_mean_a = matched.loc[perm_labels == group_a, value_col].mean()
                perm_mean_b = matched.loc[perm_labels == group_b, value_col].mean()
                perm_diff = float(perm_mean_a - perm_mean_b)
                perm_diffs.append(perm_diff)
                resample_rows.append(
                    {
                        motif_col: motif_name,
                        "draw_type": "permuted",
                        "resample_iteration": resample_idx,
                        "permutation_iteration": perm_idx,
                        "difference": perm_diff,
                        "group_a_mean": float(perm_mean_a),
                        "group_b_mean": float(perm_mean_b),
                        "matched_introns_per_group": matched_per_group,
                    }
                )

        obs_array = np.asarray(obs_diffs, dtype=float)
        perm_array = np.asarray(perm_diffs, dtype=float)
        obs_mean = float(obs_array.mean()) if obs_array.size else np.nan
        perm_mean = float(perm_array.mean()) if perm_array.size else np.nan
        perm_sd = float(perm_array.std(ddof=1)) if perm_array.size > 1 else np.nan
        pvalue = (
            float((np.sum(np.abs(perm_array) >= abs(obs_mean)) + 1) / (perm_array.size + 1))
            if perm_array.size
            else np.nan
        )
        z_score = (
            float((obs_mean - perm_mean) / perm_sd)
            if np.isfinite(perm_sd) and perm_sd > 0
            else np.nan
        )

        summary_rows.append(
            {
                motif_col: motif_name,
                "group_col": group_col,
                "group_a": group_a,
                "group_b": group_b,
                "value_col": value_col,
                "group_a_total_introns": int((sub["_comparison_group"] == group_a).sum()),
                "group_b_total_introns": int((sub["_comparison_group"] == group_b).sum()),
                "usable_length_bins": len(usable_bins),
                "matched_introns_per_group": matched_per_group,
                "matched_total_introns": 2 * matched_per_group,
                "observed_mean_difference": obs_mean,
                "observed_median_difference": float(np.median(obs_array)),
                "observed_ci_low": float(np.quantile(obs_array, 0.025)),
                "observed_ci_high": float(np.quantile(obs_array, 0.975)),
                "permuted_mean_difference": perm_mean,
                "permuted_sd": perm_sd,
                "difference_z": z_score,
                "pvalue": pvalue,
                "higher_group": group_a if obs_mean >= 0 else group_b,
                "length_bins": ", ".join(f"{bin_name}:{matched_per_bin[bin_name]}" for bin_name in usable_bins),
            }
        )

    summary = pd.DataFrame(summary_rows)
    resamples = pd.DataFrame(resample_rows)
    if not summary.empty:
        summary["qvalue"] = benjamini_hochberg(summary["pvalue"])
        summary = summary.sort_values([motif_col, "qvalue", "pvalue"]).reset_index(drop=True)
    return GroupComparisonBundle(summary=summary, resamples=resamples)


def plot_enrichment_heatmap(
    results: pd.DataFrame,
    seq_id_col: str = "EVENT",
    value_col: str = "hit_count_z",
    cmap: str = "viridis",
) -> plt.Figure:
    plot_df = results.copy()
    order_col = "Intron_Position" if "Intron_Position" in plot_df.columns else seq_id_col
    row_order = (
        plot_df[[seq_id_col, order_col]]
        .drop_duplicates()
        .sort_values(order_col, na_position="last")[seq_id_col]
        .tolist()
    )
    pivot = plot_df.pivot(index=seq_id_col, columns="motif", values=value_col).reindex(row_order)

    fig_width = max(6, 1.2 * max(1, pivot.shape[1]))
    fig_height = max(3, 0.6 * max(1, pivot.shape[0]) + 1.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(pivot, cmap=cmap, annot=True, fmt=".2f", linewidths=0.5, ax=ax)
    ax.set_title(f"Intron-by-motif enrichment heatmap ({value_col})")
    ax.set_xlabel("Motif")
    ax.set_ylabel("Intron")
    fig.tight_layout()
    return fig


def plot_hit_positions(
    hits: pd.DataFrame,
    selected_introns: pd.DataFrame,
    motif_name: str,
    seq_id_col: str = "EVENT",
) -> plt.Figure | None:
    motif_hits = hits[hits["motif"] == motif_name].copy()
    if motif_hits.empty:
        return None

    order_col = "Intron_Position" if "Intron_Position" in selected_introns.columns else seq_id_col
    ordered_introns = (
        selected_introns[[seq_id_col, order_col]]
        .drop_duplicates()
        .sort_values(order_col, na_position="last")[seq_id_col]
        .astype(str)
        .tolist()
    )
    y_positions = {seq_id: i for i, seq_id in enumerate(ordered_introns)}
    motif_hits["y"] = motif_hits[seq_id_col].astype(str).map(y_positions)

    fig_height = max(3, 0.45 * max(1, len(ordered_introns)) + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    ax.scatter(
        motif_hits["relative_position"],
        motif_hits["y"],
        c=motif_hits["score"],
        cmap="magma",
        s=60,
        edgecolor="black",
        linewidth=0.3,
    )
    ax.set_yticks(range(len(ordered_introns)))
    ax.set_yticklabels(ordered_introns)
    ax.set_xlabel("Relative intron position")
    ax.set_ylabel("Intron")
    ax.set_xlim(-0.02, 1.02)
    ax.set_title(f"Hit positions for {motif_name}")
    fig.tight_layout()
    return fig


def export_analysis(
    bundle: AnalysisBundle,
    out_dir: str | Path,
    prefix: str,
) -> dict[str, Path]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    outputs = {
        "selected_introns": out_path / f"{prefix}_selected_introns.csv",
        "thresholds": out_path / f"{prefix}_thresholds.csv",
        "results": out_path / f"{prefix}_results.csv",
        "pairwise": out_path / f"{prefix}_pairwise.csv",
        "hits": out_path / f"{prefix}_hits.csv",
    }

    bundle.selected_introns.to_csv(outputs["selected_introns"], index=False)
    bundle.thresholds.to_csv(outputs["thresholds"], index=False)
    bundle.results.to_csv(outputs["results"], index=False)
    bundle.hits.to_csv(outputs["hits"], index=False)
    if not bundle.pairwise.empty:
        bundle.pairwise.to_csv(outputs["pairwise"], index=False)

    return outputs


def _infer_group_columns(meta: pd.DataFrame | None) -> list[str]:
    if meta is None or meta.empty:
        return ["<all>"]
    preferred = [col for col in ["GENE", "gene_name", "transcript_id", "EVENT"] if col in meta.columns]
    rest = [col for col in meta.columns if col not in preferred]
    return ["<all>"] + preferred + rest


def _sanitize_prefix(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "all_introns"


def build_notebook_app(
    repo_root: str | Path,
    default_fasta: str | Path,
    default_metadata: str | Path | None,
    default_meme: str | Path,
    default_output_dir: str | Path,
):
    import ipywidgets as widgets

    repo_root = Path(repo_root)
    state: dict[str, Any] = {"metadata": None, "motifs": None}

    fasta_path = widgets.Text(value=str(default_fasta), description="FASTA", layout=widgets.Layout(width="100%"))
    metadata_path = widgets.Text(
        value="" if default_metadata is None else str(default_metadata),
        description="Metadata",
        layout=widgets.Layout(width="100%"),
    )
    meme_path = widgets.Text(value=str(default_meme), description="MEME", layout=widgets.Layout(width="100%"))
    output_dir = widgets.Text(
        value=str(default_output_dir),
        description="Output",
        layout=widgets.Layout(width="100%"),
    )

    seq_id_col = widgets.Dropdown(options=["EVENT"], value="EVENT", description="Seq ID")
    group_col = widgets.Dropdown(options=["<all>"], value="<all>", description="Group by")
    group_value = widgets.Text(value="", description="Group value", placeholder="e.g. AASDH")

    motif_filter = widgets.Text(value="", description="Motif filter", placeholder="substring match")
    motif_select = widgets.SelectMultiple(
        options=[],
        value=(),
        rows=10,
        description="Motifs",
        layout=widgets.Layout(width="100%"),
    )

    n_null = widgets.IntText(value=500, description="Null draws")
    threshold_samples = widgets.IntText(value=50_000, description="Thr draws")
    window_fpr = widgets.FloatText(value=0.01, description="Window FPR")
    random_seed = widgets.IntText(value=7, description="Seed")
    null_mode = widgets.Dropdown(options=["shuffle", "iid"], value="shuffle", description="Null model")
    pairwise_stat = widgets.Dropdown(
        options=["hit_count", "hit_density", "max_score"],
        value="hit_count",
        description="Compare by",
    )
    save_outputs = widgets.Checkbox(value=True, description="Write CSV outputs")

    refresh_button = widgets.Button(description="Refresh inputs", button_style="info")
    run_button = widgets.Button(description="Run analysis", button_style="success")
    status_out = widgets.Output()
    results_out = widgets.Output()

    def refresh_motif_options(*_args):
        motifs = state.get("motifs") or {}
        current = set(motif_select.value)
        names = list(motifs.keys())
        query = motif_filter.value.strip().lower()
        if query:
            names = [name for name in names if query in name.lower()]
        motif_select.options = names
        kept = tuple(name for name in names if name in current)
        if kept:
            motif_select.value = kept
        elif names:
            motif_select.value = tuple(names[: min(2, len(names))])

    def refresh_inputs(_button=None):
        with status_out:
            status_out.clear_output()
            try:
                meta = None
                if metadata_path.value.strip():
                    meta = load_metadata_table(metadata_path.value.strip())
                motifs = parse_meme_minimal(meme_path.value.strip())
                state["metadata"] = meta
                state["motifs"] = motifs

                if meta is not None and not meta.empty:
                    seq_opts = meta.columns.tolist()
                    seq_id_col.options = seq_opts
                    seq_id_col.value = "EVENT" if "EVENT" in seq_opts else seq_opts[0]
                else:
                    seq_id_col.options = ["EVENT"]
                    seq_id_col.value = "EVENT"

                group_col.options = _infer_group_columns(meta)
                if "GENE" in group_col.options:
                    group_col.value = "GENE"
                refresh_motif_options()

                print(f"Loaded {len(motifs):,} motifs from {meme_path.value}")
                if meta is not None:
                    print(f"Loaded metadata table with {len(meta):,} rows and {len(meta.columns)} columns")
                else:
                    print("No metadata table loaded. The notebook will analyze all FASTA entries as provided.")
            except Exception:
                traceback.print_exc()

    def run_analysis(_button=None):
        with results_out:
            results_out.clear_output()
            try:
                display(
                    Markdown(
                        "### Pipeline logic\n"
                        "1. Parse the MEME motifs from plain-text MEME minimal format.\n"
                        "2. Load intron FASTA sequences and optionally merge metadata.\n"
                        "3. For each selected motif, estimate a score threshold from pooled background windows.\n"
                        "4. For each intron, compare observed motif signal against length-matched random sequences.\n"
                        "5. Report empirical p-values, BH-adjusted q-values, and pairwise intron comparisons."
                    )
                )

                table = build_sequence_table(
                    fasta_path=fasta_path.value.strip(),
                    metadata_path=metadata_path.value.strip() or None,
                    seq_id_col=seq_id_col.value,
                )
                selected = subset_sequence_table(
                    table=table,
                    group_col=group_col.value,
                    group_value=group_value.value,
                )
                if selected.empty:
                    raise ValueError("The chosen group returned no introns.")

                chosen_motifs = list(motif_select.value)
                if not chosen_motifs:
                    raise ValueError("Please select at least one motif.")

                bundle = analyze_sequences(
                    sequence_table=selected,
                    motifs=state["motifs"],
                    motif_names=chosen_motifs,
                    seq_id_col=seq_id_col.value,
                    n_null=n_null.value,
                    null_mode=null_mode.value,
                    window_fpr=window_fpr.value,
                    threshold_samples=threshold_samples.value,
                    random_seed=random_seed.value,
                    pairwise_stat=pairwise_stat.value,
                )

                intron_cols = [col for col in [seq_id_col.value, "GENE", "Intron_Position", "sequence_length", "gc_content"] if col in bundle.selected_introns.columns]
                display(Markdown("### Selected introns"))
                display(bundle.selected_introns[intron_cols])

                display(Markdown("### Motif thresholds"))
                display(bundle.thresholds)

                display(Markdown("### Per-intron enrichment results"))
                result_cols = [
                    col
                    for col in [
                        "motif",
                        seq_id_col.value,
                        "GENE",
                        "Intron_Position",
                        "sequence_length",
                        "hit_count",
                        "null_mean_hit_count",
                        "hit_count_excess",
                        "hit_count_z",
                        "hit_count_pvalue",
                        "hit_count_qvalue",
                        "max_score",
                        "null_mean_max_score",
                        "max_score_z",
                        "max_score_pvalue",
                        "max_score_qvalue",
                    ]
                    if col in bundle.results.columns
                ]
                display(bundle.results[result_cols])

                if not bundle.pairwise.empty:
                    display(Markdown(f"### Pairwise intron comparisons ({pairwise_stat.value})"))
                    display(bundle.pairwise)

                heatmap = plot_enrichment_heatmap(bundle.results, seq_id_col=seq_id_col.value, value_col="hit_count_z")
                display(heatmap)
                plt.show()

                for motif_name in chosen_motifs:
                    fig = plot_hit_positions(bundle.hits, bundle.selected_introns, motif_name, seq_id_col=seq_id_col.value)
                    if fig is not None:
                        display(fig)
                        plt.show()

                if save_outputs.value:
                    prefix_group = group_value.value if group_value.value.strip() else "all_introns"
                    prefix = _sanitize_prefix(prefix_group)
                    outputs = export_analysis(bundle, output_dir.value.strip(), prefix)
                    display(Markdown("### Written outputs"))
                    display(pd.Series({name: str(path) for name, path in outputs.items()}))
            except Exception:
                traceback.print_exc()

    motif_filter.observe(refresh_motif_options, names="value")
    refresh_button.on_click(refresh_inputs)
    run_button.on_click(run_analysis)

    controls = widgets.VBox(
        [
            widgets.HTML("<h3>Intron motif enrichment notebook UI</h3>"),
            fasta_path,
            metadata_path,
            meme_path,
            output_dir,
            widgets.HBox([seq_id_col, group_col, group_value]),
            motif_filter,
            motif_select,
            widgets.HBox([n_null, threshold_samples, window_fpr]),
            widgets.HBox([random_seed, null_mode, pairwise_stat]),
            save_outputs,
            widgets.HBox([refresh_button, run_button]),
            status_out,
            results_out,
        ]
    )

    refresh_inputs()
    return controls
