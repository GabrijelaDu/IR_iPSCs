from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

BIN_SIZE = 10_000


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper()


def load_introns(path_csv: Path) -> pd.DataFrame:
    if not path_csv.exists():
        raise FileNotFoundError(f"Intron source table not found: {path_csv}")

    df = pd.read_csv(path_csv, low_memory=False)
    required = {"chrom", "start", "end", "strand"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in intron table: {sorted(missing)}")

    df = df.copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "strand"]).reset_index(drop=True)

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["strand"] = df["strand"].astype(str).str.strip()
    df = df[(df["strand"].isin({"+", "-"})) & (df["end"] >= df["start"])].reset_index(drop=True)

    if df.empty:
        raise ValueError("No valid intron rows left after cleaning.")

    df["chrom_norm"] = df["chrom"].map(normalize_chrom)
    # Convert 1-based inclusive to 0-based half-open for interval arithmetic.
    df["start0"] = (df["start"] - 1).clip(lower=0)
    df["end0"] = df["end"]
    return df


def load_narrowpeak(path_peak: Path) -> pd.DataFrame:
    if not path_peak.exists():
        raise FileNotFoundError(f"Peak file not found: {path_peak}")

    cols = [
        "chrom",
        "start0",
        "end0",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
    ]
    peaks = pd.read_csv(path_peak, sep="\t", header=None, names=cols, usecols=range(10), low_memory=False)

    peaks = peaks.copy()
    peaks["start0"] = pd.to_numeric(peaks["start0"], errors="coerce")
    peaks["end0"] = pd.to_numeric(peaks["end0"], errors="coerce")
    peaks["signalValue"] = pd.to_numeric(peaks["signalValue"], errors="coerce")
    peaks["qValue"] = pd.to_numeric(peaks["qValue"], errors="coerce")
    peaks = peaks.dropna(subset=["chrom", "start0", "end0", "signalValue", "qValue"]).reset_index(drop=True)

    peaks["start0"] = peaks["start0"].astype(int)
    peaks["end0"] = peaks["end0"].astype(int)
    peaks = peaks[peaks["end0"] > peaks["start0"]].reset_index(drop=True)

    if peaks.empty:
        raise ValueError(f"No valid peaks left after cleaning: {path_peak}")

    peaks["chrom_norm"] = peaks["chrom"].map(normalize_chrom)
    return peaks


def build_peak_bin_index(peaks: pd.DataFrame, bin_size: int) -> Dict[str, Dict[int, List[int]]]:
    index: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))

    for idx, row in peaks[["chrom_norm", "start0", "end0"]].iterrows():
        chrom = row["chrom_norm"]
        start0 = int(row["start0"])
        end0 = int(row["end0"])

        start_bin = start0 // bin_size
        end_bin = (end0 - 1) // bin_size

        bins = index[chrom]
        for b in range(start_bin, end_bin + 1):
            bins[b].append(int(idx))

    return index


def merged_length(intervals: Iterable[Tuple[int, int]]) -> int:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return 0

    total = 0
    cur_s, cur_e = sorted_intervals[0]
    for s, e in sorted_intervals[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    total += cur_e - cur_s
    return int(total)


def summarize_introns_with_peaks(
    introns: pd.DataFrame,
    peaks_fwd: pd.DataFrame,
    peaks_rev: pd.DataFrame,
    idx_fwd: Dict[str, Dict[int, List[int]]],
    idx_rev: Dict[str, Dict[int, List[int]]],
    bin_size: int,
) -> pd.DataFrame:
    peak_count = np.zeros(len(introns), dtype=np.int64)
    overlap_bp = np.zeros(len(introns), dtype=np.int64)
    signal_sum = np.zeros(len(introns), dtype=np.float64)
    qvalue_max = np.full(len(introns), np.nan)

    for i, row in introns[["chrom_norm", "start0", "end0", "strand"]].iterrows():
        chrom = row["chrom_norm"]
        r_start = int(row["start0"])
        r_end = int(row["end0"])
        strand = row["strand"]

        if strand == "+":
            peaks = peaks_fwd
            idx = idx_fwd
        else:
            peaks = peaks_rev
            idx = idx_rev

        bins = idx.get(chrom)
        if not bins:
            continue

        start_bin = r_start // bin_size
        end_bin = (r_end - 1) // bin_size

        candidates: List[int] = []
        for b in range(start_bin, end_bin + 1):
            candidates.extend(bins.get(b, []))

        if not candidates:
            continue

        seen = set()
        overlaps: List[Tuple[int, int]] = []
        qvals: List[float] = []
        sigvals: List[float] = []

        for p_idx in candidates:
            if p_idx in seen:
                continue
            seen.add(p_idx)

            p = peaks.iloc[p_idx]
            p_start = int(p["start0"])
            p_end = int(p["end0"])

            ov_s = max(r_start, p_start)
            ov_e = min(r_end, p_end)
            if ov_s >= ov_e:
                continue

            overlaps.append((ov_s, ov_e))
            qvals.append(float(p["qValue"]))
            sigvals.append(float(p["signalValue"]))

        if not overlaps:
            continue

        peak_count[i] = len(overlaps)
        overlap_bp[i] = merged_length(overlaps)
        signal_sum[i] = float(np.sum(sigvals))
        qvalue_max[i] = float(np.max(qvals))

    out = introns.drop(columns=["chrom_norm", "start0", "end0"]).copy()
    out["ssDRIP_peak_count_strand_matched"] = peak_count
    out["ssDRIP_overlap_bp_strand_matched"] = overlap_bp
    out["ssDRIP_signalValue_sum_strand_matched"] = signal_sum
    out["ssDRIP_qValue_max_strand_matched"] = qvalue_max
    out["ssDRIP_has_peak_strand_matched"] = peak_count > 0
    return out


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        description=(
            "Summarize strand-matched ssDRIP narrowPeak overlaps per intron. "
            "+ introns are overlapped with forward peaks, - introns with reverse peaks."
        )
    )
    parser.add_argument(
        "--introns",
        type=Path,
        default=project_root / "resubmission" / "data" / "all_vastdb_introns_source.csv",
        help="Path to intron CSV (must include chrom,start,end,strand).",
    )
    parser.add_argument(
        "--fwd-peaks",
        type=Path,
        default=project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "peaks_macs" / "batch2_fwd_peaks.narrowPeak",
        help="Path to forward-strand narrowPeak file.",
    )
    parser.add_argument(
        "--rev-peaks",
        type=Path,
        default=project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "peaks_macs" / "batch2_rev_peaks.narrowPeak",
        help="Path to reverse-strand narrowPeak file.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_root / "resubmission" / "data" / "all_vastdb_introns_ssDRIP_peaks_batch2.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=BIN_SIZE,
        help="Genomic bin size used for overlap acceleration (default: 10000).",
    )
    args = parser.parse_args()

    introns = load_introns(args.introns)
    peaks_fwd = load_narrowpeak(args.fwd_peaks)
    peaks_rev = load_narrowpeak(args.rev_peaks)

    idx_fwd = build_peak_bin_index(peaks_fwd, args.bin_size)
    idx_rev = build_peak_bin_index(peaks_rev, args.bin_size)

    print(f"Loaded {len(introns):,} introns from {args.introns}")
    print(f"Loaded {len(peaks_fwd):,} forward peaks from {args.fwd_peaks}")
    print(f"Loaded {len(peaks_rev):,} reverse peaks from {args.rev_peaks}")

    out_df = summarize_introns_with_peaks(
        introns=introns,
        peaks_fwd=peaks_fwd,
        peaks_rev=peaks_rev,
        idx_fwd=idx_fwd,
        idx_rev=idx_rev,
        bin_size=args.bin_size,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)

    n_peaks = int(out_df["ssDRIP_has_peak_strand_matched"].sum())
    print(f"Introns with strand-matched peak overlap: {n_peaks:,}/{len(out_df):,}")
    print(f"Wrote per-intron ssDRIP peak summary to: {args.out}")


if __name__ == "__main__":
    main()
