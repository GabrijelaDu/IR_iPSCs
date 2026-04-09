from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

BIN_SIZE = 10_000


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper()


def load_regions(regions_csv: Path) -> pd.DataFrame:
    if not regions_csv.exists():
        raise FileNotFoundError(f"Region table not found: {regions_csv}")

    df = pd.read_csv(regions_csv, low_memory=False)
    required = {"chrom", "start", "end"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in region table: {sorted(missing)}")

    df = df.copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end"]).reset_index(drop=True)

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["end"] >= df["start"]].reset_index(drop=True)

    if df.empty:
        raise ValueError("No valid region rows left after cleaning.")

    df["chrom_norm"] = df["chrom"].map(normalize_chrom)
    return df


def build_bin_index(regions: pd.DataFrame, bin_size: int) -> Dict[str, Dict[int, List[int]]]:
    index: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))

    for idx, row in regions[["chrom_norm", "start", "end"]].iterrows():
        chrom = row["chrom_norm"]
        start = int(row["start"])
        end = int(row["end"])

        start_bin = start // bin_size
        end_bin = end // bin_size

        bins = index[chrom]
        for b in range(start_bin, end_bin + 1):
            bins[b].append(int(idx))

    return index


def iter_cov_lines(cov_gz: Path) -> Iterable[str]:
    with gzip.open(cov_gz, "rt") as handle:
        for line in handle:
            yield line


def summarize_cov_per_region(
    regions: pd.DataFrame,
    index: Dict[str, Dict[int, List[int]]],
    cov_gz: Path,
    bin_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    meth_calls = np.zeros(len(regions), dtype=np.int64)
    unmeth_calls = np.zeros(len(regions), dtype=np.int64)
    cpg_sites = np.zeros(len(regions), dtype=np.int64)

    starts = regions["start"].to_numpy()
    ends = regions["end"].to_numpy()

    for line in iter_cov_lines(cov_gz):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            fields = line.strip().split()
            if len(fields) < 6:
                continue

        chrom = normalize_chrom(fields[0])
        try:
            pos = int(fields[1])
            methylated = int(fields[4])
            unmethylated = int(fields[5])
        except ValueError:
            continue

        bins = index.get(chrom)
        if not bins:
            continue

        bin_id = pos // bin_size
        candidates = bins.get(bin_id)
        if not candidates:
            continue

        for region_idx in candidates:
            if starts[region_idx] <= pos <= ends[region_idx]:
                meth_calls[region_idx] += methylated
                unmeth_calls[region_idx] += unmethylated
                cpg_sites[region_idx] += 1

    return meth_calls, unmeth_calls, cpg_sites


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        description=(
            "Summarize CpG methylation per region from a Bismark coverage file (.cov.gz). "
            "Outputs methylated/unmethylated counts, number of CpG sites hit per region, "
            "and weighted methylation percentage."
        )
    )
    parser.add_argument(
        "--regions",
        type=Path,
        default=project_root / "resubmission" / "data" / "all_vastdb_introns_source.csv",
        help="Path to region CSV (must include chrom,start,end).",
    )
    parser.add_argument(
        "--cov",
        type=Path,
        default=(
            project_root
            / "ipsc"
            / "data_raw"
            / "PRJNA608890"
            / "WGBS"
            / "SRR11185407"
            / "SRR11185407_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
        ),
        help="Path to Bismark coverage file (.cov.gz).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_root / "resubmission" / "data" / "all_vastdb_introns_source_cpg_methylation.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=BIN_SIZE,
        help="Genomic bin size used for overlap acceleration (default: 10000).",
    )
    args = parser.parse_args()

    regions = load_regions(args.regions)
    lookup = build_bin_index(regions, bin_size=args.bin_size)

    print(f"Loaded {len(regions):,} regions from {args.regions}")
    print(f"Using Bismark coverage file: {args.cov}")

    meth_calls, unmeth_calls, cpg_sites = summarize_cov_per_region(
        regions=regions,
        index=lookup,
        cov_gz=args.cov,
        bin_size=args.bin_size,
    )

    total_calls = meth_calls + unmeth_calls

    out_df = regions.drop(columns=["chrom_norm"]).copy()
    out_df["CpG_sites_in_region"] = cpg_sites
    out_df["CpG_methylated_calls"] = meth_calls
    out_df["CpG_unmethylated_calls"] = unmeth_calls
    out_df["CpG_total_calls"] = total_calls

    weighted_pct = np.where(total_calls > 0, 100.0 * meth_calls / total_calls, np.nan)
    out_df["CpG_methylation_pct_weighted"] = weighted_pct

    # Mean per-site methylation is approximated from aggregated calls and equals weighted_pct
    # for this table; the explicit column is included for downstream readability.
    out_df["CpG_methylation_pct"] = out_df["CpG_methylation_pct_weighted"]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)

    n_with_sites = int((cpg_sites > 0).sum())
    print(f"Regions with at least one CpG site: {n_with_sites:,}/{len(out_df):,}")
    print(
        f"Global calls: methylated={int(meth_calls.sum()):,}, "
        f"unmethylated={int(unmeth_calls.sum()):,}, "
        f"total={int(total_calls.sum()):,}"
    )
    print(f"Wrote per-region CpG methylation summary to: {args.out}")


if __name__ == "__main__":
    main()
