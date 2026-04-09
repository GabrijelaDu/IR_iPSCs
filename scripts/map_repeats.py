from pathlib import Path

import numpy as np
import pandas as pd


def merged_covered_bp(starts: np.ndarray, ends: np.ndarray) -> int:
    if len(starts) == 0:
        return 0
    order = np.argsort(starts, kind="mergesort")
    starts = starts[order]
    ends = ends[order]
    covered = 0
    cur_start, cur_end = starts[0], ends[0]
    for start_val, end_val in zip(starts[1:], ends[1:]):
        if start_val <= cur_end + 1:
            if end_val > cur_end:
                cur_end = end_val
        else:
            covered += cur_end - cur_start + 1
            cur_start, cur_end = start_val, end_val
    covered += cur_end - cur_start + 1
    return int(covered)


def load_vastdb_introns(path_vastdb: Path) -> pd.DataFrame:
    if not path_vastdb.exists():
        raise FileNotFoundError(
            f"VASTDB annotation file not found at {path_vastdb}. "
            "Please get one at https://vastdb.crg.eu/downloads/hg38/EVENT_INFO-hg38.tab.gz."
        )

    vast_annot = pd.read_csv(path_vastdb, sep="\t")
    vast_annot = vast_annot[lambda df: df["EVENT"].str.contains("HsaINT", na=False)]

    introns_source = (
        vast_annot[["GENE", "EVENT", "COORD_o", "FULL_CO"]]
        .assign(
            chrom=vast_annot["COORD_o"].str.split(":").str[0],
            start=vast_annot["COORD_o"].str.split(":").str[1].str.split("-").str[0].astype(int),
            end=vast_annot["COORD_o"].str.split(":").str[1].str.split("-").str[1].astype(int),
            strand=vast_annot["FULL_CO"].str.split(":").str[2],
        )
        .drop(columns=["COORD_o", "FULL_CO"], inplace=False)
        .reset_index(drop=True)
    )

    return introns_source


def load_repeatmasker(path_repeat: Path) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    if not path_repeat.exists():
        raise FileNotFoundError(f"RepeatMasker file not found at {path_repeat}")

    repeats = pd.read_csv(
        path_repeat,
        sep=r"\s+",
        skiprows=3,
        header=None,
        usecols=[4, 5, 6, 10],
        names=["chrom", "rep_start", "rep_end", "repeat_type"],
        compression="gzip",
        engine="python",
    )

    repeats["rep_start"] = pd.to_numeric(repeats["rep_start"], errors="coerce")
    repeats["rep_end"] = pd.to_numeric(repeats["rep_end"], errors="coerce")
    repeats = repeats.dropna(subset=["chrom", "rep_start", "rep_end"]).copy()
    repeats["rep_start"] = repeats["rep_start"].astype(int)
    repeats["rep_end"] = repeats["rep_end"].astype(int)

    rep_by_chr: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, group in repeats.groupby("chrom", sort=False):
        group = group.sort_values("rep_start", kind="mergesort")
        rep_by_chr[chrom] = (
            group["rep_start"].to_numpy(),
            group["rep_end"].to_numpy(),
            group["repeat_type"].astype(str).to_numpy(),
        )
    return rep_by_chr


def compute_repeat_scores(
    introns_source: pd.DataFrame,
    rep_by_chr: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> pd.DataFrame:
    repeat_pct = np.full(len(introns_source), np.nan, dtype=float)
    repeat_types = np.empty(len(introns_source), dtype=object)

    for i, (chrom, start, end) in enumerate(
        introns_source[["chrom", "start", "end"]].itertuples(index=False, name=None)
    ):
        data = rep_by_chr.get(chrom)
        if data is None:
            repeat_pct[i] = 0.0
            repeat_types[i] = np.nan
            continue

        rep_start, rep_end, rep_type = data
        right = np.searchsorted(rep_start, end, side="right")
        if right == 0:
            repeat_pct[i] = 0.0
            repeat_types[i] = np.nan
            continue

        candidates = rep_end[:right] >= start
        if not np.any(candidates):
            repeat_pct[i] = 0.0
            repeat_types[i] = np.nan
            continue

        ov_start = np.maximum(rep_start[:right][candidates], start)
        ov_end = np.minimum(rep_end[:right][candidates], end)
        valid = ov_start <= ov_end
        if not np.any(valid):
            repeat_pct[i] = 0.0
            repeat_types[i] = np.nan
            continue

        covered_bp = merged_covered_bp(ov_start[valid], ov_end[valid])
        intron_len = end - start + 1
        repeat_pct[i] = 100.0 * covered_bp / intron_len

        types = pd.unique(rep_type[:right][candidates][valid])
        repeat_types[i] = ", ".join(sorted(types)) if len(types) else np.nan

    out = introns_source.copy()
    out["Repeat_Percentage"] = repeat_pct
    out["Repeat_Types"] = repeat_types
    return out


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    dir_data = project_root / "resubmission" / "data"

    path_vastdb = dir_data / "external" / "EVENT_INFO-hg38.tab"
    path_repeat = dir_data / "external" / "hg38.fa.out.gz"
    path_out = dir_data / "all_vastdb_introns_repeatmasker_scores.csv"

    introns_source = load_vastdb_introns(path_vastdb)
    rep_by_chr = load_repeatmasker(path_repeat)
    introns_scored = compute_repeat_scores(introns_source, rep_by_chr)

    path_out.parent.mkdir(parents=True, exist_ok=True)
    introns_scored.to_csv(path_out, index=False)
    print(f"Wrote {len(introns_scored):,} rows to {path_out}")


if __name__ == "__main__":
    main()