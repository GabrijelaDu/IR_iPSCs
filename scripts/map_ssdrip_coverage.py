from pathlib import Path
import importlib
import re

import pandas as pd


def get_tqdm():
    try:
        return importlib.import_module("tqdm").tqdm
    except ModuleNotFoundError:
        def tqdm_fallback(iterable, total=None, desc=None):
            return iterable
        return tqdm_fallback


def load_regions(path_regions: Path) -> pd.DataFrame:
    if not path_regions.exists():
        raise FileNotFoundError(f"Region file not found: {path_regions}")

    regions = pd.read_csv(path_regions, low_memory=False)
    required = {"EVENT", "chrom", "start", "end", "strand"}
    missing = required.difference(regions.columns)
    if missing:
        raise ValueError(f"Region file missing required columns: {sorted(missing)}")

    regions = regions.copy()
    regions["EVENT"] = regions["EVENT"].astype(str)
    regions["start"] = pd.to_numeric(regions["start"], errors="coerce")
    regions["end"] = pd.to_numeric(regions["end"], errors="coerce")
    regions = regions.dropna(subset=["EVENT", "chrom", "start", "end", "strand"])
    regions["start"] = regions["start"].astype(int)
    regions["end"] = regions["end"].astype(int)
    regions["strand"] = regions["strand"].astype(str).str.strip()
    regions = regions[regions["end"] >= regions["start"]].reset_index(drop=True)
    regions = regions.drop_duplicates(subset=["EVENT"], keep="first").reset_index(drop=True)
    return regions


def parse_batch_rep(sample_dir: Path) -> tuple[str, str]:
    match = re.search(r"batch(\d+)_rep(\d+)", sample_dir.name)
    if match is None:
        raise ValueError(
            f"Could not parse batch/rep from folder name: {sample_dir.name}. "
            "Expected pattern like '...batch1_rep2'."
        )
    return f"batch{match.group(1)}", f"rep{match.group(2)}"


def find_strand_bams(sample_name: str, tmp_dir: Path) -> dict[str, Path]:
    bam_fwd = tmp_dir / f"{sample_name}_Aligned.sortedByCoord.out_duplRm_fwd.bam"
    bam_rev = tmp_dir / f"{sample_name}_Aligned.sortedByCoord.out_duplRm_rev.bam"
    if not bam_fwd.exists():
        raise FileNotFoundError(f"Missing forward-strand BAM: {bam_fwd}")
    if not bam_rev.exists():
        raise FileNotFoundError(f"Missing reverse-strand BAM: {bam_rev}")
    return {"fwd": bam_fwd, "rev": bam_rev}


def count_reads_per_intron_for_replicate(strand_bams: dict[str, Path], regions: pd.DataFrame) -> list[int]:
    try:
        pysam = importlib.import_module("pysam")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Missing required dependency 'pysam'. Install it with `pip install pysam`."
        ) from exc

    counts: list[int] = []
    tqdm_fn = get_tqdm()

    with (
        pysam.AlignmentFile(str(strand_bams["fwd"]), "rb") as bam_fwd,
        pysam.AlignmentFile(str(strand_bams["rev"]), "rb") as bam_rev,
    ):
        intron_iter = regions[["chrom", "start", "end", "strand"]].itertuples(index=False, name=None)
        for chrom, start_1based, end_1based, strand in tqdm_fn(
            intron_iter,
            total=len(regions),
            desc="Introns",
        ):
            if strand == "+":
                bam = bam_fwd
            elif strand == "-":
                bam = bam_rev
            else:
                counts.append(0)
                continue

            start_0based = max(0, int(start_1based) - 1)
            end_exclusive = int(end_1based)

            try:
                n_reads = bam.count(
                    str(chrom),
                    start_0based,
                    end_exclusive,
                    read_callback="all",
                )
                counts.append(int(n_reads))
            except ValueError:
                counts.append(0)

    return counts


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]

    path_regions = project_root / "resubmission" / "data" / "all_vastdb_introns_source.csv"
    tmp_dir = project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "tmp"
    sample_dirs = [
        project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "ssDRIP_iPSC_batch1_rep1",
        project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "ssDRIP_iPSC_batch1_rep2",
        project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "ssDRIP_iPSC_batch2_rep1",
        project_root / "ipsc" / "data_raw" / "PRJNA608890" / "ssDRIP" / "ssDRIP_iPSC_batch2_rep2",
    ]

    out_dir = project_root / "resubmission" / "data"
    path_out = out_dir / "all_vastdb_introns_ssDRIP_reads.csv"

    tqdm_fn = get_tqdm()
    regions = load_regions(path_regions)

    out_df = regions[["EVENT", "chrom", "start", "end", "strand"]].copy()
    for sample_dir in tqdm_fn(sample_dirs, total=len(sample_dirs), desc="Replicates"):
        if not sample_dir.exists():
            raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

        batch, replicate = parse_batch_rep(sample_dir)
        sample_name = sample_dir.name
        strand_bams = find_strand_bams(sample_name, tmp_dir)
        reads = count_reads_per_intron_for_replicate(strand_bams, regions)
        col_name = f"total_reads_{batch}_{replicate}"
        out_df[col_name] = reads

    out_dir.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(path_out, index=False)
    print(f"Wrote per-intron read counts to {path_out}")


if __name__ == "__main__":
    main()
