from pathlib import Path
import importlib

import numpy as np
import pandas as pd


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


def compute_phylop_scores(introns_source: pd.DataFrame, path_cons: Path) -> pd.DataFrame:
    if not path_cons.exists():
        raise FileNotFoundError(f"phyloP bigWig file not found at {path_cons}")

    try:
        pybigwig_module = importlib.import_module("pyBigWig")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Missing required dependency 'pyBigWig'. Install it in your environment, e.g. `pip install pyBigWig`."
        ) from exc

    try:
        tqdm_fn = importlib.import_module("tqdm").tqdm
    except ModuleNotFoundError:
        def tqdm_fn(iterable, total=None):
            return iterable

    records = []
    with pybigwig_module.open(str(path_cons)) as bw_phylop:
        for event, chrom, start, end in tqdm_fn(
            introns_source[["EVENT", "chrom", "start", "end"]].itertuples(index=False, name=None),
            total=introns_source.shape[0],
        ):
            try:
                values = bw_phylop.values(chrom, int(start), int(end))
            except RuntimeError:
                values = []

            arr = np.asarray(values, dtype=float) if values is not None else np.array([], dtype=float)
            if arr.size == 0:
                phylop_mean = np.nan
                has_nan = True
            else:
                phylop_mean = float(np.nanmean(arr))
                has_nan = bool(np.isnan(arr).any())

            records.append({"idx": event, "phylop_mean": phylop_mean, "has_nan": has_nan})

    return pd.DataFrame(records)


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    dir_data = project_root / "resubmission" / "data"

    path_vastdb = dir_data / "external" / "EVENT_INFO-hg38.tab"
    path_cons = dir_data / "external" / "hg38.phyloP100way.bw"

    path_out_all = dir_data / "all_vastdb_introns_phylop_scores.csv"
    path_out_clean = dir_data / "all_introns_phylop_scores.csv"

    introns_source = load_vastdb_introns(path_vastdb)
    phylop_df = compute_phylop_scores(introns_source, path_cons)

    path_out_all.parent.mkdir(parents=True, exist_ok=True)
    phylop_df.to_csv(path_out_all, index=False)
    phylop_df.drop(columns=["has_nan"]).to_csv(path_out_clean, index=False)

    print(f"Wrote {len(phylop_df):,} rows to {path_out_all}")
    print(f"Wrote {len(phylop_df):,} rows to {path_out_clean}")


if __name__ == "__main__":
    main()