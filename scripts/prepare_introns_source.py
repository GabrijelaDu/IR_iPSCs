from pathlib import Path

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


def add_intron_gene_features(introns_source: pd.DataFrame) -> pd.DataFrame:
    out = introns_source.copy()
    out["Total_Introns"] = out.groupby("GENE")["EVENT"].transform("size")
    out["Intron_Position"] = pd.NA

    plus_mask = out["strand"].eq("+")
    minus_mask = out["strand"].eq("-")
    other_mask = ~(plus_mask | minus_mask)

    out.loc[plus_mask, "Intron_Position"] = (
        out.loc[plus_mask]
        .sort_values(["GENE", "start", "end"], ascending=[True, True, True])
        .groupby("GENE")
        .cumcount()
        .add(1)
    )

    out.loc[minus_mask, "Intron_Position"] = (
        out.loc[minus_mask]
        .sort_values(["GENE", "end", "start"], ascending=[True, False, False])
        .groupby("GENE")
        .cumcount()
        .add(1)
    )

    out.loc[other_mask, "Intron_Position"] = (
        out.loc[other_mask]
        .sort_values(["GENE", "start", "end"], ascending=[True, True, True])
        .groupby("GENE")
        .cumcount()
        .add(1)
    )
    out["Intron_Position"] = out["Intron_Position"].astype("Int64")
    return out


def load_gene_types(path_gtf: Path) -> pd.DataFrame:
    if not path_gtf.exists():
        raise FileNotFoundError(
            f"GENCODE annotation file not found at {path_gtf}. "
            "Please get one at https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz."
        )

    gene_annot = pd.read_csv(path_gtf, low_memory=False)
    gene_annot = (
        gene_annot[gene_annot["feature"] == "gene"]
        .dropna(subset=["gene_name"])[["gene_name", "gene_type"]]
        .drop_duplicates(subset=["gene_name"], keep="first")
        .reset_index(drop=True)
    )
    return gene_annot


def add_gene_type(introns_source: pd.DataFrame, gene_annot: pd.DataFrame) -> pd.DataFrame:
    out = introns_source.merge(
        gene_annot,
        left_on="GENE",
        right_on="gene_name",
        how="left",
    ).drop(columns=["gene_name"])
    return out


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    dir_data = project_root / "resubmission" / "data"

    path_vastdb = dir_data / "external" / "EVENT_INFO-hg38.tab"
    path_gtf = dir_data / "external" / "gencode.v44.basic.annotation.csv.gz"
    path_out = dir_data / "all_vastdb_introns_source.csv"

    introns_source = load_vastdb_introns(path_vastdb)
    introns_source = add_intron_gene_features(introns_source)
    gene_annot = load_gene_types(path_gtf)
    introns_source = add_gene_type(introns_source, gene_annot)

    path_out.parent.mkdir(parents=True, exist_ok=True)
    introns_source.to_csv(path_out, index=False)
    print(f"Wrote {len(introns_source):,} rows to {path_out}")


if __name__ == "__main__":
    main()




