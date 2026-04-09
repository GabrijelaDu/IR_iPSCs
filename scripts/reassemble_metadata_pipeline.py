#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Mapping

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "resubmission" / "data"

DEFAULT_INPUTS = {
    "base_metadata": DATA_DIR / "combined_dataTPMno100.csv",
    "introns_source": DATA_DIR / "all_vastdb_introns_source.csv",
    "repeat_masker": DATA_DIR / "all_vastdb_introns_repeatmasker_scores.csv",
    "phylop": DATA_DIR / "all_vastdb_introns_phylop_scores.csv",
    "ssdrip_peaks": DATA_DIR / "all_vastdb_introns_ssDRIP_peaks.csv",
    "ssdrip_reads": DATA_DIR / "all_vastdb_introns_ssDRIP_reads.csv",
    "hl_gene": DATA_DIR / "hl" / "hl_genewise_filtered.csv",
    "hl_intron": DATA_DIR / "hl" / "hl_intronwise_filtered.csv",
    "speckles": DATA_DIR / "external" / "speckles_enrichment.csv",
    "cpg": DATA_DIR / "all_vastdb_introns_source_cpg_methylation.csv",
}

DEFAULT_OUTPUTS = {
    "metadata_selected": DATA_DIR / "metadata_selected.csv",
}

STABILITY_LABEL_MAP = {
    "no IR": "noIR",
    "unstable2h": "Unstable_2h",
    "stable 4h": "Stable_4h",
    "unstable30min": "Unstable_30min",
    "unstable 4h": "Unstable_4h",
}

SELECTED_COLUMNS = [
    "GENE",
    "EVENT",
    "LENGTH",
    "Total_Introns",
    "Intron_Position",
    "Repeat_Percentage",
    "Repeat_Types",
    "GC_Content",
    "phylop_mean",
    "PIR_Nuc_baseline",
    "TPM_Nuc_baseline",
    "TPM_Cyt_baseline",
    "Nuclear_Enrichment",
    "hl_gated_gwise",
    "hl_gated_intwise",
    "hl_gated_intwise_scaled",
    "hl_gated_intwise_scaled_logged",
    "stability",
    "CpG_total",
    "CpG_methylated_perc",
    "Speckle_Enrichment",
    "ssDRIP_peaks",
    "ssDRIP_qValue",
    "ssDRIP_cov",
]


def resolve_path(path_like: str | Path) -> Path:
    path = Path(path_like)
    return path if path.is_absolute() else REPO_ROOT / path


def build_config(
    *,
    inputs: Mapping[str, str | Path] | None = None,
    outputs: Mapping[str, str | Path] | None = None,
) -> Dict[str, Path]:
    config: Dict[str, Path] = {}
    for key, default in DEFAULT_INPUTS.items():
        config[key] = resolve_path(inputs[key]) if inputs and key in inputs else default
    for key, default in DEFAULT_OUTPUTS.items():
        config[key] = resolve_path(outputs[key]) if outputs and key in outputs else default
    return config


def validate_config(config: Mapping[str, Path], verbose: bool = True) -> List[str]:
    missing: List[str] = []
    for key, path in config.items():
        if key in DEFAULT_OUTPUTS:
            continue
        ok = path.exists()
        if verbose:
            print(f"[{'OK' if ok else 'MISSING'}] {key}: {path}")
        if not ok:
            missing.append(key)
    return missing


def normalize_stability_values(tab: pd.DataFrame) -> pd.DataFrame:
    out = tab.copy()
    if "stability" in out.columns:
        out["stability"] = out["stability"].replace(STABILITY_LABEL_MAP)
    return out


def normalize_base_metadata(tab_meta_gabi: pd.DataFrame) -> pd.DataFrame:
    tab = tab_meta_gabi.copy()
    tab.columns = tab.columns.str.replace('.x', '_PIR', regex=False).str.replace('.y', '_TPM', regex=False)
    drop_cols = ['X', 'non_na_count', 'COORD', 'FullCO', 'COMPLEX', 'chromosome', 'start', 'end', 'Sequence']
    tab = tab.drop(columns=[col for col in drop_cols if col in tab.columns])
    tab = tab.drop_duplicates(subset='EVENT')
    tab = normalize_stability_values(tab)
    tab = tab.drop(columns=tab.columns[tab.columns.str.startswith(('Nuc_A', 'Nuc_B', 'Cyt_A', 'Cyt_B'))])
    tab = tab.drop(columns=tab.columns[tab.columns.str.contains('_30_|_2h_|_4h_')])
    tab['Nuclear_Enrichment'] = np.log2(tab['Nuc_UT_TPM_mean'] / tab['Cyt_UT_TPM_mean'].replace(0, np.nan))
    return tab


def append_introns_source(tab: pd.DataFrame, path_introns_source: Path) -> pd.DataFrame:
    introns_source = pd.read_csv(path_introns_source)
    return tab.merge(
        introns_source[['GENE', 'EVENT', 'gene_type', 'Total_Introns', 'Intron_Position']],
        on=['GENE', 'EVENT'],
        how='left',
    )


def append_repeat_masker(tab: pd.DataFrame, path_repeat_masker: Path) -> pd.DataFrame:
    tab_repeat_masker = pd.read_csv(path_repeat_masker)
    return tab.merge(
        tab_repeat_masker[['GENE', 'EVENT', 'Repeat_Percentage', 'Repeat_Types']],
        on=['GENE', 'EVENT'],
        how='left',
    )


def append_phylop(tab: pd.DataFrame, path_phylop: Path) -> pd.DataFrame:
    tab_phylop = pd.read_csv(path_phylop)
    return tab.merge(
        tab_phylop[['idx', 'phylop_mean']],
        left_on='EVENT',
        right_on='idx',
        how='left',
    ).drop(columns=['idx'])


def append_ssdrip(tab: pd.DataFrame, path_peaks: Path, path_reads: Path) -> pd.DataFrame:
    drip_pks = pd.read_csv(path_peaks).rename(
        columns={
            'ssDRIP_overlap_bp_strand_matched': 'ssDRIP_peaks',
            'ssDRIP_qValue_max_strand_matched': 'ssDRIP_qValue',
        }
    )
    out = tab.merge(drip_pks[['EVENT', 'ssDRIP_peaks', 'ssDRIP_qValue']], on='EVENT', how='left')

    drip_cov = pd.read_csv(path_reads)
    total_reads_cols = drip_cov.columns[drip_cov.columns.str.startswith('total_reads')]
    col_sums = drip_cov[total_reads_cols].sum(axis=0).replace(0, np.nan)
    drip_cov[total_reads_cols] = drip_cov[total_reads_cols].div(col_sums, axis=1) * 1e6
    length_kb = (drip_cov['end'] - drip_cov['start']).replace(0, np.nan) / 1e3
    drip_cov[total_reads_cols] = drip_cov[total_reads_cols].div(length_kb, axis=0)
    avg_batch2_cols = drip_cov.columns[drip_cov.columns.str.startswith('total_reads_batch2')]
    drip_cov['ssDRIP_cov'] = np.log2(drip_cov[avg_batch2_cols].mean(axis=1) + 1)
    drip_cov = drip_cov.drop(columns=total_reads_cols)
    return out.merge(drip_cov[['EVENT', 'ssDRIP_cov']], on='EVENT', how='left')


def append_half_lives(tab: pd.DataFrame, path_hl_gene: Path, path_hl_intron: Path) -> pd.DataFrame:
    hl_gwise = pd.read_csv(path_hl_gene)
    hl_gwise_sel = hl_gwise.query('Compartment == "Nuc"')[['GENE', 'hl_gated']].rename(columns={'hl_gated': 'hl_gated_gwise'})
    hl_iwise = pd.read_csv(path_hl_intron)
    hl_iwise_sel = hl_iwise.query('Compartment == "Nuc"')[['EVENT', 'hl_gated']].rename(columns={'hl_gated': 'hl_gated_intwise'})
    out = tab.merge(hl_gwise_sel, on='GENE', how='left').merge(hl_iwise_sel, on='EVENT', how='left')
    out['hl_gated_intwise_scaled'] = out['hl_gated_intwise'] / out['hl_gated_gwise']
    out['hl_gated_intwise_scaled_logged'] = np.log10(out['hl_gated_intwise_scaled'])
    return out


def append_speckles(tab: pd.DataFrame, path_speckles: Path) -> pd.DataFrame:
    tab_speck = pd.read_csv(path_speckles).query('baseMean > 10')
    return tab.merge(
        tab_speck[['gene_name', 'log2FoldChange']].rename(columns={'gene_name': 'GENE', 'log2FoldChange': 'Speckle_Enrichment'}),
        on='GENE',
        how='left',
    )


def append_cpg(tab: pd.DataFrame, path_cpg: Path) -> pd.DataFrame:
    tab_cpg = pd.read_csv(path_cpg).rename(columns={'CpG_sites_in_region': 'CpG_total', 'CpG_methylation_pct': 'CpG_methylated_perc'})
    return tab.merge(tab_cpg[['EVENT', 'CpG_total', 'CpG_methylated_perc']], on='EVENT', how='left')


def finalize_selected_metadata(tab: pd.DataFrame) -> pd.DataFrame:
    out = tab.rename(columns={'Nuc_UT_mean': 'PIR_Nuc_baseline', 'Nuc_UT_TPM_mean': 'TPM_Nuc_baseline', 'Cyt_UT_TPM_mean': 'TPM_Cyt_baseline'})
    out['PIR_Nuc_baseline'] = pd.to_numeric(out['PIR_Nuc_baseline'].replace('#DIV/0!', np.nan), errors='coerce').astype(float)
    out = out[SELECTED_COLUMNS].copy()
    out.index = out['EVENT']
    out.index.name = 'EVENT'
    return out.drop(columns='EVENT')


def assemble_metadata(config: Mapping[str, Path]) -> pd.DataFrame:
    tab = pd.read_csv(config['base_metadata'])
    tab = normalize_base_metadata(tab)
    tab = append_introns_source(tab, config['introns_source'])
    tab = append_repeat_masker(tab, config['repeat_masker'])
    tab = append_phylop(tab, config['phylop'])
    tab = append_ssdrip(tab, config['ssdrip_peaks'], config['ssdrip_reads'])
    tab = append_half_lives(tab, config['hl_gene'], config['hl_intron'])
    tab = append_speckles(tab, config['speckles'])
    tab = append_cpg(tab, config['cpg'])
    return finalize_selected_metadata(tab)


def summarize_selected_metadata(tab_selected: pd.DataFrame) -> pd.Series:
    return pd.Series({
        'rows': int(tab_selected.shape[0]),
        'columns': int(tab_selected.shape[1]),
        'missing_hl_gene': int(tab_selected['hl_gated_gwise'].isna().sum()),
        'missing_hl_intron': int(tab_selected['hl_gated_intwise'].isna().sum()),
        'missing_speckles': int(tab_selected['Speckle_Enrichment'].isna().sum()),
        'missing_cpg': int(tab_selected['CpG_total'].isna().sum()),
        'missing_ssdrip_cov': int(tab_selected['ssDRIP_cov'].isna().sum()),
    })


def summarize_stability_labels(tab_selected: pd.DataFrame) -> pd.Series:
    return tab_selected['stability'].value_counts(dropna=False).rename('count')


def save_selected_metadata(tab_selected: pd.DataFrame, path_output: Path) -> None:
    path_output.parent.mkdir(parents=True, exist_ok=True)
    tab_selected.to_csv(path_output, index=True)


def run_pipeline(*, inputs: Mapping[str, str | Path] | None = None, outputs: Mapping[str, str | Path] | None = None, save: bool = True) -> pd.DataFrame:
    config = build_config(inputs=inputs, outputs=outputs)
    missing = validate_config(config, verbose=True)
    if missing:
        raise SystemExit(f"Missing required inputs: {', '.join(missing)}")
    tab_selected = assemble_metadata(config)
    if save:
        save_selected_metadata(tab_selected, config['metadata_selected'])
        print(f"Saved metadata_selected to {config['metadata_selected']}")
    return tab_selected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Reassemble metadata_selected.csv for the resubmission workflow.')
    parser.add_argument('--output', default=str(DEFAULT_OUTPUTS['metadata_selected']))
    parser.add_argument('--validate-only', action='store_true')
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = build_config(outputs={'metadata_selected': args.output})
    missing = validate_config(config, verbose=True)
    if missing:
        raise SystemExit(1)
    if args.validate_only:
        return
    tab_selected = assemble_metadata(config)
    save_selected_metadata(tab_selected, config['metadata_selected'])
    print(summarize_selected_metadata(tab_selected))


if __name__ == '__main__':
    main()
