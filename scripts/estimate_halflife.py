#!/usr/bin/env python3
"""
Half-Life Estimation from Transcription-Pause RNA-seq

Estimates intron and gene-level half-lives from time-series RNA-seq data.
Handles NaN values, distinguishes stable transcripts from decay, and applies
quality filtering.

Usage:
    python estimate_halflife.py --input <path_to_intron_hl.tab> --output <output_dir>
"""

import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import linregress


def calculate_empirical_stable_range(data_df, compartment=None, method='mad', k=1.0, sample_size=5000):
    """
    Calculate stable_range boundaries empirically from the data.
    
    Parameters:
    -----------
    data_df : pd.DataFrame
        dataFrame with columns: 'id', 'quant', 'Time'
    compartment : str, optional
        if provided, calculate range only for this compartment ('Cyt' or 'Nuc')
    method : str
        method to use: 'mad' (Median Absolute Deviation)
    k : float
        scaling factor for MAD
    sample_size : int
        number of introns to sample for faster computation
        
    Returns:
    --------
    tuple: (lower_bound, upper_bound)
    """
    # filter by compartment if specified
    if compartment is not None:
        data_subset = data_df[data_df['id'].str.endswith(f'_{compartment}')].copy()
    else:
        data_subset = data_df.copy()
    
    if method == 'mad':
        sample_ids = data_subset['id'].unique()[:sample_size]
        data_sample = data_subset[data_subset['id'].isin(sample_ids)].copy()
        data_sample = data_sample[~data_sample['quant'].isna()]
        
        valid_counts = data_sample.groupby('id').size()
        valid_ids = valid_counts[valid_counts >= 4].index
        data_sample = data_sample[data_sample['id'].isin(valid_ids)]
        
        medians = data_sample.groupby('id')['quant'].transform('median')
        data_sample['abs_dev'] = np.abs(data_sample['quant'] - medians)
        mad_per_id = data_sample.groupby('id')['abs_dev'].median()
        median_mad = mad_per_id.median()
        
        lower = 1.0 - k * median_mad
        upper = 1.0 + k * median_mad
        return (max(0.0, lower), upper)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'mad'")


def estimate_hl_final(quantities, timepoints, epsilon=1e-10, min_valid_points=3,
                      stable_range=(0.7, 1.3), stable_hl_hours=12):
    """
    Robust half-life estimation that distinguishes:
    1. stable transcripts (very long hl) → set to stable_hl_hours
    2. technical noise (positive slope within stable range) → set to stable_hl_hours
    3. true increases (positive slope outside range) → negative hl (filtered later)
    4. regular decay → estimate hl normally
    
    returns: pd.Series with 'hl', 'hl_error', 'n_valid', and 'is_stable' flag
    """
    quantities = np.array(quantities)
    timepoints = np.array(timepoints)
    
    # remove NaN values
    valid_mask = ~np.isnan(quantities)
    quantities_clean = quantities[valid_mask]
    timepoints_clean = timepoints[valid_mask]
    
    n_valid = len(quantities_clean)
    if n_valid < min_valid_points:
        return pd.Series({'hl': np.nan, 'hl_error': np.nan, 'n_valid': n_valid, 'is_stable': False})
    
    # replace zeros with epsilon
    quantities_clean = quantities_clean.copy()
    quantities_clean[quantities_clean == 0] = epsilon
    log_quants = np.log(quantities_clean)
    
    try:
        slope, intercept, r_value, p_value, std_err = linregress(timepoints_clean, log_quants)
        k = -slope
        
        # zero-slope case → no decay, likely stable
        if np.abs(k) < 1e-10:
            return pd.Series({
                'hl': stable_hl_hours * 60,
                'hl_error': 0,
                'n_valid': n_valid,
                'is_stable': True
            })
        
        # positive slope (k < 0, would give negative hl)
        if k < 0:
            min_val, max_val = stable_range
            original_clean = quantities[valid_mask]
            if np.all((original_clean >= min_val) & (original_clean <= max_val)):
                # technical noise - label as stable
                return pd.Series({
                    'hl': stable_hl_hours * 60,
                    'hl_error': 0,
                    'n_valid': n_valid,
                    'is_stable': True
                })
            else:
                # true increase → return negative hl (filtered later)
                hl = np.log(2) / k
                hle = (np.log(2) / k**2) * std_err
                return pd.Series({
                    'hl': hl,
                    'hl_error': hle,
                    'n_valid': n_valid,
                    'is_stable': False
                })
        
        # regular decay (positive k)
        hl = np.log(2) / k
        hle = (np.log(2) / k**2) * std_err
        
        return pd.Series({
            'hl': hl,
            'hl_error': hle,
            'n_valid': n_valid,
            'is_stable': False
        })
    
    except Exception as e:
        return pd.Series({'hl': np.nan, 'hl_error': np.nan, 'n_valid': n_valid, 'is_stable': False})


def append_hl_filter_enhanced(hl_tab, cr_filter='TP_post0_avg + 0.05 < 1', 
                               hle_quantile_filter=.75, min_valid_points=6):
    """enhanced filtering with additional quality checks."""
    hl_error_filt_gl = (
        hl_tab
            .assign(dummy = abs(hl_tab['hl_error']) / abs(hl_tab['hl']))
            .groupby('Compartment')['dummy']
            .quantile(hle_quantile_filter)
    )
    hl_tab['dummy'] = abs(hl_tab['hl_error']) / abs(hl_tab['hl'])
    hl_tab = hl_tab.merge(hl_error_filt_gl, on='Compartment', suffixes=('', '_quantile'))
    hl_tab['TP_post0_avg'] = hl_tab[['T1', 'T2', 'T3']].mean(axis=1)
    hl_tab = hl_tab.assign(
        cr_pass = hl_tab.eval(cr_filter),
        hl_pass = hl_tab.eval('hl > 0'),
        hle_pass = hl_tab.eval('dummy <= dummy_quantile'),
        valid_points_pass = hl_tab['n_valid'] >= min_valid_points,
        finite_pass = np.isfinite(hl_tab['hl'])
    )
    hl_tab.drop(columns=['TP_post0_avg', 'dummy', 'dummy_quantile'], inplace=True)
    return hl_tab


def process_intronwise_halflife(tab_hl, output_dir):
    """process and estimate intron-level half-lives."""
    print("\n=== INTRONWISE Half-Life Estimation ===")
    
    # prepare intronwise data
    tmp = tab_hl[['EVENT']].join(
        tab_hl.filter(regex=r'_1\.x$')
    ).melt(
        id_vars='EVENT',
        var_name='dummy',
        value_name='quant'
    )
    tmp[['Compartment', 'Replicate', 'Timepoint']] = tmp['dummy'].str.split('_', expand=True)[[0, 1, 2]]
    tmp.drop(columns=['dummy'], inplace=True)
    
    # extract T0 values
    tmp_t0 = tmp[tmp['Timepoint'] == 'UT'].pivot(
        index=['EVENT', 'Compartment'],
        columns='Replicate',
        values='quant'
    ).reset_index()
    tmp_t0['NA_count'] = tmp_t0[['A', 'B']].isna().sum(axis=1)
    tmp_t0['T0'] = tmp_t0[['A', 'B']].mean(axis=1)
    tmp_t0 = tmp_t0[(tmp_t0['NA_count'] == 0) & (tmp_t0['T0'] > 0)]
    
    # pivot to post-T0 timepoints and normalize
    tmp_post = tmp[tmp['Timepoint'] != 'UT']
    tmp_post = tmp_post.pivot(
        index=['EVENT', 'Compartment', 'Replicate'],
        columns='Timepoint',
        values='quant'
    ).reset_index().rename(columns={'30': 'T1', '2h': 'T2', '4h': 'T3'}).merge(
        tmp_t0[['EVENT', 'Compartment', 'T0']],
        how='right',
        on=['EVENT', 'Compartment']
    )
    
    tmp_post_norm = tmp_post[['EVENT', 'Compartment', 'Replicate', 'T0', 'T1', 'T2', 'T3']].copy()
    tmp_post_norm['T1'] = tmp_post_norm['T1'] / tmp_post_norm['T0']
    tmp_post_norm['T2'] = tmp_post_norm['T2'] / tmp_post_norm['T0']
    tmp_post_norm['T3'] = tmp_post_norm['T3'] / tmp_post_norm['T0']
    tmp_post_norm['T0'] = 1
    
    # prepare for fitting
    tmp_fit = tmp_post_norm.assign(
        id = lambda df_: df_['EVENT'] + '_' + df_['Compartment'],
    ).drop(columns=['Replicate', 'EVENT', 'Compartment']).melt(
        id_vars='id',
        value_vars=['T0', 'T1', 'T2', 'T3'],
        var_name='Timepoint',
        value_name='quant'
    ).assign(
        Time = lambda df_: df_['Timepoint'].map({'T0': 0, 'T1': 30, 'T2': 120, 'T3': 280})
    )
    
    print(f"Prepared data for {tmp_fit['id'].nunique()} introns")
    
    # calculate empirical stable ranges per compartment
    stable_range_cyt = calculate_empirical_stable_range(
        tmp_fit, compartment='Cyt', method='mad', k=1.5, sample_size=5000
    )
    stable_range_nuc = calculate_empirical_stable_range(
        tmp_fit, compartment='Nuc', method='mad', k=1.5, sample_size=5000
    )
    
    print(f"Empirical stable ranges (MAD k=1.5):")
    print(f"  Cytoplasmic: ({stable_range_cyt[0]:.3f}, {stable_range_cyt[1]:.3f})")
    print(f"  Nuclear:     ({stable_range_nuc[0]:.3f}, {stable_range_nuc[1]:.3f})")
    
    # estimate half-lives per compartment
    results = []
    for compartment, stable_range in [('Cyt', stable_range_cyt), ('Nuc', stable_range_nuc)]:
        compartment_data = tmp_fit[tmp_fit['id'].str.endswith(f'_{compartment}')].copy()
        
        hl_df_comp = compartment_data.groupby('id')[['id', 'Timepoint', 'quant', 'Time']].apply(
            lambda group: estimate_hl_final(
                group['quant'], group['Time'], 
                min_valid_points=4,
                stable_range=stable_range,
                stable_hl_hours=12
            )
        ).reset_index()
        
        results.append(hl_df_comp)
        print(f"{compartment}: Processed {len(hl_df_comp)} introns")
    
    hl_df = pd.concat(results, ignore_index=True)
    hl_df[['EVENT', 'Compartment']] = hl_df['id'].str.split('_', expand=True)[[0, 1]]
    hl_df = hl_df[['EVENT', 'Compartment', 'hl', 'hl_error', 'n_valid', 'is_stable']]
    
    # convert to hours
    hl_df['hl'] = round(hl_df['hl'] / 60, 3)
    hl_df['hl_error'] = round(hl_df['hl_error'] / 60, 3)
    
    print(f"\nTotal introns with half-life estimates: {len(hl_df)}")
    print(f"Stable transcripts identified: {hl_df['is_stable'].sum()} ({100*hl_df['is_stable'].sum()/len(hl_df):.1f}%)")
    
    # apply quality filters
    avg_vals = tmp_post_norm.groupby(['EVENT', 'Compartment'])[['T1', 'T2', 'T3']].mean().reset_index()
    hl_df = avg_vals.merge(hl_df, how='right', on=['EVENT', 'Compartment'])
    hl_df = append_hl_filter_enhanced(hl_df, cr_filter=f'TP_post0_avg < {stable_range_nuc[1]}', min_valid_points=4)
    
    hl_df_filt = hl_df[
        (hl_df['cr_pass']) & (hl_df['hl_pass']) & 
        (hl_df['hle_pass']) & (hl_df['valid_points_pass']) & 
        (hl_df['finite_pass']) & (hl_df['T1'] > 0)
    ]
    
    # gate half-life to 12 hours
    hl_df_filt = hl_df_filt.assign(hl_gated = lambda df_: df_['hl'].clip(upper=12))
    
    print(f"\n=== Filtering Summary ===")
    print(f"Total introns: {len(hl_df)}")
    print(f"High-quality filtered: {len(hl_df_filt)} ({100*len(hl_df_filt)/len(hl_df):.1f}%)")
    print(f"\nBy compartment:")
    print(hl_df_filt.groupby('Compartment').size())
    
    # save results
    hl_df.to_csv(os.path.join(output_dir, 'hl_intronwise.csv'), header=True, index=False)
    hl_df_filt.to_csv(os.path.join(output_dir, 'hl_intronwise_filtered.csv'), header=True, index=False)
    
    print(f"\nSaved results:")
    print(f"  - {os.path.join(output_dir, 'hl_intronwise.csv')} ({len(hl_df)} introns)")
    print(f"  - {os.path.join(output_dir, 'hl_intronwise_filtered.csv')} ({len(hl_df_filt)} introns)")
    
    return hl_df, hl_df_filt


def process_genewise_halflife(tab_hl, output_dir):
    """process and estimate gene-level half-lives."""
    print("\n=== GENEWISE Half-Life Estimation ===")
    
    # prepare genewise data
    tmp = tab_hl[['GENE']].drop_duplicates(subset=['GENE']).join(
        tab_hl.filter(regex=r'_1\.y$')
    ).melt(
        id_vars='GENE',
        var_name='dummy',
        value_name='quant'
    )
    tmp[['Compartment', 'Replicate', 'Timepoint']] = tmp['dummy'].str.split('_', expand=True)[[0, 1, 2]]
    tmp.drop(columns=['dummy'], inplace=True)
    
    # get T0 values
    tmp_t0 = tmp[tmp['Timepoint'] == 'UT'].pivot(
        index=['GENE', 'Compartment'],
        columns='Replicate',
        values='quant'
    ).reset_index()
    tmp_t0['NA_count'] = tmp_t0[['A', 'B']].isna().sum(axis=1)
    tmp_t0['T0'] = tmp_t0[['A', 'B']].mean(axis=1)
    tmp_t0 = tmp_t0[(tmp_t0['NA_count'] == 0) & (tmp_t0['T0'] > 0)]
    
    # get post-pause timepoints and normalize
    tmp_post = tmp[tmp['Timepoint'] != 'UT']
    tmp_post = tmp_post.pivot(
        index=['GENE', 'Compartment', 'Replicate'],
        columns='Timepoint',
        values='quant'
    ).reset_index().rename(columns={'30': 'T1', '2h': 'T2', '4h': 'T3'}).merge(
        tmp_t0[['GENE', 'Compartment', 'T0']],
        how='right',
        on=['GENE', 'Compartment']
    )
    
    tmp_post_norm = tmp_post[['GENE', 'Compartment', 'Replicate', 'T0', 'T1', 'T2', 'T3']].copy()
    tmp_post_norm['T1'] = tmp_post_norm['T1'] / tmp_post_norm['T0']
    tmp_post_norm['T2'] = tmp_post_norm['T2'] / tmp_post_norm['T0']
    tmp_post_norm['T3'] = tmp_post_norm['T3'] / tmp_post_norm['T0']
    tmp_post_norm['T0'] = 1
    
    # prepare for regression
    tmp_fit_gene = tmp_post_norm.assign(
        id = lambda df_: df_['GENE'] + '_' + df_['Compartment'],
    ).drop(columns=['Replicate', 'GENE', 'Compartment']).melt(
        id_vars='id',
        value_vars=['T0', 'T1', 'T2', 'T3'],
        var_name='Timepoint',
        value_name='quant'
    ).assign(
        Time = lambda df_: df_['Timepoint'].map({'T0': 0, 'T1': 30, 'T2': 120, 'T3': 280})
    )
    
    print(f"Prepared data for {tmp_fit_gene['id'].nunique()} genes")
    
    # calculate empirical ranges for genes
    stable_range_gene_cyt = calculate_empirical_stable_range(
        tmp_fit_gene, compartment='Cyt', method='mad', k=2.0, sample_size=5000
    )
    stable_range_gene_nuc = calculate_empirical_stable_range(
        tmp_fit_gene, compartment='Nuc', method='mad', k=2.0, sample_size=5000
    )
    
    print(f"Gene-level empirical stable ranges:")
    print(f"  Cytoplasmic: ({stable_range_gene_cyt[0]:.3f}, {stable_range_gene_cyt[1]:.3f})")
    print(f"  Nuclear:     ({stable_range_gene_nuc[0]:.3f}, {stable_range_gene_nuc[1]:.3f})")
    
    # estimate half-lives
    results_gene = []
    for compartment, stable_range in [('Cyt', stable_range_gene_cyt), ('Nuc', stable_range_gene_nuc)]:
        compartment_data = tmp_fit_gene[tmp_fit_gene['id'].str.endswith(f'_{compartment}')].copy()
        
        hl_df_comp = compartment_data.groupby('id')[['id', 'Timepoint', 'quant', 'Time']].apply(
            lambda group: estimate_hl_final(
                group['quant'], group['Time'], 
                min_valid_points=4, stable_range=stable_range, stable_hl_hours=12
            )
        ).reset_index()
        
        results_gene.append(hl_df_comp)
    
    hl_df_gene = pd.concat(results_gene, ignore_index=True)
    hl_df_gene[['GENE', 'Compartment']] = hl_df_gene['id'].str.split('_', expand=True)[[0, 1]]
    hl_df_gene = hl_df_gene[['GENE', 'Compartment', 'hl', 'hl_error', 'n_valid', 'is_stable']]
    hl_df_gene['hl'] = round(hl_df_gene['hl'] / 60, 3)
    hl_df_gene['hl_error'] = round(hl_df_gene['hl_error'] / 60, 3)
    
    # apply filters
    avg_vals_gene = tmp_post_norm.groupby(['GENE', 'Compartment'])[['T1', 'T2', 'T3']].mean().reset_index()
    hl_df_gene = avg_vals_gene.merge(hl_df_gene, how='right', on=['GENE', 'Compartment'])
    hl_df_gene = append_hl_filter_enhanced(hl_df_gene, cr_filter=f'TP_post0_avg < {stable_range_gene_nuc[1]}', min_valid_points=4)
    
    hl_df_gene_filt = hl_df_gene[
        (hl_df_gene['cr_pass']) & (hl_df_gene['hl_pass']) & 
        (hl_df_gene['hle_pass']) & (hl_df_gene['valid_points_pass']) & 
        (hl_df_gene['finite_pass']) & (hl_df_gene['T1'] > 0)
    ]
    
    # gate half-life to 12 hours
    hl_df_gene_filt = hl_df_gene_filt.assign(hl_gated = lambda df_: df_['hl'].clip(upper=12))
    
    # save
    hl_df_gene.to_csv(os.path.join(output_dir, 'hl_genewise.csv'), header=True, index=False)
    hl_df_gene_filt.to_csv(os.path.join(output_dir, 'hl_genewise_filtered.csv'), header=True, index=False)
    
    print(f"\nGene-level results:")
    print(f"  Total: {len(hl_df_gene)}")
    print(f"  Filtered: {len(hl_df_gene_filt)}")
    print(f"  Saved to {os.path.join(output_dir, 'hl_genewise*.csv')}")
    
    return hl_df_gene, hl_df_gene_filt


def main():
    parser = argparse.ArgumentParser(
        description='Estimate half-lives from transcription-pause RNA-seq data'
    )
    parser.add_argument(
        '--input', '-i', 
        required=True,
        help='Path to intron_hl.tab file'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    # create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print(f"Loading data from {args.input}...")
    tab_hl = pd.read_csv(args.input, sep='\t')
    
    # append strand information to intron coordinates
    tab_hl['COORD'] = tab_hl['COORD'] + ':' + tab_hl['FullCO'].str.rsplit(':', n=3).str[-1]
    tab_hl = tab_hl.drop_duplicates(subset=['GENE', 'EVENT'], keep='first')
    
    print(f"Loaded {len(tab_hl)} unique events")
    print(f"Output directory: {args.output}")
    
    # process intron-level half-lives
    hl_intron, hl_intron_filt = process_intronwise_halflife(tab_hl, args.output)
    
    # process gene-level half-lives
    hl_gene, hl_gene_filt = process_genewise_halflife(tab_hl, args.output)
    
    print("\n=== Complete ===")
    print(f"All results saved to {args.output}")


if __name__ == '__main__':
    main()
