#!/usr/bin/env python3
"""
sample sequences by stratified sampling and compute MFE with RNAfold.
"""

import argparse
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from ViennaRNA import RNA


def run_rnafold_wrapper(item):
    """wrapper for parallel execution."""
    return run_rnafold(item[0], item[1])


def run_rnafold(seq_id, sequence):
    """run RNAfold for a single sequence and return MFE."""
    try:
        fc = RNA.fold_compound(sequence)
        structure, mfe = fc.mfe()
        return seq_id, structure, float(mfe)
    except Exception as e:
        print(f"error processing {seq_id}: {e}")
        return seq_id, None, None


def main():
    parser = argparse.ArgumentParser(description='stratified sampling + RNAfold MFE')
    parser.add_argument('--input', required=True, help='input CSV with metadata')
    parser.add_argument('--column', required=True, help='column for stratified sampling')
    parser.add_argument('--sample-size', type=int, required=True, help='samples per group')
    parser.add_argument('--fasta', default='resubmission/data/all_vastdb_introns_50fks.fasta', help='FASTA file')
    parser.add_argument('--output', required=True, help='output CSV path')
    parser.add_argument('--index-col', default=0, type=int, help='index column in CSV (default: 0)')
    parser.add_argument('--workers', type=int, default=8, help='parallel workers')
    parser.add_argument('--seed', type=int, default=42, help='random seed')
    
    args = parser.parse_args()
    
    # load metadata
    print(f"loading metadata from {args.input}...")
    df = pd.read_csv(args.input, index_col=args.index_col)
    
    # stratified sampling
    print(f"sampling {args.sample_size} per group from column '{args.column}'...")
    sampled = df.groupby(args.column, group_keys=False).apply(
        lambda x: x.sample(min(len(x), args.sample_size), random_state=args.seed),
        include_groups=False
    )
    sampled_ids = set(sampled.index.astype(str))
    print(f"sampled {len(sampled_ids)} sequences from {df[args.column].nunique()} groups")
    
    # load sequences
    print(f"loading sequences from {args.fasta}...")
    seqs = {}
    for record in SeqIO.parse(args.fasta, 'fasta'):
        seq_id = record.id  #.split('_')[0]  # extract EVENT id
        if seq_id in sampled_ids:
            seqs[seq_id] = str(record.seq)
    
    print(f"matched {len(seqs)} sequences")
    
    # compute MFE in parallel
    print(f"computing MFE with RNAfold using {args.workers} workers...")
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(run_rnafold_wrapper, seqs.items()))
    
    # compile results
    results_df = pd.DataFrame(results, columns=['EVENT', 'structure', 'MFE'])
    results_df = results_df[results_df['MFE'].notna()]  # filter failed
    
    # merge with metadata
    results_df = results_df.set_index('EVENT').join(df[[args.column]], how='left')
    
    # save
    results_df.to_csv(args.output)
    print(f"saved {len(results_df)} results to {args.output}")
    print(f"failed: {len(seqs) - len(results_df)}")


if __name__ == '__main__':
    main()
