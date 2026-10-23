#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd

# Column names based on your BLAST command:
# -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

parser = argparse.ArgumentParser()
parser.add_argument('--input', required= True, help= "Blast Query file")
parser.add_argument('--output', required=True, help= "Output TSV")
args = parser.parse_args()

def main():
    col_names = [
    'qseqid',    # query sequence ID
    'sseqid',    # subject (database) sequence ID
    'pident',    # percentage of identical matches
    'length',    # alignment length
    'mismatch',  # number of mismatches
    'gapopen',   # number of gap openings
    'qstart',    # start of alignment in query
    'qend',      # end of alignment in query
    'sstart',    # start of alignment in subject
    'send',      # end of alignment in subject
    'evalue',    # expect value
    'bitscore',  # bit score
    'qlen',      # query sequence length
    'slen'       # subject sequence length
    ]

    df = pd.read_csv(
    args.input,
    sep='\t',
    header=None,
    names=col_names
    )

    split_cols = df['qseqid'].str.split('=', expand=True)
    df['read_id'] = split_cols[0]
    df['mapping_target'] = split_cols[1]
    df['sample_id'] = split_cols[2].str.replace(';size', '', regex=False)
    df['read_size'] = split_cols[3].astype(int)

    # sort df by sseqid, read_size, then pident in order before subsetting sample with highest abundance
    df_sorted = df.sort_values(by=['sseqid', 'read_size', 'pident'], ascending=[True, False, False])
    df_top = df_sorted.groupby('sseqid').first().reset_index()

    # print(df_top[['sseqid', 'read_id', 'mapping_target', 'read_size', 'pident']])
    # print(df_top.dtypes)

    df_wide = df_top.pivot(index='sample_id', columns='sseqid', values='pident')
    df_wide = df_wide.reset_index()

    # boolean conditions for specimen identification
    conditions = [
        (df_wide['Paratyphi_60_FR-Sal_RN5_01-0170'] >= 90) & (df_wide['Paratyphi_63_FR-2009K-1733'] < 90) | (df_wide['typhi_82_FR-2009K-1733'] < 90),
        (df_wide['Paratyphi_60_FR-Sal_RN5_01-0170'] >= 90) & (df_wide['Paratyphi_63_FR-2009K-1733']  >= 90) & (df_wide['typhi_82_FR-2009K-1733'] >= 90),
        (df_wide['Paratyphi_60_FR-Sal_RN5_01-0170'] < 90) & (df_wide['Paratyphi_63_FR-2009K-1733'] >= 90) & (df_wide['typhi_82_FR-2009K-1733'] >= 90)
    ]
    call = ["salmonella species unknown", "S. enterica", "retest"]
    df_wide['call'] = np.select(conditions, call, default = 'negative')
    df_wide.to_csv(args.output, sep = '\t', index = False)
    # print(df_wide)

if __name__ == '__main__':
    main()
