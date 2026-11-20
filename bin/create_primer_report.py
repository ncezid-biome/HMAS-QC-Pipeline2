#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd

# Column names based on your BLAST command:
# -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

parser = argparse.ArgumentParser()
parser.add_argument('--input', required= True, help= "Blast Query file")
parser.add_argument('--translation_table', required= True, help= "id_to_md5.tsv")
parser.add_argument('--output', required=True, help= "Output TSV")
args = parser.parse_args()

def load_translation_table(path):
    """Return a dict mapping md5_hash -> original_id."""
    md5_to_id = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            orig_id, md5 = line.split("\t")
            md5_to_id[md5] = orig_id
    return md5_to_id

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

    # translate md5-hashed subject IDs back to the original reference IDs
    md5_to_id = load_translation_table(args.translation_table)

    unmapped = ~df['sseqid'].isin(md5_to_id)
    if unmapped.any():
        missing = sorted(df.loc[unmapped, 'sseqid'].unique())
        raise SystemExit(f"ERROR: {len(missing)} sseqid value(s) not found in "
                         f"{args.translation_table}, for example: {missing[:3]}")

    df['sseqid'] = df['sseqid'].map(md5_to_id)

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
        (df_wide['Paratyphi_ABC_18mer_nucl_sam_60_FR-Sal_RN5_01-0170_HGAP_assembled-ampl1'] >= 90) & (df_wide['Paratyphi_ABC_18mer_nucl_sam_63_FR-2009K-1733-M347-15-032_assembled-ampl1'] < 90) | (df_wide['typhi_nucl_18mer_mapping_sam_82_FR-2009K-1733-M347-15-032_assembled-ampl1'] < 90),
        (df_wide['Paratyphi_ABC_18mer_nucl_sam_60_FR-Sal_RN5_01-0170_HGAP_assembled-ampl1'] >= 90) & (df_wide['Paratyphi_ABC_18mer_nucl_sam_63_FR-2009K-1733-M347-15-032_assembled-ampl1']  >= 90) & (df_wide['typhi_nucl_18mer_mapping_sam_82_FR-2009K-1733-M347-15-032_assembled-ampl1'] >= 90),
        (df_wide['Paratyphi_ABC_18mer_nucl_sam_60_FR-Sal_RN5_01-0170_HGAP_assembled-ampl1'] < 90) & (df_wide['Paratyphi_ABC_18mer_nucl_sam_63_FR-2009K-1733-M347-15-032_assembled-ampl1'] >= 90) & (df_wide['typhi_nucl_18mer_mapping_sam_82_FR-2009K-1733-M347-15-032_assembled-ampl1'] >= 90)
    ]
    call = ["salmonella species unknown", "S. enterica", "retest"]
    df_wide['call'] = np.select(conditions, call, default = 'negative')
    df_wide.to_csv(args.output, sep = '\t', index = False)
    # print(df_wide)

if __name__ == '__main__':
    main()
