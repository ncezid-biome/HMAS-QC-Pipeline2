#!/usr/bin/env python3

import pandas as pd
import datetime
import argparse
import os, sys

'''
This script will generat a mothur_equivalent full format count_table file.
The only difference is that it does have the 2nd column (total) of mothur's count table, and the name of
the first columns is 'seq' instead of 'Representative_Sequence'
'''

def parse_argument():
    # note
    # the script can be run as: python3 make_count_table.py -o fina.count_table -m output.match.final.txt 
    #
    parser = argparse.ArgumentParser(prog = 'make_count_table.py')
    parser.add_argument('-o', '--output_file', metavar = '', required = True, help = 'Specify output file name')
    parser.add_argument('-m', '--match_file', metavar = '', required = True, help = 'Specify the matched file')

    return parser.parse_args()


def get_sample_primer(row):
    sep = '='
    sample = row.split(sep)[2]
    primer = row.split(sep)[1]

    return f"{sample}.{primer}"
    
def main():
    
    args = parse_argument()
    
    # file_dir = '/scicomp/home-pure/qtl7/test/hmas_test/024_demulx_data/output'
    # file_dir = args.output_dir
    # matched_file = 'output.match.final.txt'
    matched_file = args.match_file
    output_file = args.output_file
    #fasta_file = 'output.fasta'
    
    if os.path.getsize(matched_file) <= 0:
        print (f"{matched_file} is empty" \
            	", indicating the raw reads are most likely bad" \
                 ", and all generated seqs has abundance less than 10")
        sys.exit()

    # df = pd.read_csv(f'{file_dir}/{matched_file}', sep='\t', names = ['target', 'query'])
    df = pd.read_csv(f'{matched_file}', sep='\t', names = ['target', 'query'])
    total_seqs = df['query'].to_list()
    unique_seqs = list(set(df['target'].to_list()))

    df['query'] = df.groupby(['target'])['query'].transform(lambda x: ','.join(x))
    print (df.head())
    print (df.shape)

    #dictionary (key: unique representative seqs; value: list of its represented seqs)
    seq_seq_dict = dict(df.values)
    print (len(seq_seq_dict))
    # for key in list(seq_seq_dict.keys())[:2]:
    #     print (seq_seq_dict[key])

    #dictionary (key: unique representative seqs; value: list of the sample.primer pairs of its represented seqs)
    seq_sp_dict = {key:[get_sample_primer(seq) for seq in seq_seq_dict[key].split(',')] for key in seq_seq_dict}
    print (len(seq_sp_dict))
    # for key in list(seq_sp_dict.keys())[:2]:
    #     print (seq_sp_dict[key])

    # seq_sp_dict.values()
    # def get_set_union(set_list):
    #     return set.union(*set_list)

    #list of sets (of sample-primer for each unique seqs)
    sp_sets = [set(item) for item in seq_sp_dict.values()]
    print (len(sp_sets))
    for ss in sp_sets[:2]:
        print (ss)

    #this is the total sample-primer pairs in the count_table
    sps = list(set.union(*sp_sets))
    sps.sort()
    print (len(sps))


    print (f"before start...{datetime.datetime.now()}")
    #making the rows of the count_table, which is a dictionary of:
    #key: unique seqs; value: list of counts for each sample-primer pair
    occurrence = {seq: [seq_sp_dict[seq].count(sp) for sp in sps] 
                for seq in unique_seqs}

    print (len(occurrence))
    print (unique_seqs[:20])
    # print (occurrence.keys())
    for key in list(occurrence.keys())[:20]:
        print (key, occurrence[key])

    df = pd.DataFrame.from_dict(occurrence, orient='index')
    df.columns = sps




    # df = pd.read_pickle(f"{file_dir}/024_not_mothur_count_table.pkl")


    df.reset_index(inplace=True)
    df.rename(columns={'index':'seq'}, inplace=True)
    df['seq'] = df['seq'].apply(str)

    print (df.head().iloc[:5,:5])


    print (df.shape)
    print (f"after done... {datetime.datetime.now()}")

    # pd.to_pickle(df, f"{file_dir}/024_not_mothur_count_table_2.pkl") 
    pd.to_pickle(df, f"{output_file}.pkl") 

if __name__ == "__main__":
    main()