#!/usr/bin/env python3

import pandas as pd
import utilities
import argparse
from Bio import SeqIO


def parse_argument():

    parser = argparse.ArgumentParser(prog = 'create_report.py')
    parser.add_argument('-s', '--sample', metavar = '', required = True, help = 'Specify sample name')
    parser.add_argument('-c', '--count_table', metavar = '', required = True, help = 'Specify count table file')
    parser.add_argument('-p', '--primers', metavar = '', required = True, help = 'Specify oligos/(primer) file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    parser.add_argument('-q', '--primer_stats', metavar = '', required = True, help = 'Specify primer_stats file')
    parser.add_argument('-f', '--fasta', metavar = '', required = True, help = 'Specify fasta file')
    parser.add_argument('-l', '--read_length', metavar = '', required = True, help = 'Specify output read_length file')
    return parser.parse_args()


def generate_primer_stats(sample, count_file, primer_stats):
    '''
    this method calculates primer stats for the given sample 
    and generate a report (tsv file) 
    #                            primer1 primer2 primer3
    #sample(abundance count)     10      20      30

    Parameters
    ----------
    sample: sample name
    count_file: original count_table file
    primer_stats: the output file

    Returns: DataFrame
    ----------
    '''   
    raw_df = pd.read_csv(count_file, sep='\t')
    raw_df.drop(['seq'], axis=1, inplace=True)
    
    #convert to the format of:
    #                            primer1 primer2 primer3
    #sample(abundance count)     10      20      30
    report_df = raw_df.sum().to_frame(name=sample).T

    # to update the column name
    # from 2014K_0979.OG0000348primerGroup3 to OG0000348primerGroup3
    def rename_column(col):
        return col.split('.', 1)[1] if '.' in col else col

    # Rename columns using a function directly
    report_df = report_df.rename(columns=rename_column)

    report_df.to_csv(f'{primer_stats}', sep='\t')
    
    return report_df


def generate_read_length(sample, fasta_file, output):
    '''
    this method calculates read length stats for the given sample 
    and generate a report (tsv file) 
    #                            num_seqs   min_len avg_len max_len
    #sample                      10         150 175 200

    Parameters
    ----------
    sample: sample name
    fasta_file: the fasta file name
    output: the output read_length file name

    Returns: DataFrame
    ----------
    '''   
    # Read the sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Calculate the statistics
    num_seqs = len(sequences)
    lengths = [len(seq.seq) for seq in sequences]
    min_length = min(lengths)
    avg_length = round(sum(lengths) / num_seqs, 1)
    max_length = max(lengths)

    stats = {
        "num_seqs": [num_seqs],
        "avg_len": [avg_length],
        "min_len": [min_length],
        "max_len": [max_length]
    }
    
    # Convert to DataFrame
    df = pd.DataFrame(stats, index=[sample])
    df.to_csv(f'{output}', sep='\t')
    
    return df


def report(sample, count_file, oligos_file, report_file):
    '''
    this method calculates metrics like: Mean read depth, # of failed primer pairs and 
    generate a report (csv file) for state public health lab
    for one single sample

    Parameters
    ----------
    sample: sample name
    count_file: original count_table file
    oligos_file: oligo file (which contains the primer information)
    report_file: the output file name

    Returns: DataFrame
    ----------
    '''   
    raw_df = pd.read_csv(count_file, sep='\t')
    raw_df.drop(['seq'], axis=1, inplace=True)
    
    #convert to the format of:
    #                            primer1 primer2 primer3
    #sample(abundance count)     10      20      30
    report_df = raw_df.sum().to_frame(name=sample).T
    
    oligo_primers = utilities.Primers(oligos_file)
    total_primer_count = len(oligo_primers.pnames)
    
    #this is failed primer pairs (common denominator) for all the samples
    all_failed_pp_count = total_primer_count - len(raw_df.columns)
    
    #including those empty cell while calculating the mean
    report_df['mean'] = report_df.fillna(0).mean(axis=1).apply(lambda x:(round(x,1)))
    report_df['failed_pp'] = report_df.isna().sum(axis=1).apply(lambda x: x+all_failed_pp_count)
    report_df['perc_successful_pp'] = round(1 - report_df['failed_pp']/total_primer_count, 3)
    # report_df['failed_pp'] = report_df['failed_pp'].apply(lambda x: f" {x} / {total_primer_count}")
    report_df['good_pp'] = report_df['failed_pp'].apply(lambda x: f" {total_primer_count-x} / {total_primer_count}")
    
    report_df = report_df[['mean','perc_successful_pp','good_pp']]
    #add 'Mean read depth across entire run' as a separate line at the end
    ### don't need this for this single sample report
    # report_df.loc['Mean read depth across entire run'] = [int(round(report_df['mean'].mean(),0)),None,None]
    
    report_columns = [f'Mean read depth',
                      f'% of successful primer-pairs\n(has at least 2 amplicons in the sample)',
                      f'# of primer pairs with more than 2 amplicons mapping\nover total primer-pairs']
    report_df.columns = report_columns
    
    report_df.to_csv(f'{report_file}')
    
    return report_df
    
    
if __name__ == "__main__":

    args = parse_argument()
    report(args.sample, args.count_table, args.primers, args.output)
    generate_primer_stats(args.sample, args.count_table, args.primer_stats)
    generate_read_length(args.sample, args.fasta, args.read_length)