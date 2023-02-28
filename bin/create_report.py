#!/usr/bin/env python3

import pandas as pd
import utilities
import argparse


def parse_argument():

    parser = argparse.ArgumentParser(prog = 'create_report.py')
    parser.add_argument('-s', '--sample', metavar = '', required = True, help = 'Specify sample name')
    parser.add_argument('-c', '--count_table', metavar = '', required = True, help = 'Specify count table file')
    parser.add_argument('-p', '--primers', metavar = '', required = True, help = 'Specify oligos/(primer) file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    return parser.parse_args()


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
    None
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
    report_df['failed_pp'] = report_df['failed_pp'].apply(lambda x: f" {x} / {total_primer_count}")
    
    report_df = report_df[['mean','perc_successful_pp','failed_pp']]
    #add 'Mean read depth across entire run' as a separate line at the end
    ### don't need this for this single sample report
    # report_df.loc['Mean read depth across entire run'] = [int(round(report_df['mean'].mean(),0)),None,None]
    
    report_columns = [f'Mean read depth',
                      f'% of successful primer-pairs\n(has at least 2 amplicons in the sample)',
                      f'# of primer pairs with less than 2 amplicons mapping\nover total primer-pairs']
    report_df.columns = report_columns
    
    report_df.to_csv(f'{report_file}.csv')
    
    return report_df
    
    
if __name__ == "__main__":

    args = parse_argument()
    report(args.sample, args.count_table, args.primers, args.output)