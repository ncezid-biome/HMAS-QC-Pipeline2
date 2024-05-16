#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import re

'''
This script reads all report.csv (for each sample) and concantenate them into one single report.csv. It will 
also add a 'Mean read depth across entire run' row at the end.

This script requires all report.csv file names be passed in a concatenated string as a command line argument.
'''

def list_files_in_folder(folder_path):
    """
    Return a list of all files in the specified folder.
    """
    files = []
    #MF6-3-M3235-23-017_S39_L001_R2_001.fastq.gz will turn into
    #MF6-3-M3235-23-017
    pattern = r'_S[0-9]+_L[0-9]+_R.*\.fastq\.gz'

    #traverse through the specified folder and all its subdirectories
    for root, _, filenames in os.walk(folder_path):
        for file in filenames:
            if os.path.isfile(os.path.join(root, file)):
                # Extract base name without matched pattern
                base_name = extract_base_name(file, pattern)
                if base_name:
                    files.append(base_name)
    return files


def extract_base_name(file_name, pattern):
    """
    Extracts the base name of the file without the matched pattern.
    """
    match = re.search(pattern, file_name)
    if match:
        return file_name[:match.start()]
    else:
        print(f'WARNING: not in the file format expected: {file_name}')

def parse_argument():

    parser = argparse.ArgumentParser(prog = 'create_report.py')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    #passed in as a string (space delimited) from command line
    parser.add_argument('-p', '--reports', metavar = '', required = True, help = 'Specify reports')
    parser.add_argument('-i', '--folder_path', metavar = '', required = True, help = 'Specify folder path for fasta.gz files')
    
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_argument()
    
    #read each report csv file as a df
    df_list = [pd.read_csv(report, index_col=[0]) for report in args.reports.split()]
        
    report_df = pd.concat(df_list)
    original_samples = list_files_in_folder(args.folder_path)
    noshow_samples = set(original_samples) - set(report_df.index)

    #add 'Mean read depth across entire run' as a separate line at the end
    report_df.loc['Mean read depth across entire run'] = [int(round(report_df['Mean read depth'].mean(),0)),None,None]
    
    for sample in noshow_samples:
        report_df.loc[f'{sample}'] = [None,None,None]

    report_df.to_csv(f"{args.output}")