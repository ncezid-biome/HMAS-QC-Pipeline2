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

def parse_argument():

    parser = argparse.ArgumentParser(prog = 'create_report.py')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')    #passed in as a string (space delimited) from command line
    parser.add_argument('-p', '--reports', metavar = '', required = True, help = 'Specify reports')
    
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_argument()
    
    #read each report csv file as a df
    df_list = [pd.read_csv(report, index_col=[0]) for report in args.reports.split()]
        
    report_df = pd.concat(df_list)
    report_df.to_csv(f"{args.output}", sep='\t')