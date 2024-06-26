#!/usr/bin/env python3

import argparse
import csv

'''
This script reads all report.csv (for each sample) and concantenate them into one single report.csv. It will 
also add a 'Mean read depth across entire run' row at the end.

This script requires all report.csv file names be passed in a concatenated string as a command line argument.
'''

# Function to parse the input and extract values
def parse_input(input_text):
    lines = input_text.strip().split('\n')

    total_reads = 0
    discarded_reads = 0

    for line in lines:

        if "sequences kept" in line:
            parts = line.split()

            # Extracting and validating the required values
            total_reads = parts[0]
            discarded_reads = parts[-3]

            return [total_reads, discarded_reads]

    return [total_reads, discarded_reads]

# Function to write values to a CSV file
def write_to_csv(sample_name, values, output_file):
    header = ['Sample name', 'Total reads', 'Discarded reads']
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow([sample_name] + values)

def parse_argument():

    parser = argparse.ArgumentParser(prog = 'parse_qfilter_log.py')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    parser.add_argument('-p', '--log', metavar = '', required = True, help = 'Specify log file')
    parser.add_argument('-s', '--sample', metavar = '', required = True, help = 'Specify sample name')
    
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_argument()
    
    # Read from the given input file
    with open(args.log, 'r') as file:
        input_text = file.read()

    # Write the values to a CSV file
    write_to_csv(args.sample, parse_input(input_text), args.output)
