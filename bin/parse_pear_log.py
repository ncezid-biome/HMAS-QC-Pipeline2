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
    
    assembled_reads = 0
    discarded_reads = 0
    unassembled_reads = 0

    for line in lines:
        if "Assembled reads ...................:" in line:
            parts = line.split(':')[-1].split('/')
            assembled_reads = int(parts[0].strip().replace(',', ''))
        elif "Discarded reads ...................:" in line:
            discarded_reads = int(line.split(':')[-1].split('/')[0].strip().replace(',', ''))
        elif "Not assembled reads ...............:" in line:
            unassembled_reads = int(line.split(':')[-1].split('/')[0].strip().replace(',', ''))

    # return [total_reads, assembled_reads, discarded_reads, unassembled_reads]
            return [assembled_reads, discarded_reads, unassembled_reads]

# Function to write values to a CSV file
def write_to_csv(sample_name, values, output_file):
    header = ['Sample name', 'Assembled reads', 'Discarded reads', 'Un-assembled reads']
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerow([sample_name] + values)

def parse_argument():

    parser = argparse.ArgumentParser(prog = 'parse_pear_log.py')
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
