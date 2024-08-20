#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import re
import yaml
import utilities

def make_primer_stats_yaml(output_file, primer_stats, oligos_file):
    '''
    this method generates a custom content yaml file specific for the multiqc report
    this yaml file is for the primer pair performance report

    Parameters
    ----------
    output_file: String, output file (yaml) name
    primer_stats: String, concatenated names of each primer_stats csv file (per sample)
    oligos_file: String, oligo file name(which contains the primer information)

    Returns: None
    ----------
    '''   
    def create_df(primer_stats, oligos_file):
        oligo_primers = utilities.Primers(oligos_file)
        # our primer panel as a dictionary, its value is a list of 3 items
        # running total: int
        # (min read count, sample name), tuple
        # (max read count, sample name), tuple
        primer_dict = {primer: [] for primer in oligo_primers.pnames}

        #read each report csv file as a df
        df_list = [pd.read_csv(report, index_col=[0], sep='\t') for report in primer_stats.split()]
        # each df has the following format (it has only 1 row)
        #                            primer1 primer2 primer3
        #sample(abundance count)     10      20      30
        for df in df_list:
            df_dict = df.to_dict(orient='list')
            sample = df.index.to_list()[0]
            #go through each primer in our original primer panel
            #check if it exists in the current primer_stats
            for primer in primer_dict:
                if primer in df_dict:
                    read_count = df_dict[primer][0]
                    if len(primer_dict[primer]) <= 0: #first timer
                        primer_dict[primer].append(read_count) #running total
                        primer_dict[primer].append((read_count,sample)) #min 
                        primer_dict[primer].append((read_count,sample)) #max
                    else:
                        primer_dict[primer][0] = primer_dict[primer][0] + read_count
                        if read_count < primer_dict[primer][1][0]:
                            primer_dict[primer][1] = (read_count,sample)
                        elif read_count > primer_dict[primer][2][0]:
                            primer_dict[primer][2] = (read_count,sample)

        #go through our primer panel again
        #if one primer has no associated reads at all, mark as 'n/a'
        #otherwise, calculate the mean reads count, and conver the min/max tuples to String
        for primer in primer_dict:
            if len(primer_dict[primer]) <= 0:
                primer_dict[primer].extend(['n/a']*3)
            else:
                primer_dict[primer][0] = primer_dict[primer][0]/len(df_list)
                primer_dict[primer][1] = f'{primer_dict[primer][1][0]} / {primer_dict[primer][1][1]}'
                primer_dict[primer][2] = f'{primer_dict[primer][2][0]} / {primer_dict[primer][2][1]}'

        df = pd.DataFrame(primer_dict).transpose()
        df.columns = ['p_col1', 'p_col2', 'p_col3']

        return df

    # Create headers dictionary
    headers = {
        'p_col1': {
            'title': 'average read count',
            'description': 'mean reads count per primer pair across all samples',
            'format': '{:,.1f}',
        },
        'p_col2': {
            'title': 'minimum read count',
            'description': 'minimum reads count for this primer pair across all samples (if there is a tie, only show the first one)',
            # 'format': '{:,.3f}',
            "scale": False
        },
        'p_col3': {
            'title': 'maxmium read count',
            'description': 'maximum reads count for this primer pair across all samples (if there is a tie, only show the first one)',
        },
    }

    # Convert the DataFrame to the required format
    data_yaml = create_df(primer_stats, oligos_file).to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': 'primer_report',
        'section_name': 'Primer performance report',
        'description': 'reads count report per primer pair across all samples in this run',
        'plot_type': 'table',
        'pconfig': {
            'id': 'primer_report',
            'sort_rows': False,
            'col1_header': 'Primer Name',
            "no_violin": True,
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(output_file, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)

        
def make_read_length_yaml(output, read_length, noshow_samples):
    '''
    this method generates a custom content yaml file specific for the multiqc report
    this yaml file is for the final combined reads length report

    Parameters
    ----------
    output: String, output file (yaml) name
    read_length: String, concatenated names of each read_length csv file
    noshow_samples: List, a list of sample names which does not generate any valid sequences

    Returns: None
    ----------
    '''   
    df_list = [pd.read_csv(report, index_col=[0], sep='\t') for report in read_length.split()]
    report_df = pd.concat(df_list)
    report_df.columns = ['l_col0','l_col1','l_col2','l_col3', 'l_col4']
    report_df.fillna('n/a', inplace=True)
    for sample in noshow_samples:
        report_df.loc[f'{sample}'] = [0, 0,'n/a', 'n/a', 'n/a']


    # Create headers dictionary
    headers = {
        'l_col0': {
            'title': 'total reads(non-unique) count',
            'description': 'total high quality reads count per sample across all primer pairs',
            "format": "{:,.0f}",
        },
        'l_col1': {
            'title': 'total reads(unique) count',
            'description': 'total high quality unique reads count per sample across all primer pairs',
            "format": "{:,.0f}",
        },
        'l_col2': {
            'title': 'average read length',
            'description': 'mean reads length per sample across all primer pairs',
            'format': '{:,.1f}',
            "scale": False,
        },
        'l_col3': {
            'title': 'minimum read length',
            'description': 'minimum reads length per sample across all primer pairs',
            "scale": False,
            "format": "{:,.0f}"
        },
        'l_col4': {
            'title': 'maximum read length',
            'description': 'maximum reads length per sample across all primer pairs',
            "scale": False,
            "format": "{:,.0f}"
        },
    }

    # Convert the DataFrame to the required format
    data_yaml = report_df.to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': 'read_length_report',
        'section_name': 'Sample read length report',
        'description': 'reads length report per sample across all primer pairs in this run',
        'plot_type': 'table',
        'pconfig': {
            'id': 'read_length_report',
            'sort_rows': False,
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(output, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)


def make_report_yaml(output_file, data_df):
    '''
    this method generates a custom content yaml file specific for the multiqc report
    this yaml file is for the final combined hmas summary report

    Parameters
    ----------
    output_file: output file (yaml) name
    data_df: data part of the yaml file in the format of dataframe

    Returns: None
    ----------
    '''   
    # Create headers dictionary
    headers = {
        'col1': {
            'title': 'Mean read depth',
            'description': 'we include only reads with at least 2 sequence count',
            'format': '{:,.1f}',
            "scale": False
        },
        'col2': {
            'title': '% of successful primer-pairs',
            'description': 'primer pairs with at lease 2 amplicons in the sample',
            'format': '{:,.3f}',
            "scale": False
        },
        'col3': {
            'title': '# of successful primer-pairs',
            'description': 'primer pairs with at lease 2 amplicons in the sample',
        },
    }

    # Convert the DataFrame to the required format
    data_yaml = data_df.to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': 'hmas_run_report',
        'section_name': 'HMAS run report',
        'description': 'combined summary report for all samples in this run',
        'plot_type': 'table',
        'pconfig': {
            'id': 'hmas_run_report',
            'sort_rows': False
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(output_file, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)


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
    parser.add_argument('-y', '--yaml', metavar = '', required = True, help = 'Specify output mqc report file')
    #passed in as a string (space delimited) from command line
    parser.add_argument('-p', '--reports', metavar = '', required = True, help = 'Specify reports')
    parser.add_argument('-i', '--folder_path', metavar = '', required = True, help = 'Specify folder path for fasta.gz files')
    
    parser.add_argument('-z', '--pyaml', metavar = '', required = True, help = 'Specify output primer_stats mqc report file')
    parser.add_argument('-q', '--primer_stats', metavar = '', required = True, help = 'Specify input primer_stats')  
    parser.add_argument('-l', '--primers', metavar = '', required = True, help = 'Specify primers')  

    parser.add_argument('-x', '--lyaml', metavar = '', required = True, help = 'Specify output read_length mqc report file')
    parser.add_argument('-r', '--read_length', metavar = '', required = True, help = 'Specify input read_length file')  

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

    # now generate report yaml file for MultiQC report
    # change column name
    report_df.columns = ['col1','col2','col3']
    #update empty cell to n/a
    report_df.fillna('n/a', inplace=True)
    make_report_yaml(args.yaml, report_df)
    make_primer_stats_yaml(args.pyaml, args.primer_stats, args.primers)
    make_read_length_yaml(args.lyaml, args.read_length, noshow_samples)