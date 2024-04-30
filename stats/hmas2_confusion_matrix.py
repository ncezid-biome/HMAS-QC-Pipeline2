import os
import configparser
import subprocess
from concurrent.futures import ProcessPoolExecutor
import argparse
import csv
import pandas as pd

# Define the fixed section for the config.ini file
fixed_section = 'file_inputs'


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input file folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify the combined confusion_matrix file')
    parser.add_argument('-r', '--reference', metavar = '', required = True, help = 'reference_fasta file')
    parser.add_argument('-m', '--metasheet', metavar = '', required = True, help = 'metasheet file')
    parser.add_argument('-p', '--mapping', metavar = '', required = True, help = 'sample isolates mapping file')
    parser.add_argument('-s', '--script', metavar = '', required = True, help = 'python confusio_matrix script')
    return parser.parse_args()


def get_folders(args):
    
    # Get all folder names in the input folder
    folder_names = [f for f in os.listdir(args.input) if os.path.isdir(os.path.join(args.input, f))]
    
    return folder_names


def get_fixed_options(args):
    
    # Define the fixed options for the config.ini file
    fixed_options = {
        'reference_fasta': args.reference,
        'metasheet': args.metasheet,
        'mapping_file': args.mapping,
        # if necessary to generate SPHL report
        'report_file': '',  # CAN NOT set it to None !!
        # specify logging flag to output raw_df/blast_df (set it to DEBUG:10 to output raw_df/blast_df)
        # default is INFO: 20, MUST use the int value
        'logging_flag': '20', # CAN NOT use int here !!
        'blast_file': ''
    }
    
    return fixed_options

# Define a function that creates a config.ini file and runs another script with it
def process_folder(parameters):
    
    folder_name, fixed_options, args = parameters
    
    # Create a new configparser instance
    config = configparser.ConfigParser()
    config[fixed_section] = {} #this step seems to be a must ?

    # Set the fixed options for the config.ini file
    for key, value in fixed_options.items():
        config.set(fixed_section, key, value)

    # Set the customized options for the config.ini file based on the folder name
    # count_table might be saved in temp folder
    count_table_file = f'{args.input}/{folder_name}/{folder_name}.final.count_table'
    if not os.path.isfile(count_table_file):
        count_table_file = f'{args.input}/{folder_name}/temp/{folder_name}.final.count_table'
        
    config.set(fixed_section, 'count_table', count_table_file)
    config.set(fixed_section, 'unique_fasta', f'{args.input}/{folder_name}/{folder_name}.final.unique.fasta')
    
    # create a sample.csv for this folder
    csv_file_path = f'{args.input}/{folder_name}/{folder_name}.sample.csv'
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow([f'{folder_name}'])
        
    config.set(fixed_section, 'sample_list', csv_file_path)
    config.set(fixed_section, 'output_file', f'{args.input}/{folder_name}/confusion_matrix_{folder_name}')

    # Save the config.ini file in the current folder
    with open(f'{args.input}/{folder_name}/config.ini', 'w') as configfile:
        config.write(configfile)

    # Execute another Python script with the created config.ini file
    subprocess.run(['python3', args.script, '-c', f'{args.input}/{folder_name}/config.ini'])
    confusion_matrix_df = pd.read_csv(f'{args.input}/{folder_name}/confusion_matrix_{folder_name}', index_col=[0])

    # Delete the config.ini file and the sample.csv file
    os.remove(f'{args.input}/{folder_name}/config.ini')
    os.remove(csv_file_path)
    
    return confusion_matrix_df

if __name__ == "__main__":
    
    args = parse_argument()
    folder_names = get_folders(args)
    fixed_options = get_fixed_options(args)
    
    # # Create a ThreadPoolExecutor with max_workers set to the number of folders
    #I tried to run multiple processes in parallel, but ran into errors related to blast ?
    #so we run them one by one for now
    df_list = []
    for folder_name in folder_names:
        #skip those isolates without any valid unique sequences
        if not os.path.exists(f'{args.input}/{folder_name}/{folder_name}.final.unique.fasta') or os.path.getsize(f'{args.input}/{folder_name}/{folder_name}.final.unique.fasta') <= 0:
            continue
        parameters = (folder_name, fixed_options, args)
        df_list.append(process_folder(parameters))
        
    #concat all individual confusion_matrix output into a single one
    confusion_matrix = pd.concat(df_list)
    confusion_matrix.to_csv(f'{args.output}')






# count_table = /scicomp/home-pure/qtl7/test/hmas_test/024_demulx_0mis_data/output_samplebase_24_M3235_23_008/Ctrl-SEL/Ctrl-SEL.final.count_table
# unique_fasta = /scicomp/home-pure/qtl7/test/hmas_test/024_demulx_0mis_data/output_samplebase_24_M3235_23_008/Ctrl-SEL/Ctrl-SEL.final.unique.fasta
# sample_list = /scicomp/home-pure/qtl7/HMAS_QC_Pipeline/M3235-23-008_IRP_enrichment.sample.csv
# reference_fasta = /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/19isolates/primersearch/2011K_0222_extractedAmplicons.fasta
# metasheet = /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/19isolates/primersearch/2011K_0222_metasheet.csv
# mapping_file = /scicomp/home-pure/qtl7/HMAS_QC_Pipeline/M3235-23-008_IRP_enrichment_sample_isolates_mapping.csv
