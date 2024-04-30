#!/usr/bin/env python3

import argparse, os, sys
import pandas as pd
from pathlib import Path
import logging
from configparser import ConfigParser
import run_blast as blast
import utilities
import settings


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(settings.log_handler)

# a list of control samples, which we will exclude from the analysis
# control_list = ['2014K_0979',
#                 'Water1']
control_list = []


def get_oligo_primer_dict():
    '''
    this method generates an instance of Primers class and a dictionary with primer name being the key
    and a default value 0

    Returns a dictionary

    '''
    oligo_primers = utilities.Primers(settings.OLIGO_FILE)
    primer_dict = {key: 0 for key in oligo_primers.pseqs}

    return primer_dict

def sample_primer_df_to_dict(df):
    '''
    this method converts a 2 column('sample', 'primer') dataframe into a dictionary, where
    sample is the key and primer is the value for that key

    sample    primer
    08_0810   OG0000294primerGroup8
    08_0810   OG0000296primerGroup7
    AR_0409   OG0000297primerGroup3
    AR_0409   OG0000298primerGroup4

    Parameters
    ----------
    df: the 2 columns dataframe

    Returns a dictionary
    {'08_0810':'OG0000294primerGroup8,OG0000296primerGroup7','AR_0409':'OG0000297primerGroup3,OG0000298primerGroup4'}

    '''
    df['primer'] = df.groupby(['sample'])['primer'].transform(lambda x: ','.join(x))
    return dict(df.values)

def map_sample_to_primer(csv_file):
    '''
    this method maps all the predicted primer pairs to the sample they're related with
    
    Parameters
    ----------
    csv_file: the path to the cvs file (metalsheet table, which should have 3 columns:
    #seq_id, primer, sample)

    Returns a dictionary {sample:[list of primers]}
    
    '''
    if csv_file:
        df = pd.read_csv(csv_file)
        df = df[['sample','primer']] # prep the df to be able to use sample_primer_df_to_dict
        sample_primer_dict = sample_primer_df_to_dict(df)
        # convert the dictionary value into a list
        new_sample_primer_dict = {key:value.split(',') for key,value in sample_primer_dict.items() if value}
        
        return new_sample_primer_dict
    
    
def csv_to_dict(csv_file):
    '''
    This method reads a standard csv file, uses the 1st column as index to creates a dataframe
    Then it converts the dataframe into a dictionary of dictionaries, with the index as key to 
    the first dictionary and column names of the dataframe as the keys for the 2nd dictionary(s)
    Ex.   this cvs file
    seq_id  primer  sample
    seq1    p1      s1
    seq2    p2      s2
    will turn into this structure:
    {'seq1':{'primer':'p1','sample':'s1'},
     'seq2':{'primer':'p2','sample':'s2'}}
    
    Parameters
    ----------
    cvs_file: the path to the cvs file 

    Returns a dictionary of dictionaries
    '''
    df = pd.read_csv(csv_file, index_col=0)
    map_dict = df.T.to_dict()

    return map_dict


def map_sample_to_isolate(map_file):
    '''
    this method reads mapping file (sample-to-isolate) in csv format and return a dictionary of it 
    key: sample (String)
    value: a list of corresponding isolates 
    (the list has no None value in it, and it's ensured to be a non-empty list)
    
    Parameters
    ----------
    map_file: the path to the mapping csv file 

    Returns a dictionary
    
    '''

    if map_file:
        df = pd.read_csv(map_file, index_col=0)
        # convert dataframe to a dictionary, index being key, and row being value in the form of a list
        map_dict = df.T.to_dict("list")
        # remove nan from the list
        for key in map_dict:
            map_dict[key] = [item for item in map_dict[key] if not(pd.isnull(item))]
        # remove keys which has empty value
        new_map_dict = {key:val for key,val in map_dict.items() if len(val) > 0}
        return new_map_dict
    

def rev_map_isolate_to_sample(map_dict):
    '''
    this method will reverse the above sample-to-isolate dictionary, and will generate a isolate-to-sample dictionary
    key: isolate (String)
    value: a list of corresponding samples 
    
    Parameters
    ----------
    map_dict: the mapping dictionary (sample-to-isolates) 

    Returns a dictionary
    
    '''
    rev_map_dict = {}
    for key in map_dict:
        isolates = map_dict[key]
        for isolate in isolates:
            if isolate not in rev_map_dict:
                rev_map_dict[isolate] = [key]
            else:
                if key not in rev_map_dict[isolate]:
                    rev_map_dict[isolate].append(key)
    return rev_map_dict


# helper method to read a full format count table, and/or save it as a pickle file for faster retrieval
def read_count_table(count_file):

    if Path(f"{count_file}.pkl").is_file():
        print(f"...loading {count_file}.pkl")
        df = pd.read_pickle(f"{count_file}.pkl")
    else:
        # read count_table
        print("loading csv count file")
        df = pd.read_csv(count_file, sep='\t')
        print("done loading csv count file")
        pd.to_pickle(df, f"{count_file}.pkl") 
        print("done dumping pkl file")

    return df

def create_blast_df(blast_file, query_fasta, reference, max_hit, metasheet_file):

    # the 'primer' is like: OG0000890-OG0000890primerGroup9-2014K_0979
    bcolnames = ["seq", "primer", "query_len", "subj_len", "len_aln", "eval", "cov", "pident", "mismatch"]

    if not blast_file: #if we don't already have blast result as a text file
        blast_file = blast.blast(query_fasta, reference, os.path.basename(query_fasta), max_hit)
        
    # blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames, index_col=False, header=None)
    blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames)
    
    # de-couple the formatting 
    # use the csv_to_dict() fuction to get the dict
    # '>OG0000294-OG0000294primerGroup8-ParatyphiA rc' - the seq_id will be a key of a dictionary
    # metasheet is a 3 column csv file, which contains explicit information of 
    # which primer pair and sample does an amplicon sequence correspond to. For example:
    # seq_id,primer,sample
    # OG0002941-OG0002941primerGroup0-2014K_0324,OG0002941primerGroup0,2014K_0324
    map_dict = csv_to_dict(metasheet_file)
    blast_df['primer_new'] = [ map_dict[seq_id]['primer'] for seq_id in blast_df['primer']]
    blast_df['sample'] = [ map_dict[seq_id]['sample'] for seq_id in blast_df['primer']]
    blast_df.drop(columns='primer', inplace=True)
    blast_df.rename(columns={'primer_new':'primer'}, inplace=True)
    blast_df['sample_primer'] = blast_df['sample'] + '.' + blast_df['primer']

    return blast_df

def blast_map_sample_to_isolate(blast_df, map_dict):
    '''
    this method helps to convert the 'isolate.primer' column in blast result into its corresponding
    'sample.primer' instead. One isolate might correspond to multiple samples, and this is taken
    care of by explode() method

    Parameters
    ----------
    blast_df: the original blast result df
    map_dict: the dictionary of isolate(key) to samples (a list of all corresponding samples)

    Returns the converted blast result dataframe
    '''
    def map_sample_to_isolate(x):
        sep = '.'
        isolate = x.split(sep)[0]
        primer_pair = x.split(sep)[1]
        sample_list = map_dict[isolate]
        return [f"{sample}{sep}{primer_pair}" for sample in sample_list]

    blast_df['sample_primer'] = blast_df['sample_primer'].apply(map_sample_to_isolate)

    return blast_df.explode('sample_primer', ignore_index=True)


def merge_count_blast(df, primer_list, blast_df):
    '''
    this method filters the original full-format count table (below), so that seqs 
    all have matches in blast result and their pident == 100 & cov >= 90 (or any value in the settings.py)

    Rep_Seq Sam1_PP1    Sam2_PP1    Sam1_PP2    Sam2_PP2
    Seq1    10  42  0   0
    Seq2    0   0   86  0
    Seq3    4   0   0   0

    Parameters
    ----------
    df: the original count table dataframe
    primer_list: the list of targeted sample.primer ['Sam1_PP1', 'Sam1_PP2'] etc.
    blast_df: blast result dataframe

    Returns the filtered dataframe

    '''

    #1.1 melt the df 
    df = df.melt(id_vars=['seq'],value_vars=primer_list,var_name='sample_primer_orig',value_name='count')
    # change dtype for 'count' to save memory
    df = df.astype({'count':int})
    df['sample_primer'] = df['sample_primer_orig']

    #1.2 filter and merge
    blast_df = blast_df[(blast_df['pident'] >= settings.PIDENT) & \
                        (blast_df['len_aln']/blast_df['subj_len'] >= settings.P_ALIGN) & \
                        (blast_df['cov'] >= settings.PCOV)]
    df = pd.merge(df,blast_df,on=['seq','sample_primer'])

    #1.3
    # duplicates come from the exploded blast_df, and they're all identical in 
    # the subset of ['seq','sample_primer_orig','count'], so we only need to keep one of them
    df.drop_duplicates(subset=['seq','sample_primer_orig'], inplace=True)

    #1.4 reverse the melt, and return to original format of df
    # df = df.drop(columns = ["query_len", "subj_len", "eval", "cov", "pident", "mismatch", "primer", "sample", "sample_primer_x","sample_primer_y"], errors='ignore')
    df = df[['seq','sample_primer_orig','count']]
    df = df.pivot(index='seq',columns='sample_primer_orig',values='count')

    return df


def split_samle_primer(sr, primers, sample_list, raw_idx):
    '''
    this method will convert:

    Sam1_PP1    2
    Sam2_PP1    1
    Sam1_PP2    1
    Sam2_PP2    0

    into:

        PP1 PP2
    Sam1    2   1
    Sam2    1   0

    Parameters
    ----------
    sr: the original dataframe Series
    primers: the list of intended column name ['PP1', 'PP2']
    sample_list: the list of samples we're studying
    raw_idx: index from the raw df, which is essentially the original columns from count_table

    Returns the converted dataframe

    '''

    #1. construct lists of count values(as columns) for the final df
    n_column_list = []
    _idx = sr.index
    for pp in primers:
        # check if sample.pp is in the index, otherwise make it blank
        # n_column_list.append([sr[f"{sample}.{pp}"] if f"{sample}.{pp}" in _idx else 0 for sample in sample_list])
        count_list = []
        for sample in sample_list:
            if f"{sample}.{pp}" in _idx:
                count_list.append(sr[f"{sample}.{pp}"])
            elif f"{sample}.{pp}" in raw_idx:
                count_list.append(0)
            else:
                count_list.append(None)
        n_column_list.append(count_list)

    #2. create the final_df from the list, transpose and reset the index name
    final_df = pd.DataFrame(n_column_list, index=primers).T
    final_df.index = sample_list

    return final_df

def parse_argument():
    # note
    # 1. the full format count table file is required !
    # 2. assume the sample/primer_group pair (as column name) in count table is in format of sample.primer_group
    # 3. sample_file is an one-column csv file containing the names of the sample we're analyzing
    #
    # the script can be run as: python3 parse_count_table.py -c juno.final.full.count_table -s M347-21-026-15sample.csv
    #                               -f juno.final.fasta -r juno_design_primers.fasta
    # - b / - f / - r are optional
    # parser = argparse.ArgumentParser(prog = 'pasr_count_table.py')
    # parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify count table')
    # parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    # parser.add_argument('-o', '--output', metavar = '', required = False, help = 'Specify output confusion-matrix file')
    # parser.add_argument('-b', '--blast', metavar = '', required = False, help = 'Specify blast result file')
    # parser.add_argument('-f', '--fasta', metavar = '', required = False, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    # parser.add_argument('-r', '--reference', metavar = '', required = False, help = 'Specify fasta file containing the positive control targets.')
    # parser.add_argument('-m', '--map_file', metavar = '', required = False, help = 'Specify mapping file')
    # parser.add_argument('-e', '--metasheet', metavar = '', required = False, help = 'Specify metasheet csv file')
    # parser.add_argument('-p', '--report_file', metavar = '', required = False, help = 'Specify report file name')

    # return parser.parse_args()

    parser = argparse.ArgumentParser(prog = 'pasr_count_table.py')
    parser.add_argument('-c', '--config_file', metavar = '', required = True, help = 'Specify configure file')
    config = ConfigParser()
    config.read(parser.parse_args().config_file)
    
    def dirFileExists(config,section_name,option_name):
        """
        Test if the option_name is in 'section_name' section and if it exists/readable
        """
        try:
            path = config[section_name][option_name]
            return os.path.exists(os.path.expanduser(path))
        except (KeyError): #option_name is not in 'file_inputs' section
            return False
        
    section = 'file_inputs' # we have only 1 section in the config file
    req_option_list = ['count_table','unique_fasta','reference_fasta','sample_list','metasheet','mapping_file']
    opt_option_list = ['output_file','report_file','logging_flag','blast_file']

    args = {} # dictionary holding all the inputs
    for option in req_option_list:
        if dirFileExists(config,section,option):
            args[option] = config[section][option]
        else:
            logger.error(f"{config[section][option]} does not exist")
            print (f"ERROR! {config[section][option]} does not exist")
            sys.exit(1)
    
    for option in opt_option_list:
        args[option] = config[section][option]
        
    return args

# helper method to split the column names, sort it
def get_column_list(df):

    column_list = list(set([i.split('.')[1] for i in df.index]))
    column_list.sort(key=str.lower)
    return column_list


def build_confusion_matrix(sample_list, raw_df_t, blast_df_t, map_dict, meta_sheet):
    '''
    this method flattens a list of sets into one single list and generate the occurrence count for each item
    as a dictionary

    Parameters
    ----------
    sample_list: a list of all samples in the data set
    blast_df_t: the transformed (and blast filtered) dataframe (row: primers, column: samples)
    map_dict: mapping dictionary between sample and isolate name 

    Returns
    ----------
    df_dict: a dictionary of TP/FP/TN/FN metrics (key: sample_name, value: a list of number of TP/FP/TN/FN counts )
    '''
    df_dict = {}
    for key in sample_list:
        primer_dict = get_oligo_primer_dict() # our full 2461 primer list
        # dictionary of sample:list_primer(predicted)
        sample_to_primer_dict = map_sample_to_primer(meta_sheet)
        primer_list = [] # hold the predicted primers for the sample
        
        
        
        # if map_dict and (key in map_dict): # map_dict is for mapping sample to isolates
        if key in map_dict:
            # Gut_10_3_1  Typhimurium	ParatyphiA might have 2 isolates
            for i in range(len(map_dict[key])):
                if map_dict[key][i] in sample_to_primer_dict:
                    primer_list.extend(sample_to_primer_dict[map_dict[key][i]])
                else:
                    logger.info(f"WARNING: {map_dict[key][i]} can't be found in the metasheet ")
        else:
            logger.info(f"WARNING: {key} can't be found in the sample-isolates mapping file")
                    
                    
        pred_pos_primer_set = set(primer_dict.keys()) & set(primer_list)
        pred_neg_primer_set = set(primer_dict.keys()) - pred_pos_primer_set
        
        obs_pos_primer_set = set(raw_df_t[key].dropna().index.tolist())
        obs_neg_primer_set = set(primer_dict.keys()) - obs_pos_primer_set
        
        # this is our true positives
        TP_primer_set = set(blast_df_t[key].dropna().index.tolist())
        # this is standard definition of true negatives
        TN_primer_set = obs_neg_primer_set & pred_neg_primer_set
        # FP is refined to those we falsely indentified and they're predicted to be negatives
        FP_primer_set = (obs_pos_primer_set - TP_primer_set) & pred_neg_primer_set
        # FN is made of 2 parts: those we didn't identify and those didn't meet blast threshold
        FN_primer_set = (obs_neg_primer_set - TN_primer_set) | (obs_pos_primer_set - TP_primer_set - FP_primer_set)
        
        df_dict[key] = list(map(len,(TP_primer_set,FP_primer_set,FN_primer_set,TN_primer_set)))
        
        logger.info(f"sample is: {key}")
        logger.info(f"FP are: {FP_primer_set}")
        logger.info(f"FN are: {FN_primer_set}")
        logger.info(f"TN are: {TN_primer_set}")
        logger.info(f"sanity check: {len(primer_dict)}, "
                    f"{len(primer_dict) == sum(df_dict[key])}")
        logger.info(f"sanity check 2, the following 3 sets SHOULD be emptry")
        logger.info(TP_primer_set & FP_primer_set)
        logger.info(TP_primer_set & FN_primer_set)
        logger.info(FN_primer_set & FP_primer_set)
        
    return df_dict


def create_report(raw_df, report_file):
    '''
    this method calculates metrics like: Mean read depth, # of failed primer pairs and 
    generate a report (csv file) for state public health lab

    Parameters
    ----------
    raw_df: a reformatted dataframe (row:samples, column:primer pairs) of original count_table
    report_file: the output file name

    Returns
    ----------
    None
    '''   
    report_df = raw_df.copy()
    
    total_primer_count = len(get_oligo_primer_dict())
    #this is failed primer pairs (common denominator) for all the samples
    all_failed_pp_count = total_primer_count - len(raw_df.columns)
    
    #including those empty cell while calculating the mean
    report_df['mean'] = report_df.fillna(0).mean(axis=1).apply(lambda x:(round(x,1)))
    report_df['failed_pp'] = report_df.isna().sum(axis=1).apply(lambda x: x+all_failed_pp_count)
    report_df['perc_successful_pp'] = round(1 - report_df['failed_pp']/total_primer_count, 3)
    report_df['failed_pp'] = report_df['failed_pp'].apply(lambda x: f"{x} / {total_primer_count}")
    
    report_df = report_df[['mean','perc_successful_pp','failed_pp']]
    #add 'Mean read depth across entire run' as a separate line at the end
    report_df.loc['Mean read depth across entire run'] = [int(round(report_df['mean'].mean(),0)),None,None]
    
    report_columns = [f'Mean read depth',
                      f'% of successful primer-pairs\n(has at least 10 amplicons across all samples)',
                      f'# of primer pairs with less than 10 amplicons mapping\nover total primer-pairs']
    report_df.columns = report_columns
    
    report_df.to_csv(f'{report_file}.csv') 
    
def main():

    args = parse_argument()
    
    # default is INFO:20
    if args['logging_flag']:
        logger.setLevel(int(args['logging_flag']))

    df = read_count_table(args['count_table'])

    # read sample list file
    sample_list = pd.read_csv(args['sample_list'], names = ['sample'])['sample'].tolist()
    sample_list = [sample for sample in sample_list if sample not in control_list]
    sample_list.sort(key=str.lower)
    print (sample_list)

	#filter columns to contain only sample names in the sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    seq_df = df.iloc[:,0] # keep the 1st column (sequence)
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]
    df['seq'] = seq_df
    raw_idx = df.sum().index

    raw_df = df.sum() #total abundance all high quality seqs
    raw_df.drop(['seq'], inplace=True)
    raw_df = split_samle_primer(raw_df, get_column_list(raw_df), sample_list, raw_idx)

    #generate report for SPHL
    if args['report_file']:
        create_report(raw_df, args['report_file'])

    #print out raw amplicon sequence abundance info 
    if logger.isEnabledFor(logging.DEBUG):
        raw_df['mean'] = raw_df.fillna(0).mean(axis=1)
        raw_df.to_csv('raw_df', sep='\t') # save as a tsv file

    #creat the blast df, to blast filtering all sequences
    blast_df = create_blast_df(args['blast_file'], args['unique_fasta'], args['reference_fasta'], 100, args['metasheet'])
    # mapping dictionary between sample and isolate
    map_dict = map_sample_to_isolate(args['mapping_file'])
    
    # mapping (samples-isolates) is required on 02/15/2023
    # create a reverse mapping between isolate and samples, it runs faster this way
    rev_map_dict = rev_map_isolate_to_sample(map_dict)
    # create a new blast df based on the above mapping
    blast_df = blast_map_sample_to_isolate(blast_df, rev_map_dict)

    # the filtered version of our count table df
    df = merge_count_blast(df, primer_list, blast_df)
    blast_df = df.sum()
    blast_df = split_samle_primer(blast_df, get_column_list(blast_df), sample_list, raw_idx)
    
    #print out blast filtered amplicon sequence abundance info
    if logger.isEnabledFor(logging.DEBUG):
        blast_df['mean'] = blast_df.fillna(0).mean(axis=1)
        blast_df.to_csv('blast_df', sep='\t') # save as a tsv file

    raw_df_t = raw_df.T
    blast_df_t = blast_df.T
    df_dict = build_confusion_matrix(sample_list, raw_df_t, blast_df_t, map_dict, args['metasheet'])
    df = pd.DataFrame.from_dict(df_dict, orient='index')
    df.columns = ['TP','FP','FN','TN']
    df['sensitivity (TP/P)'] = df['TP']/(df['TP'] + df['FN'])
    df['sensitivity (TP/P)'] = df['sensitivity (TP/P)'].apply(lambda x:round(x,3))
    df['specificity (TN/N)'] = (df['TN']/(df['TN'] + df['FP'])).apply(lambda x:round(x,3))
    df['precision (TP/TP+FP)'] = (df['TP']/(df['TP'] + df['FP'])).apply(lambda x:round(x,3))
    df['ACC'] = (df['TP'] + df['TN'])/(df['TP'] + df['TN'] + df['FP'] + df['FN'])
    df['ACC'] = df['ACC'].apply(lambda x:round(x,3))

    # if we need to generate output file
    if args['output_file']:
        df.to_csv(args['output_file'], sep='\t')

if __name__ == "__main__":
    main()