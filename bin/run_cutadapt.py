#!/usr/bin/env python3
import sys, os, shutil, subprocess, argparse
from datetime import datetime
import concurrent.futures
import utilities, settings
import pandas as pd
from os.path import isfile
from os import access, R_OK

def cmd_exists(command):
	"""Checks for existence of a command on user's path

	Params
	------
	command: String
		Name of the command

	Returns
	------
	Path to executable command
	
	"""
	if shutil.which(command) is None:
		print(f"{command} not found on PATH")
		sys.exit()
	else:
		return(shutil.which(command))
 
def run_cutadapt(command):
	"""
	run cutadapt command

	Params
	------
	command: list (the actual run command in form of a list)
	
	"""

	# p = subprocess.run(command) 
	p = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	if p.returncode != 0:
		# print("cutadapt could not be executed.  Error encountered.")
		# print(p.returncode)
		print(p.stdout)
  
  
def parse_argument():
    # note
    # the script can be run as: python3 run_cutadapt.py -f read1.fq.gz -r read2.fg.gz -o out_dir
    #
    parser = argparse.ArgumentParser(prog = 'run_cutadapt.py')
    parser.add_argument('-f', '--read1', metavar = '', required = True, help = 'Specify R1 read')
    parser.add_argument('-r', '--read2', metavar = '', required = True, help = 'Specify R2 read')
    parser.add_argument('-o', '--out_dir', metavar = '', required = True, help = 'Specify output folder')
    parser.add_argument('-s', '--sample', metavar = '', required = True, help = 'Specify sample name')
    
    return parser.parse_args()


def remove_primer(sample,R1_gz,R2_gz,out_dir,num_process):
	"""
	TO DO...   add descriptions here...

	Params
	------
	config: config object
		Name of the command

	new_fasta: String
		name of the new fasta file (merged)

	num-process: int
		designated number of process to run the script

	"""

	#1. prep for cutadapt commands
	cutadapt_commands = []
	primers = utilities.Primers(settings.OLIGO_FILE)
	cutadapt_cmd = cmd_exists('cutadapt')
	for key in primers.pseqs:
		fprimer = primers.pseqs[key][0]
		rc_fprimer = utilities.revcomp(fprimer)
		rc_rprimer = primers.pseqs[key][1]
		rprimer = utilities.revcomp(rc_rprimer)
  
		cutadapt_commands.append([cutadapt_cmd, '-a', 
			f'{key}=^{fprimer}...{rc_rprimer}', '-A', 
       		f'{key}=^{rprimer}...{rc_fprimer}', '-o', 
			f'{out_dir}/{sample}.{key}.1.fastq', '-p', 
   			f'{out_dir}/{sample}.{key}.2.fastq',
			# f'{sample}.{key}.1.fastq', '-p', 
   			# f'{sample}.{key}.2.fastq',
			R1_gz, R2_gz,
			f"--rename='{{id}}  adapter={{adapter_name}} {{comment}}'",
       		'--quiet', '--discard-untrimmed', '-e', '0', '-m', '1', '-j', '1'])
    
    
	# print(datetime.now())
	# print(f"start with {sample} cutadapt")
	# print(len(cutadapt_commands))


	#2. run cutadapt 
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:
		executor.map(run_cutadapt, cutadapt_commands)
  	
	# print(datetime.now())
	# print(f"done with {sample} cutadapt")



def main():
    
    args = parse_argument()
    
    R1_gz = args.read1
    R2_gz = args.read2
    out_dir = args.out_dir
    sample = args.sample
    
    remove_primer(sample,R1_gz,R2_gz,out_dir,20)
            
            
if __name__ == "__main__":
    
    # settings.OLIGO_FILE is set in my local folder, you might need to update it
    if isfile(settings.OLIGO_FILE) and access(settings.OLIGO_FILE, R_OK):
        main()
    else:
        print (f"{settings.OLIGO_FILE} does not exist or is not readable" \
            	", please check it in the settings.py file ")
        sys.exit()
    
def create_primer_revcomp(primer_file):
    """
	This method read the oligos file and extract all the primer sequences and 
	write a new primer file, consisting of 5 columns (tab delimited)
	f-primer r-primer rc-f-primer rc-r-primer primer_name

	Params
	------
	primer_file: the path for the new primer file

	"""
    primers = utilities.Primers(settings.OLIGO_FILE)
    primers_list = []
    for key in primers.pseqs:
        fprimer = primers.pseqs[key][0]
        rc_fprimer = utilities.revcomp(fprimer)
        rc_rprimer = primers.pseqs[key][1]
        rprimer = utilities.revcomp(rc_rprimer)
        primers_list.append(f"{fprimer}\t{rprimer}\t{rc_fprimer}\t{rc_rprimer}\t{key}")
        
    with open (primer_file, 'w') as f:
        f.write('\n'.join(primers_list) + '\n')
        

# primer_file = '/scicomp/home-pure/qtl7/test/hmas_test/024_demulx_0mis_data/M3235_22_024.primers'
# create_primer_revcomp(primer_file)