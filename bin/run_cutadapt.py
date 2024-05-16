#!/usr/bin/env python3
import sys, os, shutil, subprocess, argparse
import concurrent.futures
import utilities
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
    # oligo_file contains the primer information (4 columns tab delimited), Ex:
    # primer	TTATCGGGATGCCAGATCTGC	GRCGGGGACATTCTCCTCCAG	OG0003222primerGroup
    # primer / forward_seq / rc_reverse_seq / primer_name
    parser = argparse.ArgumentParser(prog = 'run_cutadapt.py')
    parser.add_argument('-f', '--read1', metavar = '', required = True, help = 'Specify R1 read')
    parser.add_argument('-r', '--read2', metavar = '', required = True, help = 'Specify R2 read')
    parser.add_argument('-o', '--out_dir', metavar = '', required = True, help = 'Specify output folder')
    parser.add_argument('-s', '--sample', metavar = '', required = True, help = 'Specify sample name')
    parser.add_argument('-p', '--oligo_file', metavar = '', required = True, help = 'Specify oligo file')
    parser.add_argument('-e', '--max_error', metavar = '', required = True, help = 'max error')
    parser.add_argument('-l', '--min_length', metavar = '', required = True, help = 'min length')
    parser.add_argument('-t', '--thread', metavar = '', required = True, help = 'thread')
    parser.add_argument('-b', '--mode', metavar = '', required = True, help = 'boolean value if reads are shorter than amplicons')
    
    
    return parser.parse_args()

def remove_primer(args):
	"""
	takes in all arguments and run cutadapt to remove primres concurrently

	Params
	------
	args: argparse.Namespace object
		holding all the arguments

	"""
	sample = args.sample
	R1_gz = args.read1
	R2_gz = args.read2
	out_dir = args.out_dir
	oligo_file = args.oligo_file
	max_error = args.max_error
	min_length = args.min_length
	thread = args.thread
	flag = args.mode.lower() == 'true' #convert string to boolean
    
	#1. prep for cutadapt commands
	# cutadapt_commands = []
	primers = utilities.Primers(oligo_file)
	cutadapt_cmd = cmd_exists('cutadapt')

	cutadapt_commands = [cutadapt_cmd, '-o',
		f'{out_dir}/{sample}.1.fastq', '-p', 
		f'{out_dir}/{sample}.2.fastq',
		R1_gz, R2_gz,
		f"--rename={{id}}  adapter={{adapter_name}}={sample} {{comment}}",
		'--quiet', '--discard-untrimmed', '-e', f'{max_error}', '-m', f'{min_length}', '-j', f'{thread}']

	# go through all our primer pairs and concatenate all the primer sequences in a single cutadapt command
	for key in primers.pseqs:
		fprimer = primers.pseqs[key][0]
		rc_fprimer = utilities.revcomp(fprimer)
		rc_rprimer = primers.pseqs[key][1] # by default we use reverse complement of reverse_primer in primer file
		rprimer = utilities.revcomp(rc_rprimer)

		if flag: # reads longer than amplicons (use a linked adapter)
			cutadapt_commands.insert(1,f'{key}=^{rprimer}...{rc_fprimer}')
			cutadapt_commands.insert(1,'-A')
			cutadapt_commands.insert(1,f'{key}=^{fprimer}...{rc_rprimer}')
			cutadapt_commands.insert(1,'-a')
		else: # reads shorter than amplicons
			cutadapt_commands.insert(1,f'{key}=^{rprimer}')
			cutadapt_commands.insert(1,'-G')
			cutadapt_commands.insert(1,f'{key}=^{fprimer}')
			cutadapt_commands.insert(1,'-g')
		

	run_cutadapt(cutadapt_commands)


            
if __name__ == "__main__":
    
    args = parse_argument()
    oligo_file = args.oligo_file
    
    if isfile(oligo_file) and access(oligo_file, R_OK):
        remove_primer(args)
    else:
        print (f"{oligo_file} does not exist or is not readable" \
            	", please check it in the settings.py file ")
        sys.exit()