#!/usr/bin/env python3

import os, sys, shutil, subprocess
from pathlib import Path

"""
This script is mostly from the blast portion of hmas_validation.py by Jessica Rowell
I modify it to make it a stand-alone script.
to run blast, simply call: blast(q_fasta, r_fasta, out_file, max_hits)
"""

def file_exists(filename):
    """Checks for existence of a file

    Params
    ------
    file: String
        Name of the file

    Returns
    ------
    Path to file

    """
    if Path(filename).is_file():
        print(f"{filename} found.")
        return Path(filename).resolve()
    else:
        print(f"{filename} could not be found.")
        sys.exit()

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
        return shutil.which(command)


def make_blast_db(reference, dbtype, makeblastdb):
    """Makes a Blast database from a set of reference amplicons

    Params
    ------
    reference: String
        Path to reference amplicon set

    dbtype: String
	Specify nucleotide or protein database

    Returns
    ------

    """
    db_name = os.path.splitext(reference)[0]
    p = subprocess.run([makeblastdb, '-in', reference, '-out', db_name, '-dbtype', dbtype, '-parse_seqids'])
    if p.returncode != 0:
        print("makeblastdb could not be executed.  Error encountered.")
        print(p.returncode)
    else:
        print("makeblastdb run successfully.")
        return db_name


def run_blast(db, fasta, outfile, blastn, maxhits=10):
    """Runs Blast on a fasta against a reference database and saves output

    Params
    ------
    ref: String
        Path to reference blast database

    fasta: String
        Specify query fasta sequence

    outfile: String
        Specify output file name

    Returns
    ------
    blast_out: String
        Name of the file containing the blast results
    """
    blast_out = outfile + '_blast.out'
    if maxhits != 10:
        try:
            maxhits + 1 - 1
        except TypeError:
            print(f"You must supply an integer value for max. number of hits. The default is 10.")
            maxhits = 10
        finally:
            maxhits = maxhits
    print(f"Running blast with {maxhits} maximum hits to generate")
    p = subprocess.run(
        [blastn, '-db', db, '-query', fasta,
         '-outfmt',
         # '6 qseqid sseqid qlen slen length evalue qcovs pident mismatch',
         # '-max_target_seqs', str(maxhits), '-out', blast_out])
         '6 qseqid sseqid qlen slen length evalue qcovs pident mismatch',
         '-max_target_seqs', str(maxhits), '-max_hsps', '1', '-out', blast_out])
    if p.returncode != 0:
        print("blastn could not be executed.  Error encountered.")
        print(p.returncode)
    else:
        print("blastn run successfully.")
        print(
            "format: query ID, subject ID, query length, subj length, e val, query cov/subj, percent identical matches, # mismatches")
        return blast_out

def blast(q_fasta, r_fasta, out_file, max_hits):

    query_fasta = file_exists(q_fasta)
    reference_fasta = file_exists(r_fasta)

    path_to_makeblastdb = cmd_exists('makeblastdb')
    path_to_blastn = cmd_exists('blastn')

    db = make_blast_db(reference_fasta, 'nucl', path_to_makeblastdb)
    return run_blast(db, query_fasta, out_file, path_to_blastn, max_hits)