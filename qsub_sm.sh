#!/bin/bash
source /etc/profile

# The above lines specify the shell that will be used to execute the script and sources the global profile settings on the system, ensuring access to environment variables and system-wide configurations. Users don't usually need to change this.

# Set the path to the standard output file for the job. Change to a different file path if desired.
#$ -o qsub_hmas.out
# Set the path to the standard error file for the job. Change to a different file path if desired.
#$ -e qsub_hmas.err
# Set the name of the job. Change to a different name that describes your job.
#$ -N qsub_hmas
# Request a parallel environment (PE) called "smp" with 8 slots/processors for the job. Adjust "8" to request a different number of processors.
#$ -pe smp 8
# Set the hard runtime limit for the job. Change to specify a different maximum runtime.
#$ -l h_rt=96:00:00
# Set the hard virtual memory limit for the job. Change to request a different amount of memory.
#$ -l h_vmem=16G
# Specify the queue where the job will be submitted. Change "all.q" to a different queue if available.
#$ -q highmem.q
# Run the script in the current working directory where the qsub command is executed using "-cwd", or specify a directory with "-wd path/to/working/directory"
#$ -cwd
# Set the mail notification preferences for the job with "-m abe -M user@cdc.gov". "abe" stands for "a" (job is aborted), "b" (job begins), and "e" (job ends).
#$ -m abe -M qtl7@cdc.gov

# ------------------------------------------------------------------------------------------------------------------

# Load conda if needed (adjust path as needed)
# ls ~/miniconda3/etc/profile.d/conda.sh  or: which conda (check where is the conda in your system)
source ~/miniconda3/etc/profile.d/conda.sh

# Activate conda environment (adjust name as needed)
conda activate hmas_mqc_07132024

nextflow run hmas2.nf \
	--primer /scicomp/groups-pure/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/step_mothur/HMAS-QC-Pipeline2/M3235_22_024.oligos \
	--reads /scicomp/groups-pure/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/prod/Incoming2/HMAS_101_KY_PRJNA23040H \
	--outdir HMAS_101_KY


