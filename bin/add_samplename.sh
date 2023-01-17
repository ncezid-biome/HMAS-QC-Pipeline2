#/bin/bash

sample_fastq=$1
sample_name=$2

sed -i -e "s/\(^@M03235.*\) .*$/\1=${sample_name}/" ${sample_fastq} 
