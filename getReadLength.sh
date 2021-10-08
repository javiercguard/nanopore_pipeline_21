#!/bin/bash

# Very basic: creates a txt files with the list of read lengths from a FASTQ file
# $1: config file
# No need to queue this one

source $1

echo $sample_name

awk 'NR % 4 == 2 {print length($0)}' ${out_dir}/fastq/${sample_name}.merged.fastq > \
	${out_dir}/fastq/${sample_name}.readLength.txt