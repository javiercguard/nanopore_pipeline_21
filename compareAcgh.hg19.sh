#!/bin/bash

# $1: config file
# No need to queue this one

source $1
al=$2

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
else
 buildInfix="hg38"
fi

echo $sample_name

python3 findSV.py -id ${sample_name} -refVersion ${build} \
	-o ${out_dir}/variant_calling/${build} -bed ${out_dir}/coverage/${build}/input/${sample_name}.${buildInfix}.acgh.bed \
	-vcf ${out_dir}/variant_calling/${build}/${sample_name}.${al}.${buildInfix}.merged.vcf