#!/bin/bash

# $1					config file
# $2 					aligner to use

# Can be called like this (grep will depend on sample naming):
# for file in $(find ./config -type f | grep "\./config/[[:digit:]]"); \
# do for al in minimap ngmlr ; do bash runNanoStat.sh $file $al ; done ; done

source $1

sbatch --error=${out_dir}/log/nanostat_${2}.err --output=${out_dir}/log/nanostat_${2}.out \
	--parsable slurmScripts/runNanoStats.slurm $@