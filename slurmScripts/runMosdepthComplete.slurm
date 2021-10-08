#!/bin/bash

#SBATCH -J MOSDEPTH
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --exclusive

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Modify the environment seen by the application. For this example we need the default modules.
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
module load miniconda/3

## Source the config file
source $1
shift

#! Set ONT enviroment
source $conda_path
mosdepth_env="/home/bd417/.conda/envs/mosdepth"
conda activate ${mosdepth_env}

#! Define variables
threads=4
aligner=$1
shift

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
 refBed=${ref_bed_37}
else
 buildInfix="hg38"
 ref=$build_38
 refBed=${ref_bed_38}
fi

bamFile=${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam
inputDir=${out_dir}/coverage/${build}/input
outputDir=${out_dir}/coverage/${build}/output

#! Make directories
## ${out_dir}/coverage/${build}/input must exist beforehand
mkdir -p ${outputDir}

#! Run commands
echo "Starting GETCOVERAGE at:"
date

while [ -n "$2" ]; do
	file=${inputDir}/"$1"
	shift
	infix=${aligner}."$1"
	shift

	echo "Running mosdepth"

	  mosdepth -t ${threads} -n ${outputDir}/${sample_name}.${infix}.tmp ${bamFile}

	echo "Finished mosdepth"

	# Remove the tmp files

	rm ${outputDir}/${sample_name}.${infix}.tmp*

done

#! Deactivate ONT environment
conda deactivate