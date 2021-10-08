#!/bin/bash

#SBATCH -J NanoStats
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

#! Set ONT enviroment
source $conda_path
nanostat_env="/home/bd417/hpc-work/.envs"
conda activate ${nanostat_env}

#! Define variables
threads=14
aligner=$2

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
outputDir=${out_dir}/nanostats/${build}

#! Make directories
mkdir -p ${outputDir}

#! Run commands
echo "Starting nanostat at:"
date

NanoStat -n ${sample_name}.${buildInfix}.${aligner}.nanostat.txt -o ${outputDir} \
	-t ${threads} --bam ${bamFile}

#! Deactivate conda environment
conda deactivate
