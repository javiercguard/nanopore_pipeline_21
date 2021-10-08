#!/bin/bash

#SBATCH -J MERGE-FASTQS 
#SBATCH --nodes=1
#SBATCH --ntasks=1
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

#! Source config file
source $1

#! Define variables 
fastq=${out_dir}/fastq/${sample_name}.fastq.fof
merged=${out_dir}/fastq/${sample_name}.merged.fastq

#! Run commands
echo "Starting merging FASTQs at:"
date

echo "Merging fastqs"
cat $fastq | xargs cat > $merged

echo "Finishing merging at:"
date
