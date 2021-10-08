#! /bin/bash

# Queues slurmScripts/runMosdepthComplete.slurm,
# which generates the handy statistics mosdepth provides
# $1					config file
# $2 					aligner to use
# $n, $n+1, and so on	BED with variants, infix

source $1

sbatch --error=${out_dir}/log/mosdepth_${2}.err --output=${out_dir}/log/mosdepth_${2}.out --parsable ng_slurm_commandsRedux/runMosdepthComplete.slurm $@