#!/bin/bash

#SBATCH -J Sniffles
#SBATCH --nodes=3
#SBATCH --ntasks=48
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

#! Enable conda usage, particularly my Sniffles env
source $conda_path

#! Define variables
aligner=$2
threads=46

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
else
 buildInfix="hg38"
 ref=$build_38
fi

bamFile="${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam"
vcfDir="${out_dir}/variant_calling/${build}/sniffles"

echo ${vcfDir}

vcfTmp="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min1.unsorted.tmp.vcf"

#! Make directories
mkdir -p "${out_dir}/variant_calling/${build}/sniffles"

#! Run commands
echo "Starting SNIFFLES at:"
date

/rds/user/bd417/hpc-work/Sniffles-master/bin/sniffles-core-1.0.12/sniffles -t ${threads} -s 1 -m ${bamFile} -v ${vcfTmp}

conda activate ${env_path}

# Sort the Sniffles output
sortedVcf=${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min1.sorted.tmp.vcf
bcftools sort ${vcfTmp} -O v -o ${sortedVcf}

# Filter by RE > 2, and whatever support was listed in the config
# These come out sorted
for support in "2" "${sr_sniffles}" ; do 

	grep -E "^#|RE=([$(seq -s '' $support 9)]{1};?|[[:digit:]]{2,};?)" \
	${sortedVcf} \
	> "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}.sorted.tmp.vcf"

done

# Compress and index the vcfs, and remove the temporary ones
for support in "1" "2" "${sr_sniffles}" ; do 

	definitiveVcf="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}.vcf.gz"

	bcftools view -O z \
	"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}.sorted.tmp.vcf" \
	-o "${definitiveVcf}"
	bcftools index "${definitiveVcf}"

	if [[ -f "${definitiveVcf}" ]]; then
		rm "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}."*".tmp.vcf"
	fi

done

echo "Finishing SNIFFLES at:"
date

#! Deactivate ONT environment
conda deactivate
