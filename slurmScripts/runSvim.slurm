#!/bin/bash

#SBATCH -J SVIM
#SBATCH --nodes=1
#SBATCH --ntasks=5
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
#!source ~/.bashrc

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
vcfDir="${out_dir}/variant_calling/${build}/svim"
vcfTmp="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.unsorted.tmp.vcf"

#! Make directories
mkdir -p "${vcfDir}"
mkdir -p "${vcfDir}/${sample_name}_svim" # Technically redundant, but more explicit

#! Run commands
echo "Starting SVIM at:"
date

conda activate /home/bd417/.conda/envs/svim

svim alignment "${vcfDir}"/"${sample_name}"_svim/ "${bamFile}" "${ref}"

conda deactivate

cp "${vcfDir}"/"${sample_name}"_svim/variants.vcf "${vcfTmp}"

conda activate ${env_path}

#! Sort the SVIM output
bcftools sort "$vcfTmp" -O v \
-o "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.sorted.tmp.vcf"

#! Filter by QUAL > support in the config. Comes out sorted.

for qual in "5" "10" ; do

	awk -v 'OFS=\t' '{if($0 ~ /^#/ || $6 >= '${qual}')print}' \
	"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.sorted.tmp.vcf" \
	> "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual${qual}.sorted.tmp.vcf"

done

#! Compress and index the vcfs, and remove the temporary ones. The naming convention make looping pointless

bcftools view -O z \
"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.sorted.tmp.vcf" \
-o "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.vcf.gz"
bcftools view -O z \
"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual5.sorted.tmp.vcf" \
-o "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual5.vcf.gz"
bcftools view -O z \
"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual10.sorted.tmp.vcf" \
-o "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual10.vcf.gz"

bcftools index "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.vcf.gz"
bcftools index "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual5.vcf.gz"
bcftools index "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual10.vcf.gz"

if [[ -f "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.qual10.vcf.gz" ]]; then
	rm "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim."*".tmp.vcf"
	# rm -r "${vcfDir}"/"${sample_name}"_svim # Comment/Uncomment to keep/delete SVIM generated files
fi

echo "Finishing SVIM at:"
date

#! Deactivate ONT environment
conda deactivate