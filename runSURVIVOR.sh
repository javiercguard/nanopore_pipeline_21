#!/bin/bash

# 
# This run SURVIVOR to merge variant callers' results and filters the output by AF (allele fraction).
# 

# 
# $1 : output dir
# $2 : genome name (either "GRCh37" or "GRCh38")
# $2..$n : paths to sample dirs, wich need to follow the structure created in nanopore_slurm_pipeline.py
# 

# If SURVIVOR is not in $PATH, uncomment and edit the line below, and comment the other one
# SURVIVOR=/path/to/SURVIVOR
SURVIVOR="SURVIVOR"

outputDir="$1"
shift
refGenome="$2"
shift

if [ $refGenome = "GRCh37" ]; then
 refInfix="hg19"
else
 refInfix="hg38"
fi

samples=$@

# For the inter-patient consensus
if [[ -f ${outputDir}/acghComparable.merged.minimap.${refInfix}.input ]] ; then rm ${outputDir}/acghComparable.merged.minimap.${refInfix}.input ; fi
# touch ${outputDir}/acghComparable.merged.minimap.${refInfix}.input
if [[ -f ${outputDir}/acghComparable.merged.ngmlr.${refInfix}.input ]] ; then rm ${outputDir}/acghComparable.merged.ngmlr.${refInfix}.input ; fi
# touch ${outputDir}/acghComparable.merged.ngmlr.${refInfix}.input

serpinPath=/rds/project/who1000/rds-who1000-wgs10k/WGS10K/data/projects/nanopore/us/analysis/PromethION/SERPINC1

# First we do the intra-patient consensus, and prepare the inter-patient input
for sample in ${samples[@]} ; do 
	# for al in "minimap" ; do
	for al in "minimap" "ngmlr"; do
		echo $sample $al
		sampleName=`echo $sample | cut -d"_" -f 1`
		# For the inter-patient consensus
		combinedSurvivorInput=${outputDir}/acghComparable.merged.$al.${refInfix}.input

		for file in `ls $sample/variant_calling/${refGenome}/{sniffles,svim}/*${refInfix}*$al*5*gz` ; do 
			gunzip < $file > `echo $file | sed 's/\.[^.]*$//'`
		done 
		(
			ls $sample/variant_calling/${refGenome}/svim/*${refInfix}*$al*5*vcf
			ls $sample/variant_calling/${refGenome}/sniffles/*${refInfix}*$al*5*vcf
		) > survivorInput.tmp # The Sniffles file has to go last, to not get overriden by SVIM BNDs
		$SURVIVOR merge survivorInput.tmp 500 1 0 0 0 1 $sample/variant_calling/${refGenome}/$sampleName.$al.${refInfix}.merged.vcf 
		rm survivorInput.tmp
		rm $sample/variant_calling/${refGenome}/{sniffles,svim}/*${refInfix}*$al*5*vcf
		
		grep -E "^#|SVLEN=-?([56789][[:digit:]]{4}|[[:digit:]]{6,});" \
		$sample/variant_calling/${refGenome}/${sampleName}.$al.${refInfix}.merged.vcf | \
		grep -E "^#|SVTYPE=(DUP|DEL)" | grep -E "^#|chr([[:digit:]]+|[XY])	" \
		> $sample/variant_calling/${refGenome}/${sampleName}.$al.merged.${refInfix}.acghComparable.vcf

		# echo "filtering by AF..."

		./filterByAf \
		$sample/variant_calling/${refGenome}/${sampleName}.$al.merged.${refInfix}.acghComparable.vcf > \
		$sample/variant_calling/${refGenome}/${sampleName}.$al.merged.${refInfix}.acghComparable.af20.vcf
	done 
done

echo "Merging patients"

Now we do the inter-patient consensus
for al in minimap ngmlr ; do
	combinedSurvivorInput=${outputDir}/acghComparable.merged.$al.${refInfix}.input
	$SURVIVOR merge ${combinedSurvivorInput} 500 1 0 0 0 1 ${outputDir}/acghComparable.$al.${refInfix}.merged.vcf
	rm ${combinedSurvivorInput}
done