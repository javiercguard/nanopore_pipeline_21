#!/usr/bin/env python3

# 
# This script takes a bed-like input with SVs (fourth column must be SVTYPE)
# And search for it in a VCF file
# 

# The vcf input needs to be sorted

import sys
import gzip
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description = "This script takes a bed-like input with SVs (fourth column must be SVTYPE), and search for them in a VCF file.")
parser.add_argument("-id", dest = "id", metavar = "ID", help = "ID of the sample, used to generate file names.")
parser.add_argument("-refVersion", dest = "refVersion", metavar = "GRChXX", help = "Reference version code, for file names.")
parser.add_argument("-o", dest = "outDir", metavar = "DIR", help = "Output dir.")
parser.add_argument("-bed", dest = "bed", metavar = "BED", help = "BED file with SVs to look for.")
parser.add_argument("-vcf", dest = "vcf", metavar = "VCF", help = "Sorted VCF file in which to look for variants.")
parser.add_argument("--min-size", dest = "minSize", metavar = "X", default = 0.7, help = "Candidate SVs must be at least (X * 100) % of searched SV's length. 0-1.")
parser.add_argument("--max-size", dest = "maxSize", metavar = "Y", default = 1.6, help = "Candidate SVs must not be at least (Y * 100) % of searched SV's length. 0-1.")
parser.add_argument("--dist-limit", dest = "distanceLimit", metavar = "N", default = 3e6, help = "Candidate SVs' coordinates must closer than N to target's coordinates. Accepts arabic numerals and e notation (1e2 for 100, etc.)")

args = parser.parse_args()

bed = args.bed
vcf = args.vcf

if args.refVersion == "GRCh37":
	refInfix = "hg19"
else:
	refInfix = "hg38"

args.minSize = float(args.minSize)
args.maxSize = float(args.maxSize)
args.distanceLimit = float(args.distanceLimit)

distanceLimit = args.distanceLimit


out = f"{args.outDir}/{args.vcf[args.vcf.rfind('/')+1:args.vcf.rfind('.vcf')]}.acghFound.vcf"

def getInfoValue (infoField, fieldName):
	result = infoField[infoField.index(fieldName) + len(fieldName) + 1:]
	result = result[0:result.index(";")] # This works even if there is not trailing semicolon
	return result

with (gzip.open(vcf, 'rt') if vcf.endswith(".gz") else open(vcf, 'rt')) as vcfFile, open(out, "w") as target:
	for line in vcfFile.readlines(): # lets print the header
		if (line[0] == "#"):
			target.write(line)
			pass
		else:
			break

moo = 0

with open(bed) as b, open(out, "a") as target:
	for sv in b.readlines():
		print(sv, file = sys.stderr, end = "")
		chrom, start, end, svType = sv.rstrip().split('\t') # coordinates of the sv we are looking for
		start, end = int(start), int(end)
		size = end - start
		foundChr = 0
		with (gzip.open(vcf, 'rt') if vcf.endswith(".gz") else open(vcf, 'rt')) as vcfFile:
			lines = vcfFile.readlines()
			for i in range(0, len(lines)):
				line = lines[i]
				if line[0] == "#":
					continue
				if not line.startswith(chrom + "\t"):
					if not foundChr:
						continue
					else:
						break
				foundChr = 1
				line = line.rstrip().split('\t')
				candidateChrom = line[0]
				candidateStart = int(line[1])
				if abs(candidateStart - start) > distanceLimit:
					continue
				info = line[7]
				candidateType = getInfoValue(info, "SVTYPE")
				if svType not in line[10] and svType not in line[9]:
					continue
				candidateLength = abs(int(getInfoValue(info, "SVLEN")))
				if not (candidateLength >= size * args.minSize and candidateLength <= size * args.maxSize):
					continue
				candidateEnd = int(getInfoValue(info, "END"))
				dist = ( abs(candidateStart - start) + abs(candidateEnd - end) ) / 2
				line[7] = line[7] + ("" if line[7].endswith(";") else ";") + f"DIST={dist};"
				target.write("\t".join(line) + "\n")
