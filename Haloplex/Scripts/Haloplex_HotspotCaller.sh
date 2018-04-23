#!/bin/bash

GATKv3_8="../../../Software/GenomeAnalysisTK_v3.8.jar"
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
dbsnp="../../../hg38.p6/All.vcf"
hotspots="../Dependent_Files/hotspots.list"
java="../../../Software/jre1.8.0_112/bin/java"

# COMMANDLINE VARIABLES

while getopts "f:h" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		h)
			echo "Usage: $0 [-f FILEPATH (optional)] " >&2
			echo
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
			echo "	-h		display this help message"
			exit 1
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac
done

if [ -z $filepath ]; then
	filepath=`pwd`"/"
fi

if [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

# SCRIPT

mkdir ${filepath}hotspot_bam
mkdir ${filepath}hotspot_vcf

while read pfx; do
	$java -jar $GATKv3_8 -T MuTect2 -R $hg38 -I:tumor ${filepath}recal_bam/${pfx}.sorted.recal.bam --dbsnp $dbsnp --maxReadsInRegionPerSample 10000 -A Coverage -A AlleleBalance --max_alternate_alleles 40 -bamout ${filepath}hotspot_bam/${pfx}.out.bam -o ${filepath}hotspot_vcf/${pfx}_hotspot.vcf -L $hotspots
done <${filepath}samples_noLane.txt
