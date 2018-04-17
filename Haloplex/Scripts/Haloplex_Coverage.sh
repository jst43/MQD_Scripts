#!/bin/bash

#CONSTANTS
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
GATKv3_8="../../..//Software/GenomeAnalysisTK_v3.8.jar"
java="../../../Software/jre1.8.0_112/bin/java"
coverage_bed="../Dependent_Files/coverage_amplicons38.bed"

#COMMANDLINE VARIABLES
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
	filepath=`pwd`
fi

if [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

if [ ! -d ${filepath}realigned_recal_bam ]; then
	echo "Can't find ~/realigned_recal_bam/"
	exit 1
fi

# SCRIPT

mkdir ${filepath}coverage

while read nolane; do
	$java -Xmx40g -jar $GATKv3_8 -T DepthOfCoverage -R $hg38 -o ${filepath}coverage/${nolane}.coverage -I ${filepath}realigned_recal_bam/${nolane}.sorted.realigned.recal.bam -L $coverage_bed
done <${filepath}samples_noLane.txt
