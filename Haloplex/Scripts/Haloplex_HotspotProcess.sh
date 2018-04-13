#!/bin/bash

java="../../../Software/jre1.8.0_112/bin/java"
GATKv3_5="../../../Software/GenomeAnalysisTK_v3.5.jar"
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
All="../../..//hg38.p6/All.vcf"

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

mkdir recal_bam

while read pfx; do
	$java -Xmx40g -jar $GATKv3_5 -nct 20 -T BaseRecalibrator -R $hg38 -I ${filepath}tempfiles/${nolane}.sorted.bam -l info -knownSites $All -o ${filepath}tempfiles/${nolane}.sorted.table
	$java -Xmx40g -jar $GATKv3_5 -nct 20 -T PrintReads -R $hg38 -I ${filepath}tempfiles/${nolane}.sorted.bam -l INFO -BQSR ${filepath}tempfiles/${nolane}.sorted.table -o ${filepath}recal_bam/${nolane}.sorted.recal.bam
done <${filepath}samples_noLane.txt
