#!/bin/bash

GATK_v3_8="/mnt/raid/Resources/Software/GenomeAnalysisTK.jar"
GATK_v3_5="/home/joe/GenomeAnalysisTK.jar"
hg38="/mnt/raid/Resources/hg38.p6/hg38_2MergeAll.fa"
All="/mnt/raid/Resources/hg38.p6/All.vcf"
dbsnp="/mnt/raid/Resources/hg38.p6/All.vcf"
TABLE_ANNOVAR="/mnt/raid/Resources/Software/annovar/table_annovar.pl"
humandb="/mnt/raid/Resources/Software/annovar/humandb/"

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
	java -Xmx40g -jar $GATKv3_5 -nct 20 -T BaseRecalibrator -R $hg38 -I ${filepath}tempfiles/${nolane}.sorted.bam -l info -knownSites $All -o ${filepath}tempfiles/${nolane}.sorted.table
	java -Xmx40g -jar $GATKv3_5 -nct 20 -T PrintReads -R $hg38 -I ${filepath}tempfiles/${nolane}.sorted.bam -l INFO -BQSR ${filepath}tempfiles/${nolane}.sorted.table -o ${filepath}recal_bam/${nolane}.sorted.recal.bam
done <${filepath}samples_noLane.txt
