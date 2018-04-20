#!/bin/bash

#CONSTANTS
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
GATK="../../..//Software/GenomeAnalysisTK.jar"
All="../../../hg38.p6/All.vcf"
java="../../../Software/jre1.8.0_112/bin/java"

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

#SCRIPT
echo "Filepath is $filepath"
echo ""

if [ ! -d ${filepath}realigned_recal_bam ]; then
        echo "Directory realigned_recal_bam not found"
        exit 1
fi

while read lane <&3 && read nolane <&4; do
	$java -Xmx40g -jar $GATK -T HaplotypeCaller -R $hg38 -I ${filepath}realigned_recal_bam/${nolane}.sorted.realigned.recal.bam -glm BOTH --dbsnp $All -stand_call_conf 30.0 -stand_emit_conf 10.0 -A Coverage -dcov 10000 -A AlleleBalance --max_alternate_alleles 40 -o ${filepath}tempfiles/${nolane}.vcf
done 3<${filepath}samples.txt 4<${filepath}samples_noLane.txt
