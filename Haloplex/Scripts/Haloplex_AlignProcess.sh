#!/bin/bash

#CONSTANTS
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
GATKv3_5="../../../Software/GenomeAnalysisTK_v3.5.jar"
All="../../../hg38.p6/All.vcf"
PICARD="../../../Software/picard-tools-1.141/picard.jar"
java="../../../Software/jre1.8.0_112/bin/java"
LocatIt="../../../Software/LocatIt_v4.0.1.jar"
dedup_bed="../Dependent_Files/dedup_amplicons38.bed"
coverage_bed="../Dependent_Files/coverage_amplicons38.bed"

#COMMANDLINE VARIABLES
while getopts "f:th" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		t)
			trimflag='2' >&2
			echo "Trim flag has been set"
			;;
		h)
			echo "Usage: $0 [-f FILEPATH (optional)] " >&2
			echo
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
			echo "  -t		set this flag to run $0 on the output of Haloplex_Trim2.py"
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

#SCRIPT
echo "Filepath is $filepath"

if [ ! -d ${filepath}trimmed_fastq ]; then
        echo "Directory trimmed_fastq not found"
        exit 1
fi

#Create samples file without Lane information
sed 's|_L[0-9]*||' ${filepath}samples.txt > ${filepath}samples_noLane.txt

R1_name=`grep "^R1" ${filepath}fastq_types.txt | sed -e 's|R1 type:||'`
R2_name=`grep "^R2" ${filepath}fastq_types.txt | sed -e 's|R2 type:||'`
Index_name=`grep "Index" ${filepath}fastq_types.txt | sed -e 's|Index type:||'`

#Make directories
mkdir ${filepath}postprocess_files
mkdir ${filepath}dedup_data
mkdir ${filepath}realigned_recal_bam

while read lane <&3 && read nolane <&4; do
	#Alignment 1000Genomes(Hg38)
	bwa mem -R "@RG\tID:<${nolane}>\tLB:LIBRARY_NAME\tSM:<${nolane}>\tPL:ILLUMINA" $hg38 ${filepath}trimmed_fastq/${lane}_${R1_name}.trimmed${trimflag}.fastq.gz ${filepath}trimmed_fastq/${lane}_${R2_name}.trimmed${trimflag}.fastq.gz > ${filepath}postprocess_files/${nolane}.sam
	#Remove Duplicates with LocatIt
	$java -Xmx40g -jar $LocatIt -X $filepath -U -IS -OB -b $dedup_bed -o ${filepath}postprocess_files/${nolane}_Dedup ${filepath}postprocess_files/${nolane}.sam ${filepath}${lane}_${Index_name}.fastq.gz
	mv ${filepath}postprocess_files/${nolane}_Dedup.properties ${filepath}dedup_data/
	#Convert bam without duplicates in fastq file
	$java -Xmx40g -jar $PICARD SamToFastq I=${filepath}postprocess_files/${nolane}_Dedup.bam F=${filepath}postprocess_files/${nolane}_${R1_name}.fastq F2=${filepath}postprocess_files/${nolane}_${R2_name}.fastq
	#Realign with bwa mem
	bwa mem -R "@RG\tID:<${nolane}>\tLB:LIBRARY_NAME\tSM:<${nolane}>\tPL:ILLUMINA" $hg38 ${filepath}postprocess_files/${nolane}_${R1_name}.fastq ${filepath}postprocess_files/${nolane}_${R2_name}.fastq > ${filepath}postprocess_files/${nolane}_Dedup.sam
	#Create bam file, sort + index
	samtools view -bS ${filepath}postprocess_files/${nolane}_Dedup.sam > ${filepath}postprocess_files/${nolane}.bam
	samtools sort ${filepath}postprocess_files/${nolane}.bam -o ${filepath}postprocess_files/${nolane}.sorted.bam
	samtools index ${filepath}postprocess_files/${nolane}.sorted.bam
	#Realigned and Indels
	$java -Xmx40g -jar $GATKv3_5 -nt 20 -T RealignerTargetCreator -R $hg38 -I ${filepath}postprocess_files/${nolane}.sorted.bam -o ${filepath}postprocess_files/${nolane}.bam.list
	$java -Xmx40g -jar $GATKv3_5 -T IndelRealigner -R $hg38 -I ${filepath}postprocess_files/${nolane}.sorted.bam -targetIntervals ${filepath}postprocess_files/${nolane}.bam.list -o ${filepath}postprocess_files/${nolane}.sorted.realigned.bam
	#Recalibrator and quality control
	$java -Xmx40g -jar $GATKv3_5 -nct 20 -T BaseRecalibrator -R $hg38 -I ${filepath}postprocess_files/${nolane}.sorted.realigned.bam -l info -knownSites $All -o ${filepath}postprocess_files/${nolane}.sorted.realigned.table
	$java -Xmx40g -jar $GATKv3_5 -nct 20 -T PrintReads -R $hg38 -I ${filepath}postprocess_files/${nolane}.sorted.realigned.bam -l INFO -BQSR ${filepath}postprocess_files/${nolane}.sorted.realigned.table -o ${filepath}realigned_recal_bam/${nolane}.sorted.realigned.recal.bam
done 3<${filepath}samples.txt 4<${filepath}samples_noLane.txt
