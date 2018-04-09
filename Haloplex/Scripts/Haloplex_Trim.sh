#!/bin/bash

#CONSTANTS
SureCallTRIMMER="../../../Software/SurecallTrimmer_v4.0.1.jar"
java8="../../../Software/jre1.8.0_112/bin/java"

#COMMANDLINE VARIABLES
while getopts "fqh" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		q)
			qc=TRUE >&2
			;;
		h)
			echo "Usage: $0 [-f FILEPATH (optional)] [ -q (optional) ]" >&2
			echo
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
			echo "	-q		performs quality check with fastqc"
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

./Haloplex_Input.py $filepath

while read inputname <&3 && read outputname <&4; do
	mv ${filepath}${inputname} ${filepath}${outputname}
done 3<${filepath}original_names.txt 4<${filepath}replacement_names.txt

R1_name=`grep "^R1" ${filepath}fastq_types.txt | sed -e 's|R1 type:||'`
R2_name=`grep "^R2" ${filepath}fastq_types.txt | sed -e 's|R2 type:||'`

mkdir ${filepath}trimmed_fastq

#Create list of sample names including Dup1-Dup2
(cd $filepath && ls *_${R1_name}.fastq.gz > samples.txt)
sed -i "s|_${R1_name}.fastq.gz||" ${filepath}samples.txt

#Create list of file prefixes to allow renaming of SureCallTRIMMER output files
(cd $filepath && ls *_${R1_name}.fastq.gz > pretrim_names.txt)
(cd $filepath && ls *_${R2_name}.fastq.gz >> pretrim_names.txt)
sed -i 's|.fastq.gz||' ${filepath}pretrim_names.txt

while read sample; do
	#Check sequencing quality
	if [ ! -z $qc ]; then
		fastqc ${filepath}${sample}_${R1_name}.fastq.gz
		fastqc ${filepath}${sample}_${R2_name}.fastq.gz
	fi
	#SureCall processes the read sequences to trim low quality bases from the ends, remove adaptor sequences, and mask enzyme footprints (for HaloPlex).
	$java8 -Xmx40g -jar $SureCallTRIMMER -fq1 ${filepath}${sample}_${R1_name}.fastq.gz -fq2 ${filepath}${sample}_${R2_name}.fastq.gz -hs
done <${filename}samples.txt

#Rename SureCallTRIMMER output files
while read pretrim; do
	mv ${filepath}${pretrim}.*.fastq.gz ${filepath}trimmed_fastq/${pretrim}.trimmed.fastq.gz
done <${filepath}pretrim_names.txt

rm ${filepath}pretrim_names.txt
