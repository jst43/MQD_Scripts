#!/bin/bash

#CONSTANTS
SureCallTRIMMER="/media/eguz/darwin/Resources/Software/SurecallTrimmer_v3.5.1.46.jar"
java8="/media/eguz/darwin/Resources/Software/jre1.8.0_112/bin/java"

#COMMANDLINE VARIABLES
while getopts "fh" opt; do
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
echo "######################################################"
echo "#                   Haloplex_Trim                    #"
echo "#       Writen by Joe Thompson (jst43@cam.ac.uk)     #"
echo "#                                                    #"
echo "#                   November 3rd 2016                #"
echo "# This Bash script uses the following software under #"
echo "#    GNU Public license v2: fastqc                   #"
echo "#                                                    #"
echo "######################################################"
echo ""
echo "Filepath is $filepath"
echo ""

cd $filepath

mkdir trimmed_fastq

#Create list of sample names including Dup1-Dup2
ls *_L001_R1_001.fastq.gz > samples.txt
sed -i 's|_L001_R1_001.fastq.gz||' samples.txt

#Create list of file prefixes to allow renaming of SureCallTRIMMER output files
ls *_R*.fastq.gz > pretrimNames.txt
sed -i 's|.fastq.gz||' pretrimNames.txt

for k in `cat samples.txt`; do
#Check sequencing quality
	fastqc ${k}_L001_R1_001.fastq.gz
	fastqc ${k}_L001_R2_001.fastq.gz
#Trimmed: Amplicon sizes should be 190-640 bp.SureCall processes the read sequences to trim low quality bases from the ends, remove adaptor sequences, and mask enzyme footprints (for HaloPlex).
	${java8} -Xmx250g -jar ${SureCallTRIMMER} -fq1 ${filepath}/${k}_L001_R1_001.fastq.gz -fq2 ${filepath}/${k}_L001_R2_001.fastq.gz -hs
done

#Rename SureCallTRIMMER output files
for k in `cat pretrimNames.txt`; do
	mv ${k}.*.fastq.gz ./trimmed_fastq/${k}.trimmed.fastq.gz;
done
