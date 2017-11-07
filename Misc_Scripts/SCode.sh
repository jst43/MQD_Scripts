#!/bin/bash

#Most of the fastq files that are to be run through the pipeline have filenames with a variable 'S Code' (i.e. AITL_49_S784_L001_R1_001.fastq.gz)
#This S Code is variable as the number appearing after 'S' is unique to each fastq file, making processing of multiple files difficult.
#This script removes the S Code from all fastq/fastq.gz files in a directory

#CONSTANTS
extension="_001.fastq"

#COMMANDLINE VARIABLES
while getopts "fgl:ih" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		g)
			extension+=".gz"
			;;
		l)
			lane="_L00" >&2
			lane+=$OPTARG >&2
			;;
		i)
			haloplex="TRUE" >&2
			;;	
		h)
			echo "Usage: $0 [-f FILEPATH (optional)] " >&2
			echo
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
			echo "	-g		specifies .fastq.gz extension. The default extension is .fastq"
			echo "	-l		sequencing lane number, usually 1-4"
			echo "	-i		performs $0 on the Index fastq files used by Haloplex"
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

r1_extension=${lane}_R1${extension}

#SCRIPT
echo "Filepath is $filepath"
echo ""
echo "File extension is $r1_extension"

ls *${r1_extension} > Scode.txt
sed -i "s|${r1_extension}||" Scode.txt
sed -i 's|_S[0-9][0-9]*||' Scode.txt
for k in `cat Scode.txt`; do
        mv ${k}_S*${r1_extension} ${k}${r1_extension}
        mv ${k}_S*${lane}_R2${extension} ${k}${lane}_R2${extension}
	if [ $haloplex == "TRUE" ]; then
		mv ${k}_S*${lane}_I1${extension} ${k}${lane}_I1${extension}
		mv ${k}_S*${lane}_I2${extension} ${k}${lane}_I2${extension}
	fi
done

rm Scode.txt

echo "S Code removed"
