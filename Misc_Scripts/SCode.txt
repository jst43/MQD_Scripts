#!/bin/bash

#Most of the fastq files that are to be run through the pipeline have filenames with a variable 'S Code' (i.e. AITL_49_S784_L001_R1_001.fastq.gz)
#This S Code is variable as the number appearing after 'S' is unique to each fastq file, making processing of multiple files difficult.
#This script removes the S Code from all fastq/fastq.gz files in a directory

#CONSTANTS
extension="_L001_R1_001.fastq"

#COMMANDLINE VARIABLES
while getopts "fgh" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		g)
			extension+=".gz"
			;;
		h)
			echo "Usage: $0 [-f FILEPATH (optional)] " >&2
			echo
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
			echo "	-g		specifies .fastq.gz extension. The default extension is .fastq"
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
echo "File extension is $extension"

ls *${extension} > Scode.txt
sed -i "s|${extension}||" Scode.txt
sed -i 's|_S[0-9][0-9]*||' Scode.txt
for k in `cat Scode.txt`; do
        mv ${k}_S*_L001_R1_001.fastq.gz ${k}_L001_R1_001.fastq.gz
        mv ${k}_S*_L001_R2_001.fastq.gz ${k}_L001_R2_001.fastq.gz
done

rm Scode.txt

echo "S Code removed"
