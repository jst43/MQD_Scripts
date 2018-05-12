#!/bin/bash

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
	filepath=`pwd`"/"
fi

if [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

#SCRIPT

mv ${filepath}hotspot/bam ${filepath}Output/
mv ${filepath}Output/bam ${filepath}Output/hotspot_bam
mv ${filepath}realigned_recal_bam ${filepath}Output/

mkdir ${filepath}generated_data

dataname=`basename $filepath`

for i in `ls ${filepath}`; do
	if [ $i != 'Output' ]; then
		if [ $i != 'generated_data' ]; then
			mv $i ${filepath}generated_data/
		fi
	fi
done

tar cvfj ${filepath}${dataname}.bz2 ${filepath}generated_data/

rm -r ${filepath}generated_data

