#!/bin/bash

#COMMANDLINE VARIABLES
while getopts "f:th" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		t)
			trimflag='-t' >&2
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

./Haloplex_Trim.sh -f $filepath

if [ $trimflag == '-t' ]; then
	./Haloplex_Trim2.py -f $filepath
fi

./Haloplex_AlignProcess.sh -f $filepath $trimflag

./Haloplex_HotspotProcess.sh -f $filepath

./Haloplex_SNVCaller.sh -f $filepath & ./Haloplex_IndelCaller.sh -f $filepath & ./Haloplex_HotspotCaller.sh -f $filepath

./Haloplex_Coverage.sh -f $filepath & ./Haloplex_CNVCaller.R $filepath

./Haloplex_Annotate.sh -f $filepath

./Haloplex_MergeTSV.R $filepath

./Haloplex_FilterTSV.R $filepath

./Haloplex_Clean.sh -f $filepath

