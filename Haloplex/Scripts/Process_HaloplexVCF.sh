#!/bin/bash

#CONSTANTS
VCFHEAD="/mnt/raid/Resources/MQD_Scripts/Haloplex/Dependent_Files/Haloplex_VCFHead.txt"
FilterVCF="/mnt/raid/Resources/MQD_Scripts/Haloplex/Scripts/Filter_HaloplexVCF.R"

#COMMANDLINE VARIABLES

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

while getopts ":f:h" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		h)
			echo "Usage: $0 [-f FILEPATH ] " >&2
			echo
			echo "	-f		filepath of VCF to be processed"
			echo "	-h		display this help message"
			exit 1
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument" >&2
			exit 1
			;;
	esac
done

if [ -z $filepath ]; then
	echo "This script requires a filepath to a VCF"
	exit 1
fi

if [ ! -e $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

#Make check that file is a VCF file

#SCRIPT

tail -n +2 $filepath > body.txt
cat $VCFHEAD body.txt > $filepath
rm body.txt

Rscript $FilterVCF $filepath

sed -i 's|GERP.._RS|GERP++_RS|' $filepath
