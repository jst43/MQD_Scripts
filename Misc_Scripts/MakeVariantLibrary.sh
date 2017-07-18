#!/bin/bash

#COMMANDLINE VARIABLES
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

while getopts ":p:s:fh" opt; do
	case $opt in
		p)
			pipeline=$OPTARG >&2
			if [ ${pipeline} == "SNV" ];
			then
			        filestring="_SNVs.myanno.hg38_multianno.txt"
			elif [ ${pipeline} == "Pindel" ];
			then
			        filestring="_DSI.myanno.hg38_multianno.txt"
			fi
			;;
		s)
			sample_prefix=$OPTARG >&2
			;;
		f)
			filepath=$OPTARG >&2
			;;
		h)
			echo "Usage: $0 [ -p PIPELINE ] [ -s PREFIX ] [ -f FILEPATH ]" >&2
			echo
			echo "	-p		pipeline argument; choice of SNV or Pindel"
			echo "	-s		sample prefix, eg. MP, DLBCL"
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
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
	filepath=`pwd`
fi

if [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

#SCRIPT
echo "Pipeline is $pipeline"
echo ""
echo "Sample Prefix is $sample_prefix"
echo ""
echo "Filepath is $filepath"

cd $filepath

ls ${sample_prefix}*${filestring} > ${sample_prefix}_librarynames_${pipeline}.txt
for i in `cat ${sample_prefix}_librarynames_${pipeline}.txt`; do
	sed -i 's|^|'"${i}"'\t|' $i
done

cat ${sample_prefix}*${filestring} > ${sample_prefix}_Library_${pipeline}_hg38_multianno.txt
sed -i 's|'"${filestring}"'||' ${sample_prefix}_Library_${pipeline}_hg38_multianno.txt

if [ -z $sample_prefix ];
	then
		mv ${sample_prefix}_librarynames_${pipeline}.txt All_librarynames_${pipeline}.txt
		mv ${sample_prefix}_Library_${pipeline}_hg38_multianno.txt All_Library_${pipeline}_hg38_multianno.txt
fi
