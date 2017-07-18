#!/bin/bash

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

while getopts ":p:s:fh" opt; do
	case $opt in
		p)
			Pipeline=$OPTARG >&2
			echo "Pipeline is $OPTARG" >&2
			if [ ${Pipeline} == "SNV" ];
			then
			        FileString="_SNVs.myanno.hg38_multianno.txt"
			elif [ ${Pipeline} == "Pindel" ];
			then
			        FileString="_DSI.myanno.hg38_multianno.txt"
			fi
			;;
		s)
			Sample_Prefix=$OPTARG >&2
			echo "Sample Prefix is $OPTARG" >&2
			;;
		f)
			Filepath=$OPTARG >&2
			echo "Filepath is $OPTARG" >&2
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

echo "Filepath is $filepath"

cd $Filepath

ls ${Sample_Prefix}*${FileString} > ${Sample_Prefix}_librarynames_${Pipeline}.txt
for i in `cat ${Sample_Prefix}_librarynames_${Pipeline}.txt`; do
	sed -i 's|^|'"${i}"'\t|' $i
done

cat ${Sample_Prefix}*${FileString} > ${Sample_Prefix}_Library_${Pipeline}_hg38_multianno.txt
sed -i 's|'"${FileString}"'||' ${Sample_Prefix}_Library_${Pipeline}_hg38_multianno.txt

if [ -z $Sample_Prefix ];
	then
		mv ${Sample_Prefix}_librarynames_${Pipeline}.txt All_librarynames_${Pipeline}.txt
		mv ${Sample_Prefix}_Library_${Pipeline}_hg38_multianno.txt All_Library_${Pipeline}_hg38_multianno.txt
fi
