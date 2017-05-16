#!/bin/bash

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

while getopts ":p:s:f:h" opt; do
	case $opt in
		p)
			Pipeline=$OPTARG >&2
			echo "Pipeline is $OPTARG" >&2
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
			echo "Usage: $0 [-f] FILEPATH" >&2
			echo
			echo "	-p		pipeline argument; choice of SNV or Pindel"
			echo "	-s		sample prefix, eg. MP, DLBCL"
			echo "	-f		filepath to directory containing fastq.gz files"
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

if [ ${Pipeline} == "SNV" ];
then
	FileString="_SNVs.myanno.hg38_multianno.txt"
elif [ ${Pipeline} == "Pindel" ];
then
	FileString="_DSI.myanno.hg38_multianno.txt"
fi

cd $Filepath

ls ${Sample_Prefix}*${FileString} > ${Sample_Prefix}_libraryHALOnames_${Pipeline}.txt
for i in `cat ${Sample_Prefix}_libraryHALOnames_${Pipeline}.txt`; do
	sed -i 's|^|'"${i}"'\t|' ${i}
done

cat ${Sample_Prefix}*${FileString} > ${Sample_Prefix}_LibraryHALO_${Pipeline}_hg38_multianno.txt
sed -i 's|'"${FileString}"'||' ${Sample_Prefix}_LibraryHALO_${Pipeline}_hg38_multianno.txt
