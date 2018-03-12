#!/bin/bash

#COMMANDLINE VARIABLES
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

while getopts ":p:s:fr:h" opt; do
	case $opt in
		p)
			pipeline=$OPTARG >&2
			if [ ${pipeline} == "SNV" ];
			then
			        filestring="_SNVs.myanno.hg38_multianno.txt"
			elif [ ${pipeline} == "Pindel" ];
			then
			        filestring="_DSI.myanno.hg38_multianno.txt"
			elif [ ${pipeline} == "Hotspot" ];
			then
				filestring="_hotspot.myanno.hg38_multianno.txt"
			fi
			;;
		s)
			sample_prefix=$OPTARG >&2
			;;
		f)
			filepath=$OPTARG >&2
			;;
		r)
			rename_prefix=$OPTARG >&2
			;;
		h)
			echo "Usage: $0 [ -p PIPELINE ] [ -s PREFIX ] [ -f FILEPATH ]" >&2
			echo
			echo "	-p		pipeline argument; choice of SNV or Pindel"
			echo "	-s		sample prefix, eg. MP, DLBCL, if none given then file prefix"
			echo "			will be 'All'"
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "			if no filepath is given, $0 will use the current directory"
			echo "	-r		Specifies prefix name if not sample prefix or 'All'"
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

mkdir tempLib

cp ${sample_prefix}*${filestring} ./tempLib/
cd ./tempLib/

ls ${sample_prefix}*${filestring} > ${sample_prefix}_librarynames_${pipeline}.txt
for i in `cat ${sample_prefix}_librarynames_${pipeline}.txt`; do
	sed -i 's|^|'"${i}"'\t|' $i
done

cat `head -n 1 ${sample_prefix}_librarynames_${pipeline}.txt` > ${sample_prefix}_Library_${pipeline}_hg38multianno.txt
for i in `tail -n +2 ${sample_prefix}_librarynames_${pipeline}.txt`; do
	tail -n +2 $i > ${i}.tmp && mv ${i}.tmp $i
	cat $i >> ${sample_prefix}_Library_${pipeline}_hg38multianno.txt
done

sed -i 's|'"${filestring}"'||' ${sample_prefix}_Library_${pipeline}_hg38multianno.txt
sed -i '1!b;s|^.*\tChr|Sample\tChr|' ${sample_prefix}_Library_${pipeline}_hg38multianno.txt

mv ${sample_prefix}_Library_${pipeline}_hg38multianno.txt ../
cd ../
rm -r tempLib/

if [ -z $sample_prefix ]; then
	if [ $rename_prefix ]; then
		mv ${sample_prefix}_Library_${pipeline}_hg38multianno.txt ${rename_prefix}_Library_${pipeline}_hg38multianno.txt
	else
		mv ${sample_prefix}_Library_${pipeline}_hg38multianno.txt All_Library_${pipeline}_hg38multianno.txt
	fi
fi
