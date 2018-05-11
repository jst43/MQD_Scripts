#!/bin/bash

#CONSTANTS
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
hg38Dict="../../../hg38.p6/hg38_2MergeAll.dict"
TABLE_ANNOVAR="../../../Software/annovar/table_annovar.pl"
humandb="../../../Software/annovar/humandb"
PICARD="../../../Software/picard-tools-1.141/picard.jar"
header="../header.txt"
blacklist_bed="../Dependent_Files/blacklist.bed"

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
echo "Filepath is $filepath"
echo ""

if [ ! -d ${filepath}realigned_recal_bam ]; then
        echo "Directory realigned_recal_bam not found"
        exit 1
fi

#Make directories
mkdir ${filepath}pindel
mkdir ${filepath}pindel/tempfiles

#Create list of name samples
(cd ${filepath}realigned_recal_bam/ && ls *.recal.bam > ../samplesPindel.txt)
sed -i 's|.sorted.realigned.recal.bam||' ${filepath}samplesPindel.txt

while read k; do
	ls ${filepath}realigned_recal_bam/${k}.sorted.realigned.recal.bam > ${filepath}pindel/tempfiles/${k}_pindelinput.txt
	sed -i 's|\(.\+\)\.sorted.realigned.recal.bam|\1\.sorted.realigned.recal.bam\t350\t\1|' ${filepath}pindel/tempfiles/${k}_pindelinput.txt
	#Run Pindel
	pindel -f $hg38 -i ${filepath}pindel/tempfiles/${k}_pindelinput.txt -c ALL -o ${filepath}pindel/tempfiles/${k}
	#Create vcf files of deletion and insertion out of pindel
	pindel2vcf -p ${filepath}pindel/tempfiles/${k}_D -r $hg38 -R GRCh38 -d 201312 -G -v ${filepath}pindel/tempfiles/${k}_D.vcf
	pindel2vcf -p ${filepath}pindel/tempfiles/${k}_SI -r $hg38 -R GRCh38 -d 201312 -G -v ${filepath}pindel/tempfiles/${k}_SI.vcf
	bgzip -c ${filepath}pindel/tempfiles/${k}_D.vcf > ${filepath}pindel/tempfiles/${k}_D.vcf.gz
	tabix -p vcf ${filepath}pindel/tempfiles/${k}_D.vcf.gz
	bgzip -c ${filepath}pindel/tempfiles/${k}_SI.vcf > ${filepath}pindel/tempfiles/${k}_SI.vcf.gz
	tabix -p vcf ${filepath}pindel/tempfiles/${k}_SI.vcf.gz
	#Merging Deletion and Insertion
	vcf-concat ${filepath}pindel/tempfiles/${k}_D.vcf.gz ${filepath}pindel/tempfiles/${k}_SI.vcf.gz > ${filepath}pindel/tempfiles/${k}_DSI.vcf
	bedtools intersect -v -a ${filepath}pindel/tempfiles/${k}_DSI.vcf -b $blacklist_bed -header > ${filepath}pindel/tempfiles/${k}_DSI.bedfiltered.vcf
	java -jar $PICARD SortVcf I=${filepath}pindel/tempfiles/${k}_DSI.bedfiltered.vcf O=${filepath}pindel/tempfiles/${k}_DSI.bedfiltered.sorted.vcf SD=${hg38Dict}
done <${filepath}samplesPindel.txt
