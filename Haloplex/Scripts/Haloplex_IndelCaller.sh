#!/bin/bash

#CONSTANTS
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
hg38Dict="../../../hg38.p6/hg38_2MergeAll.dict"
TABLE_ANNOVAR="../../../Software/annovar/table_annovar.pl"
humandb="../../../Software/annovar/humandb"
PICARD="../../../Software/picard-tools-1.141/picard.jar"
header="../header.txt"
blacklist_bed="../blacklist.bed"

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
	filepath=`pwd`
fi

if [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

#SCRIPT
echo "Filepath is $filepath"
echo ""

if [ ! -d ${filepath}recal_bam ]; then
        echo "Directory recal_bam not found"
        exit 1
fi

#Make directories
mkdir ${filepath}pindeltemp

#Create list of name samples
(cd ${filepath}realigned_recal_bam/ && ls *.recal.bam > ../samplesPindel.txt)
sed -i 's|.sorted.realigned.recal.bam||' ${filepath}samplesPindel.txt

#Create files: pindelinput.txt
mkdir ${filepath}pindel_anno

while read k; do
	ls ${filepath}realigned_recal_bam/${k}.sorted.realigned.recal.bam > ${filepath}pindeltemp/${k}_pindelinput.txt
	sed -i 's|\(.\+\)\.sorted.realigned.recal.bam|\1\.sorted.realigned.recal.bam\t350\t\1|' ${filepath}pindeltemp/${k}_pindelinput.txt
	#Run Pindel
	pindel -f $hg38 -i ${filepath}pindeltemp/${i}_pindelinput.txt -c ALL -o ${filepath}pindeltemp/${i}
	#Create vcf files of deletion and insertion out of pindel
	pindel2vcf -p ${filepath}pindeltemp/${i}_D -r $hg38 -R GRCh38 -d 201312 -G -v ${filepath}pindeltemp/${i}_D.vcf
	pindel2vcf -p ${filepath}pindeltemp/${i}_SI -r $hg38 -R GRCh38 -d 201312 -G -v ${filepath}pindeltemp/${i}_SI.vcf
	#Merging both duplicates and filtering vcf
        bgzip -c ${filepath}pindeltemp/${i}_D.vcf > ${filepath}pindeltemp/${i}_D.vcf.gz
        tabix -p vcf ${filepath}pindeltemp/${i}_D.vcf.gz
	#vcftools_0.1.13 --vcf ${filepath}pindeltemp/${i}_SI.vcf --max-missing 1 --recode --out ${filepath}pindeltemp/${i}_SI_F1
	mv ${filepath}pindeltemp/${i}_SI_F1.recode.vcf ${filepath}pindeltemp/${i}_SI.vcf
        bgzip -c ${filepath}pindeltemp/${i}_SI.vcf > ${filepath}pindeltemp/${i}_SI.vcf.gz
        tabix -p vcf ${filepath}pindeltemp/${i}_SI.vcf.gz
	#Merging Deletion and Insertion
	vcf-concat ${filepath}pindeltemp/${i}_D.vcf.gz ${filepath}pindeltemp/${i}_SI.vcf.gz > ${filepath}pindeltemp/${i}_DSI.vcf
	#vcftools_0.1.13 --vcf ${filepath}pindeltemp/${i}_DSI.vcf --min-meanDP 20 --recode --out ${filepath}pindeltemp/${k}_DSI_F2
	#bedtools intersect -v -a ${filepath}pindeltemp/${k}_DSI_F2.recode.vcf -b $blacklist_bed -header > ${filepath}pindeltemp/${k}_bedfiltered.vcf
	mv ${filepath}pindeltemp/${i}_bedfiltered.vcf ${filepath}pindeltemp/${i}_DSI.vcf
	java -jar $PICARD SortVcf I=${filepath}pindeltemp/${i}_DSI.vcf O=${filepath}pindeltemp/${i}_DSI.sorted.vcf SD=${hg38Dict}
#	less ${filepath}pindeltemp/${i}_DSI.sorted.vcf | grep -v "#" > ${filepath}pindeltemp/${i}_headless.vcf
#       cat $header ${filepath}pindeltemp/${i}_headless.vcf > ${filepath}pindeltemp/${i}_NewHead.vcf
done <${filepath}samplesPindel.txt
