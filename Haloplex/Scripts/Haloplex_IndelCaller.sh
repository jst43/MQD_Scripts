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

cd $filepath

if [ ! -d recal_bam ]; then
        echo "Directory recal_bam not found"
        exit 1
fi

#Make directories
mkdir pindeltemp
mkdir pindel_anno

#Create list of name samples
ls ../recal_bam/*.recal.bam > samplesPindel.txt
sed -i 's|../recal_bam/||' samplesPindel.txt
sed -i 's|.sorted.realigned.recal.bam||' samplesPindel.txt
#Create files: pindelinput.txt
for k in `cat samplesPindel.txt`; do
	ls ../recal_bam/${k}.sorted.realigned.recal.bam > ${k}_pindelinput.txt
	sed -i 's|\(.\+\)\.sorted.realigned.recal.bam|\1\.sorted.realigned.recal.bam\t350\t\1|' ${k}_pindelinput.txt
done

for i in `cat samplesPindel.txt`; do
	pindel -f ${hg38} -i ${i}_pindelinput.txt -c ALL -o ${i}
#Create vcf files of deletion and insertion out of pindel
	pindel2vcf -p ${i}_D -r ${hg38} -R GRCh38 -d 201312 -G -v ${i}_D.vcf
	pindel2vcf -p ${i}_SI -r ${hg38} -R GRCh38 -d 201312 -G -v ${i}_SI.vcf
#Merging both duplicates and filtering vcf
        bgzip -c ${i}_D.vcf > ${i}_D.vcf.gz
        tabix -p vcf ${i}_D.vcf.gz
	vcftools_0.1.13 --vcf ${i}_SI.vcf --max-missing 1 --recode --out ${i}_SI_F1
	mv ${i}_SI_F1.recode.vcf ${i}_SI.vcf
        bgzip -c ${i}_SI.vcf > ${i}_SI.vcf.gz
        tabix -p vcf ${i}_SI.vcf.gz
#Merging Deletion and Insertion
	vcf-concat ${i}_D.vcf.gz ${i}_SI.vcf.gz > ${i}_DSI.vcf
	vcftools_0.1.13 --vcf ${i}_DSI.vcf --min-meanDP 20 --recode --out ${k}_DSI_F2
	bedtools intersect -v -a ${k}_DSI_F2.recode.vcf -b $blacklist_bed -header > ${k}_bedfiltered.vcf
	mv ${i}_bedfiltered.vcf ${i}_DSI.vcf
	java -jar ${PICARD} SortVcf I=${i}_DSI.vcf O=${i}_DSI.sorted.vcf SD=${hg38Dict}
#	less ${i}_DSI.sorted.vcf | grep -v "#" > ${i}_headless.vcf
#        cat $header ${i}_headless.vcf > ${i}_NewHead.vcf
done
