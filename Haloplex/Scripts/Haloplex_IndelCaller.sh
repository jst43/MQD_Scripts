#!/bin/bash

#CONSTANTS
hg38="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.fa"
hg38Dict="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.dict"
TABLE_ANNOVAR="/media/eguz/darwin/Resources/Software/annovar/table_annovar.pl"
humandb="/media/eguz/darwin/Resources/Software/annovar/humandb"
PICARD="/media/eguz/darwin/Resources/Software/picard-tools-1.141/picard.jar"

#COMMANDLINE VARIABLES
while getopts "fh" opt; do
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
echo "######################################################"
echo "#                Haloplex_IndelsCaller               #"
echo "#      Writen by Joe Thompson (jst43@cam.ac.uk)      #"
echo "#                                                    #"
echo "#                  July 11th 2017                    #"
echo "# This Bash script uses the following software under #"
echo "#         GNU Public license v2: Pindel, vcftools    #"
echo "#                Picard and Annovar.                 #"
echo "#                                                    #"
echo "######################################################"
echo ""
echo "Filepath is $filepath"
echo ""

cd $filepath

#BEFORE TO START: In the folder should be included sorted.bam files and header.txt
#Create list of name samples
ls ./recal_bam/*.sorted.realigned.recal.bam > samplesPindel.txt
sed -i 's|.sorted.realigned.recal.bam||' samplesPindel.txt
#Create files: pindelinput.txt
for k in `cat samplesPindel.txt`; do
	ls ./recal_bam/${k}.sorted.realigned.recal.bam > ${k}_pindelinput.txt
	sed -i 's|\(.\+\)\.sorted.realigned.recal.bam|\1\.sorted.realigned.recal.bam\t350\t\1|' ${k}_pindelinput.txt
done

#Run Pindel
for i in `cat samplesPindel.txt`; do
	pindel -f ${hg38} -i ${i}_pindelinput.txt -c ALL -o ${i}
#Create vcf files of deletion and insertion out of pindel
	pindel2vcf -p ${i}_D -r ${hg38} -R GRCh38 -d 201312 -G -v ${i}_D.vcf
	pindel2vcf -p ${i}_SI -r ${hg38} -R GRCh38 -d 201312 -G -v ${i}_SI.vcf
#Merging both duplicates and filtering vcf
        bgzip -c ${i}_D.vcf > ${i}_D.vcf.gz
        tabix -p vcf ${i}_D.vcf.gz
        bgzip -c ${i}_SI.vcf > ${i}_SI.vcf.gz
        tabix -p vcf ${i}_SI.vcf.gz
        vcftools_0.1.13 --vcf ${i}_SI.vcf --max-missing 1 --recode --out ${i}_SI_F1
#Merging Deletion and Insertion
	vcf-concat ${i}_D.vcf.gz ${i}_SI.vcf.gz > ${i}_D+SI.vcf
	java -jar ${PICARD} SortVcf I=${i}_D+SI.vcf O=${i}_D+SI.sorted.vcf SD=${hg38Dict}
	less ${i}_D+SI.sorted.vcf | grep -v "#" > ${i}_headless.vcf
        cat header.txt ${i}_headless.vcf > ${i}_NewHead.vcf
#Annotation hg38
	perl ${TABLE_ANNOVAR} ${i}_NewHead.vcf ${humandb} -buildver hg38 -out ${i}_D+SI.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
done

echo "lympSeq_Indels analysis was completed."
