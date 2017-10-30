#!/bin/bash
hg38="/mnt/raid/Resources/hg38.p6/hg38_2MergeAll.fa"
hg38Dict="/mnt/raid/Resources/hg38.p6/hg38_2MergeAll.dict"
TABLE_ANNOVAR="/mnt/raid/Resources/Software/annovar/table_annovar.pl"
humandb="/mnt/raid/Resources/Software/annovar/humandb"
PICARD="/mnt/raid/Resources/Software/picard-tools-1.141/picard.jar"
header="/mnt/raid/Resources/MQD_Scripts/Fluidigm/Dependent_Files/header.txt"

while getopts "fh" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		h)
			echo "Usage: $0 [-f FILEPATH (optional)]" >&2
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
echo "#               Fluidigm Indel Caller                #"
echo "#    Writen by Joe Thompson (jst43@cam.ac.uk)        #"
echo "#                                                    #"
echo "#                  July 11th 2017                    #"
echo "# This Bash script uses the following software under #"
echo "#    GNU Public license v2: fastqc, pindel,          #"
echo "#               vcftools and Annovar.                #"
echo "#                                                    #"
echo "######################################################"
echo ""
echo "Filepath is $filepath"
echo ""

cd $filepath

#Make Directories
mkdir pindeltemp
mkdir pindel_anno

cd pindeltemp

#Create list of name samples
ls ../recal_bam/*.recal.bam > samplesPindel.txt
sed -i 's|../recal_bam/||' samplesPindel.txt
sed -i 's|.sorted.realigned.recal.bam||' samplesPindel.txt

#Create files: pindelinput.txt
for k in `cat samplesPindel.txt`; do
        ls ../recal_bam/${k}.sorted.realigned.recal.bam > ${k}_pindelinput.txt
	sed -i 's|\(.\+\)\.sorted.realigned.recal.bam|\1\.sorted.realigned.recal.bam\t350\t\1|' ${k}_pindelinput.txt
done

#Run Pindel
for i in `cat samplesPindel.txt`; do
        pindel -f ${hg38} -i ${i}_pindelinput.txt -c ALL -o ${i}
#Create vcf files of deletion and insertion out of pindel
        pindel2vcf -p ${i}_D -r ${hg38} -R GRCh38 -d 201312 -G -v ${i}_D.vcf
        pindel2vcf -p ${i}_SI -r ${hg38} -R GRCh38 -d 201312 -G -v ${i}_SI.vcf
done

echo "vcf files for deletions and insertions were made."
#create list of name sample no-duplicates
ls *-Dup1_D.vcf > names.txt
sed -i 's|_D.vcf||' names.txt
sed -i 's|-Dup1||' names.txt

for i in `cat names.txt`; do
#Merging both duplicates and filtering vcf
        bgzip -c ${i}-Dup1_D.vcf > ${i}-Dup1_D.vcf.gz
        tabix -p vcf ${i}-Dup1_D.vcf.gz
        bgzip -c ${i}-Dup2_D.vcf > ${i}-Dup2_D.vcf.gz
        tabix -p vcf ${i}-Dup2_D.vcf.gz
        vcf-merge ${i}-Dup1_D.vcf.gz ${i}-Dup2_D.vcf.gz > ${i}_D.vcf
        vcftools_0.1.13 --vcf ${i}_D.vcf --max-missing 1 --recode --out ${i}_D_F1
        bgzip -c ${i}-Dup1_SI.vcf > ${i}-Dup1_SI.vcf.gz
        tabix -p vcf ${i}-Dup1_SI.vcf.gz
        bgzip -c ${i}-Dup2_SI.vcf > ${i}-Dup2_SI.vcf.gz
        tabix -p vcf ${i}-Dup2_SI.vcf.gz
        vcf-merge ${i}-Dup1_SI.vcf.gz ${i}-Dup2_SI.vcf.gz > ${i}_SI.vcf
        vcftools_0.1.13 --vcf ${i}_SI.vcf --max-missing 1 --recode --out ${i}_SI_F1
        mv ${i}_D_F1.recode.vcf ${i}_D_F1.vcf
        mv ${i}_SI_F1.recode.vcf ${i}_SI_F1.vcf
#Merging Deletion and Insertion
        bgzip -c ${i}_D_F1.vcf > ${i}_D.vcf.gz
        bgzip -c ${i}_SI_F1.vcf > ${i}_SI.vcf.gz
        tabix -p vcf ${i}_D.vcf.gz
        tabix -p vcf ${i}_SI.vcf.gz
        vcf-concat ${i}_D.vcf.gz ${i}_SI.vcf.gz > ${i}_DSI.vcf
        java -jar ${PICARD} SortVcf I=${i}_DSI.vcf O=${i}_DSI.sorted.vcf SD=${hg38Dict}
        less ${i}_DSI.sorted.vcf | grep -v "#" > ${i}_headless.vcf
        cat $header ${i}_headless.vcf > ${i}_NewHead.vcf
#Annotation hg38
        perl ${TABLE_ANNOVAR} ${i}_NewHead.vcf ${humandb} -buildver hg38 -out ../pindel_anno/${i}_DSI.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
done
echo "lympSeq_Indels analysis was completed."

