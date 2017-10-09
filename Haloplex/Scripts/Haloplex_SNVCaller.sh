#!/bin/bash

#CONSTANTS
hg38="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.fa"
GATK="/media/eguz/darwin/Resources/Software/GenomeAnalysisTK.jar"
All="/media/eguz/darwin/Resources/hg38.p6/All.vcf"
TABLE_ANNOVAR="/media/eguz/darwin/Resources/Software/annovar/table_annovar.pl"
humandb="/media/eguz/darwin/Resources/Software/annovar/humandb"
PICARD="/media/eguz/darwin/Resources/Software/picard-tools-1.141/picard.jar"
java="/media/eguz/darwin/Resources/Software/jre1.8.0_112/bin/java"
LocatIt="/media/eguz/darwin/Resources/Software/LocatIt_v3.5.1.46.jar"
dedup_bed="/home/joe/MQD_Scripts/Haloplex/Dependent_Files/dedup_amplicons38.bed"
coverage_bed="/home/joe/MQD_Scripts/Haloplex/Dependent_Files/coverage_amplicons38.bed"

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
echo "#                Haloplex SNV Caller                 #"
echo "#    Writen by Joe Thompson (jst43@cam.ac.uk)        #"
echo "#                                                    #"
echo "#                  July 11th 2017                    #"
echo "# This Bash script uses the following software under #"
echo "#    GNU Public license v2: bwa, GATK, vcftools,     #"
echo "#                    and Annovar                     #"
echo "#                                                    #"
echo "######################################################"
echo ""
echo "Filepath is $filepath"
echo ""

cd $filepath

#Make directories
mkdir tempfiles
mkdir recal_bam
mkdir coverage
mkdir vcf_anno

#Change names outputs trimming software!!!!
for k in `cat samples.txt`; do
	#Alignment 1000Genomes(Hg38)
	bwa_0.7.12 mem -R "@RG\tID:<${k}>\tLB:LIBRARY_NAME\tSM:<${k}>\tPL:ILLUMINA" ${hg38} ./trimmed_fastq/${k}_L001_R1_001.trimmed.fastq.gz ./trimmed_fastq/${k}_L001_R2_001.trimmed.fastq.gz > ./tempfiles/${k}.sam
	#Remove Duplicates with LocatIt
	$java -Xmx250g -jar $LocatIt -X $filepath -U -IS -OB -b $dedup_bed -o ./tempfiles/${k}_RMD ./tempfiles/${k}.sam ${k}_L001_I2_001.fastq.gz
	#Convert bam without duplicates in fastq file
	java -Xmx250g -jar ${PICARD} SamToFastq I=./tempfiles/${k}_RMD.bam F=./tempfiles/${k}_R1.fastq F2=./tempfiles/${k}_R2.fastq
	#Realign with bwa mem
	bwa_0.7.12 mem -R "@RG\tID:<${k}>\tLB:LIBRARY_NAME\tSM:<${k}>\tPL:ILLUMINA" ${hg38} ./tempfiles/${k}_R1.fastq ./tempfiles/${k}_R2.fastq > ./tempfiles/${k}_RMD.sam
	#Create bam file, sort + index
	samtools_1.2 view -bS ./tempfiles/${k}_RMD.sam > ./tempfiles/${k}.bam
	samtools_1.2 sort ./tempfiles/${k}.bam ./tempfiles/${k}.sorted
	samtools_1.2 index ./tempfiles/${k}.sorted.bam
	#Realigned and Indels
	java -Xmx250g -jar ${GATK} -nt 20 -T RealignerTargetCreator -R ${hg38} -I ./tempfiles/${k}.sorted.bam -o ./tempfiles/${k}.bam.list
	java -Xmx250g -jar ${GATK} -T IndelRealigner -R ${hg38} -I ./tempfiles/${k}.sorted.bam -targetIntervals ./tempfiles/${k}.bam.list -o ./tempfiles/${k}.sorted.realigned.bam
	#Recalibrator and quality control
	java -Xmx250g -jar ${GATK} -nct 20 -T BaseRecalibrator -R ${hg38} -I ./tempfiles/${k}.sorted.realigned.bam -l info -knownSites ${All} -o ./tempfiles/${k}.sorted.realigned.table
	java -Xmx250g -jar ${GATK} -nct 20 -T PrintReads -R ${hg38} -I ./tempfiles/${k}.sorted.realigned.bam -l INFO -BQSR ./tempfiles/${k}.sorted.realigned.table -o ./recal_bam/${k}.sorted.realigned.recal.bam
	java -Xmx250g -jar ${GATK} -T DepthOfCoverage -R ${hg38} -o ./coverage/${k}.coverage -I ./recal_bam/${k}.sorted.realigned.recal.bam -L $coverage_bed
	#Calling variants
	java -Xmx250g -jar ${GATK} -T UnifiedGenotyper -R ${hg38} -I ./recal_bam/${k}.sorted.realigned.recal.bam -glm BOTH --dbsnp ${All} -stand_call_conf 30.0 -stand_emit_conf 10.0 -A Coverage -dcov 10000 -A AlleleBalance --max_alternate_alleles 40 -o ./tempfiles/${k}.vcf
	#Filtering	
	vcftools_0.1.13 --vcf ./tempfiles/${k}.vcf --minQ 30 --recode --out ./tempfiles/${k}_F1
	vcftools_0.1.13 --vcf ./tempfiles/${k}_F1.recode.vcf --min-meanDP 20 --recode --out ./tempfiles/${k}_F2
	#Annotation
	perl ${TABLE_ANNOVAR} ./tempfiles/${k}_F2.recode.vcf ${humandb} -buildver hg38 -out ./vcf_anno/${k}_SNVs.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f -nastring . -vcfinput
done
