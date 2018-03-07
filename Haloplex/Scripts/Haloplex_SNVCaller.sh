#!/bin/bash

#CONSTANTS
hg38="../../../hg38.p6/hg38_2MergeAll.fa"
GATK="../../..//Software/GenomeAnalysisTK.jar"
All="../../../hg38.p6/All.vcf"
TABLE_ANNOVAR="../../../Software/annovar/table_annovar.pl"
humandb="../../../Software/annovar/humandb"
PICARD="../../../Software/picard-tools-1.141/picard.jar"
java="../../../Software/jre1.8.0_112/bin/java"
LocatIt="../../../Software/LocatIt_v3.5.1.46.jar"
dedup_bed="../Dependent_Files/dedup_amplicons38.bed"
coverage_bed="../Dependent_Files/coverage_amplicons38.bed"
blacklist_bed="../Dependent_Files/blacklist.bed"

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

if [ ! -d trimmed_fastq ]; then
        echo "Directory trimmed_fastq not found"
        exit 1
fi

#Create samples file without Lane information
sed 's|_L[0-9]*||' ${filepath}samples.txt > ${filepath}samples_noLane.txt

R1_name=`grep "^R1" ${filepath}fastq_types.txt | sed -e 's|R1 type:||'`
R2_name=`grep "^R2" ${filepath}fastq_types.txt | sed -e 's|R2 type:||'`
Index_name=`grep "Index" ${filepath}fastq_types.txt | sed -e 's|Index type:||'`

#Make directories
mkdir tempfiles
mkdir recal_bam
mkdir coverage
mkdir vcf_anno

#Change names outputs trimming software!!!!
while read lane <&3 && read nolane <&4; do
	#Alignment 1000Genomes(Hg38)
	bwa mem -R "@RG\tID:<${nolane}>\tLB:LIBRARY_NAME\tSM:<${nolane}>\tPL:ILLUMINA" $hg38 trimmed_fastq/${lane}_${R1_name}.trimmed.fastq.gz trimmed_fastq/${lane}_${R2_name}.trimmed.fastq.gz > tempfiles/${nolane}.sam
	#Remove Duplicates with LocatIt
	$java -Xmx40g -jar $LocatIt -X $filepath -U -IS -OB -b $dedup_bed -o tempfiles/${nolane}_RMD tempfiles/${nolane}.sam ${lane}_${Index_name}.fastq.gz
	#Convert bam without duplicates in fastq file
	$java -Xmx40g -jar $PICARD SamToFastq I=tempfiles/${nolane}_RMD.bam F=tempfiles/${nolane}_${R1_name}.fastq F2=tempfiles/${nolane}_${R2_name}.fastq
	#Realign with bwa mem
	bwa mem -R "@RG\tID:<${nolane}>\tLB:LIBRARY_NAME\tSM:<${nolane}>\tPL:ILLUMINA" $hg38 tempfiles/${nolane}_${R1_name}.fastq tempfiles/${nolane}_${R2_name}.fastq > tempfiles/${nolane}_RMD.sam
	#Create bam file, sort + index
	samtools view -bS tempfiles/${nolane}_RMD.sam > tempfiles/${nolane}.bam
	samtools sort tempfiles/${nolane}.bam -o tempfiles/${nolane}.sorted.bam
	samtools index tempfiles/${nolane}.sorted.bam
	#Recalibrator and quality control
	$java -Xmx40g -jar $GATK -nct 20 -T BaseRecalibrator -R $hg38 -I tempfiles/${nolane}.sorted.bam -l info -knownSites $All -o tempfiles/${nolane}.sorted.table
	$java -Xmx40g -jar $GATK -nct 20 -T PrintReads -R $hg38 -I tempfiles/${nolane}.sorted.bam -l INFO -BQSR tempfiles/${nolane}.sorted.table -o recal_bam/${nolane}.sorted.recal.bam
	$java -Xmx40g -jar $GATK -T DepthOfCoverage -R $hg38 -o coverage/${nolane}.coverage -I recal_bam/${nolane}.sorted.recal.bam -L $coverage_bed
	#Calling variants
	$java -Xmx40g -jar $GATK -T HaplotypeCaller -R $hg38 -I recal_bam/${nolane}.sorted.recal.bam -glm BOTH --dbsnp $All -stand_call_conf 30.0 -stand_emit_conf 10.0 -A Coverage -dcov 10000 -A AlleleBalance --max_alternate_alleles 40 -o tempfiles/${nolane}.vcf
	#Filtering	
	vcftools --vcf tempfiles/${nolane}.vcf --minQ 30 --recode --out tempfiles/${nolane}_F1
	vcftools --vcf tempfiles/${nolane}_F1.recode.vcf --min-meanDP 20 --recode --out tempfiles/${nolane}_F2
	bedtools intersect -v -a tempfiles/${nolane}_F2.recode.vcf -b $blacklist_bed -header > tempfiles/${nolane}_bedfiltered.vcf
	#Annotation
	perl $TABLE_ANNOVAR tempfiles/${nolane}_bedfiltered.vcf $humandb -buildver hg38 -out vcf_anno/${nolane}_SNVs.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f -nastring . -vcfinput
done 3<${filepath}samples.txt 4<${filepath}samples_noLane.txt
