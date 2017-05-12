#!/bin/bash

#CONSTANTS
TRIM_GALORE="/media/eguz/darwin/Resources/Software/trim_galore_zip/trim_galore"
hg38="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.fa"
PICARD="/media/eguz/darwin/Resources/Software/picard-tools-1.141/picard.jar"
GATK="/media/eguz/darwin/Resources/Software/GenomeAnalysisTK.jar"
All="/media/eguz/darwin/Resources/hg38.p6/All.vcf"
TABLE_ANNOVAR="/media/eguz/darwin/Resources/Software/annovar/table_annovar.pl"
humandb="/media/eguz/darwin/Resources/Software/annovar/humandb"
fluidigmPrefix="/home/joe/MQD_Scripts/Fluidigm/Dependent_Files/"

#COMMANDLINE VARIABLES
if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

while getopts ":f:g:h" opt; do
	case $opt in
		f)
			filepath=$OPTARG >&2
			;;
		g)
			if [ $OPTARG == 22 ] || [ $OPTARG == 10 ] || [ $OPTARG == 8 ]; then
				genepanel=$OPTARG
				ampliconfile="${fluidigmPrefix}Amplicons/amplicons_${genepanel}Panel.bed"
				primerfile="${fluidigmPrefix}DegeneratePrimers/degeneratePrimers_${genepanel}Panel.txt"
			else
				echo "Unrecognised gene panel; current gene panels are 8, 10, and 22"
				exit 1
			fi
			;;
		h)
			echo "Usage: $0 [-f FILEPATH] [-g GENE_PANEL] " >&2
			echo
			echo "	-f		filepath to directory containing fastq.gz files"
			echo "	-g		gene panel, provides links to the amplicon and"
			echo "			degnerate primers files. Current options are"
			echo "			8, 10, and 22"
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

if [ -z $filepath ] || [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

if [ -z $genepanel ]; then
	echo "This script requires a gene panel; current gene panels are 8, 10, and 22"
	exit 1
fi

#SCRIPT
echo "######################################################"
echo "#                Fluidigm SNV Caller                 #"
echo "#    Writen by Joe Thompson (jst43@cam.ac.uk)        #"
echo "#                                                    #"
echo "#                  April 4th 2017                    #"
echo "# This Bash script uses the following software under #"
echo "#    GNU Public license v2: vim, fastqc, cutadapt,   #"
echo "#      bwa, samtools, GATK, vcftools and Annovar.    #"
echo "#                                                    #"
echo "######################################################"
echo ""
echo "Filepath is $filepath"
echo ""
echo "Gene Panel is $genepanel"

cd $filepath

#Make directories
mkdir trimmed_fastq
mkdir recal_bam
mkdir sorted_bam
mkdir coverage
mkdir vcf_anno
mkdir tempfiles

#Create list of sample names including Dup1-Dup2
ls *_L001_R1_001.fastq.gz > samples.txt
sed -i 's|_L001_R1_001.fastq.gz||' samples.txt

for k in `cat samples.txt`; do
#Check sequencing quality
#       fastqc ${k}_L001_R1_001.fastq.gz
#       fastqc ${k}_L001_R2_001.fastq.gz
#To visualize the quality results
#       firefox ${k}fastqc ${k}_L001_R1_001_fastqc.html &
#       firefox ${k}fastqc ${k}_L001_R2_001_fastqc.html &
#Trimmed
	cutadapt -b TGTAGAACCATGTCGTCAGTGT -b AGACCAAGTCTCTGCTACCGT ${k}_L001_R1_001.fastq.gz | gzip -c > ./tempfiles/${k}.1.fq.gz
	cutadapt -b TGTAGAACCATGTCGTCAGTGT -b AGACCAAGTCTCTGCTACCGT ${k}_L001_R2_001.fastq.gz | gzip -c > ./tempfiles/${k}.2.fq.gz
	perl ${TRIM_GALORE} --paired --length 50 --output_dir ./trimmed_fastq/ ./tempfiles/${k}.1.fq.gz ./tempfiles/${k}.2.fq.gz
#Alignment 1000Genomes(Hg38)
	bwa_0.7.12 mem -R "@RG\tID:<${k}>\tLB:LIBRARY_NAME\tSM:<${k}>\tPL:ILLUMINA" ${hg38} ./trimmed_fastq/${k}.1_val_1.fq.gz ./trimmed_fastq/${k}.2_val_2.fq.gz > ./tempfiles/${k}.sam
#Create bam file and clean + sort + index
	samtools_1.2 view -bS ./tempfiles/${k}.sam > ./tempfiles/${k}.bam
	samtools_1.2 sort ./tempfiles/${k}.bam ./sorted_bam/${k}.sorted
	samtools_1.2 index ./sorted_bam/${k}.sorted.bam
#Realigned and Indels
	java -Xmx250g -jar ${GATK} -nt 20 -T RealignerTargetCreator -R ${hg38} -I ./sorted_bam/${k}.sorted.bam -o ./tempfiles/${k}.bam.list -drf DuplicateRead
	java -Xmx250g -jar ${GATK} -T IndelRealigner -R ${hg38} -I ./sorted_bam/${k}.sorted.bam -targetIntervals ./tempfiles/${k}.bam.list -o ./tempfiles/${k}.sorted.realigned.bam
#Recalibrator and quality control
	java -Xmx250g -jar ${GATK} -nct 20 -T BaseRecalibrator -R ${hg38} -I ./tempfiles/${k}.sorted.realigned.bam -l info -knownSites ${All} -o ./tempfiles/${k}.sorted.realigned.table -drf DuplicateRead
	java -Xmx250g -jar ${GATK} -nct 20 -T PrintReads -R ${hg38} -I ./tempfiles/${k}.sorted.realigned.bam -l INFO -BQSR ./tempfiles/${k}.sorted.realigned.table -o ./recal_bam/${k}.sorted.realigned.recal.bam
	java -Xmx80g -jar ${GATK} -T DepthOfCoverage -R ${hg38} -o ./coverage/${k}.coverage -I ./recal_bam/${k}.sorted.realigned.recal.bam -L $ampliconfile
done

#Coverage
ls ./sorted_bam/*.sorted.bam > input_sortedbams.list
ls ./recal_bam/*.sorted.realigned.recal.bam > input_recalbams.list
java -Xmx80g -jar ${GATK} -T DepthOfCoverage -R ${hg38} -o ./coverage/Library_sortedbam.coverage -I input_sortedbams.list -L $ampliconfile
java -Xmx80g -jar ${GATK} -T DepthOfCoverage -R ${hg38} -o ./coverage/Library_recalbam.coverage -I input_recalbams.list -L $ampliconfile

#Calling variants + Filtering + Annotation
ls ./recal_bam/*-Dup1.sorted.realigned.recal.bam > names.txt
sed -i 's|-Dup1.sorted.realigned.recal.bam||' names.txt
sed -i 's|./recal_bam/||' names.txt

for i in `cat names.txt`; do
        java -Xmx250g -jar ${GATK} -nt 20 -T UnifiedGenotyper -R ${hg38} -I ./recal_bam/${i}-Dup1.sorted.realigned.recal.bam -I ./recal_bam/${i}-Dup2.sorted.realigned.recal.bam -glm BOTH -D ${All} -metrics ./tempfiles/snps.metrics -stand_call_conf 30.0 -stand_emit_conf 10.0 -dcov 1000 -A Coverage -A AlleleBalance --max_alternate_alleles 40 -o ./tempfiles/${i}.vcf -drf DuplicateRead
        vcftools_0.1.13 --vcf ./tempfiles/${i}.vcf --exclude $primerfile --recode --out ./tempfiles/${i}_RDP
        mv ./tempfiles/${i}_RDP.recode.vcf ./tempfiles/${i}.vcf
        vcftools_0.1.13 --vcf ./tempfiles/${i}.vcf --minQ 30 --recode --out ./tempfiles/${i}_F1
        vcftools_0.1.13 --vcf ./tempfiles/${i}_F1.recode.vcf --min-meanDP 50 --recode --out ./tempfiles/${i}_F2
        mv ./tempfiles/${i}_F2.recode.vcf ./tempfiles/${i}_F2.vcf
        vcftools_0.1.13 --vcf ./tempfiles/${i}_F2.vcf --max-missing 1 --recode --out ./tempfiles/${i}_F3
        echo "Duplicates of ${i} were merged"
#Annotation
        perl ${TABLE_ANNOVAR} ./tempfiles/${i}_F3.recode.vcf ${humandb} -buildver hg38 -out ./vcf_anno/${i}_SNVs.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
done

echo "Analysis Finished"
