#!/bin/bash
echo "######################################################"
echo "#                   lympHaloPlex_v0.1                #"
echo "#    Writen by Eguzkine Ochoa (eguzki8a@gmail.com)   #"
echo "#                                                    #"
echo "#                   November 3rd 2016                #"
echo "# This Bash script uses the following software under #"
echo "#    GNU Public license v2: vim, fastqc, cutadapt,   #"
echo "#      bwa, samtools, GATK, vcftools and Annovar.    #"
echo "#                                                    #"
echo "######################################################"
SureCallTRIMMER="/media/eguz/darwin/Resources/Software/SurecallTrimmer_v3.5.1.46.jar"
java8="/media/eguz/darwin/Resources/Software/jre1.8.0_112/bin/java"

##BEFORE TO START: In the folder should be fastq.gz files.
##Remove S code from all samples
ls *_L001_R1_001.fastq.gz > Scode.txt
vim -c "%s/_L001_R1_001.fastq.gz//g|wq" Scode.txt
vim -c "%s/_S\d\+//g|wq" Scode.txt
for k in `cat Scode.txt`; do
	mv ${k}_S*_L001_R1_001.fastq.gz ${k}_L001_R1_001.fastq.gz
	mv ${k}_S*_L001_R2_001.fastq.gz ${k}_L001_R2_001.fastq.gz
done
##Remove S code from Index files
ls *_L001_I1_001.fastq.gz > Scode_index.txt
vim -c "%s/_L001_I1_001.fastq.gz//g|wq" Scode_index.txt
vim -c "%s/_S\d\+//g|wq" Scode_index.txt
for k in `cat Scode_index.txt`; do
	mv ${k}_S*_L001_I1_001.fastq.gz ${k}_L001_I1_001.fastq.gz
	mv ${k}_S*_L001_I2_001.fastq.gz ${k}_L001_I2_001.fastq.gz
done

#Create list of sample names including Dup1-Dup2
ls *_L001_R1_001.fastq.gz > samples.txt
vim -c "%s/_L001_R1_001.fastq.gz//g|wq" samples.txt

#Create list of file prefixes to allow renaming of SureCallTRIMMER output files
ls *_R*.fastq.gz > pretrimNames.txt
sed -i s/.fastq.gz// pretrimNames.txt

for k in `cat samples.txt`; do
#Check sequencing quality
	fastqc ${k}_L001_R1_001.fastq.gz
	fastqc ${k}_L001_R2_001.fastq.gz
#Trimmed: Amplicon sizes should be 190-640 bp.SureCall processes the read sequences to trim low quality bases from the ends, remove adaptor sequences, and mask enzyme footprints (for HaloPlex).
	${java8} -Xmx250g -jar ${SureCallTRIMMER} -fq1 /media/eguz/darwin/Analysis/Ongoing/HALOPLEX/DATA_HiSeq/${k}_L001_R1_001.fastq.gz -fq2 /media/eguz/darwin/Analysis/Ongoing/HALOPLEX/DATA_HiSeq/${k}_L001_R2_001.fastq.gz -hs
done

#Rename SureCallTRIMMER output files
for k in `cat pretrimNames.txt`; do
	mv ${k}.*.fastq.gz ${k}.trimmed.fastq.gz;
done

