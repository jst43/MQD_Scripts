#!/bin/bash
hg38="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.fa"
GATK="/media/eguz/darwin/Resources/Software/GenomeAnalysisTK.jar"
All="/media/eguz/darwin/Resources/hg38.p6/All.vcf"
TABLE_ANNOVAR="/media/eguz/darwin/Resources/Software/annovar/table_annovar.pl"
humandb="/media/eguz/darwin/Resources/Software/annovar/humandb"
PICARD="/media/eguz/darwin/Resources/Software/picard-tools-1.141/picard.jar"

#Change names outputs trimming software!!!!
for k in `cat samples.txt`; do
	#Decompress
	gunzip ${k}_L001_R1_001.trimmed.fastq.gz
	gunzip ${k}_L001_R2_001.trimmed.fastq.gz
	#Alignment 1000Genomes(Hg38)
	bwa_0.7.12 mem -R "@RG\tID:<${k}>\tLB:LIBRARY_NAME\tSM:<${k}>\tPL:ILLUMINA" ${hg38} ${k}_L001_R1_001.trimmed.fastq ${k}_L001_R2_001.trimmed.fastq > ${k}.sam
	#Remove Duplicates with LocatIt
	/media/eguz/darwin/Resources/Software/jre1.8.0_112/bin/java -Xmx250g -jar /media/eguz/darwin/Resources/Software/LocatIt_v3.5.1.46.jar -X /media/eguz/darwin/Analysis/Ongoing/HALOPLEX/DATA -U -IS -OB -b amplicons38.bed -o ${k}_RMD ${k}.sam ${k}_L001_I2_001.fastq.gz
	#Convert bam without duplicates in fastq file
	java -Xmx250g -jar ${PICARD} SamToFastq I=${k}_RMD.bam F=${k}_R1.fastq F2=${k}_R2.fastq
	#Create bam file, sort + index
	samtools_1.2 sort ${k}_RMD.bam ${k}.sorted
	samtools_1.2 index ${k}.sorted.bam
	#Realigned and Indels
	java -Xmx250g -jar ${GATK} -nt 20 -T RealignerTargetCreator -R ${hg38} -I ${k}.sorted.rg.bam -o ${k}.bam.list
	java -Xmx250g -jar ${GATK} -T IndelRealigner -R ${hg38} -I ${k}.sorted.rg.bam -targetIntervals ${k}.bam.list -o ${k}.sorted.realigned.bam
	#Recalibrator and quality control
	java -Xmx250g -jar ${GATK} -nct 20 -T BaseRecalibrator -R ${hg38} -I ${k}.sorted.realigned.bam -l info -knownSites ${All} -o ${k}.sorted.realigned.table
	java -Xmx250g -jar ${GATK} -nct 20 -T PrintReads -R ${hg38} -I ${k}.sorted.realigned.bam -l INFO -BQSR ${k}.sorted.realigned.table -o ${k}.sorted.realigned.recal.bam
	java -Xmx250g -jar ${GATK} -T DepthOfCoverage -R ${hg38} -o ${k}.coverage -I ${k}.sorted.realigned.recal.bam -L amplicons38_2.bed
	#Calling variants
	java -Xmx250g -jar ${GATK} -T UnifiedGenotyper -R ${hg38} -I ${k}.sorted.realigned.recal.bam -glm BOTH --dbsnp ${All} -stand_call_conf 30.0 -stand_emit_conf 10.0 -A Coverage -dcov 10000 -A AlleleBalance --max_alternate_alleles 40 -o ${k}.vcf
	#Filtering	
	vcftools_0.1.13 --vcf ${k}.vcf --minQ 30 --recode --out ${k}_F1
	vcftools_0.1.13 --vcf ${k}_F1.recode.vcf --min-meanDP 50 --recode --out ${k}_F2
	#Annotation
	perl ${TABLE_ANNOVAR} ${k}.vcf ${humandb} -buildver hg38 -out ${k}_SNVs.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
done

Building result files
ls *_SNVs.myanno.hg38_multianno.txt > libraryHALOnames_SNVs.txt
for i in `cat libraryHALOnames_SNVs.txt`; do
	sed -i 's|^|${i}\t|' ${i}
done

cat *_SNVs.myanno.hg38_multianno.txt > LibraryHALO_SNVs_hg38_multianno.txt
sed -i 's|_SNVs.myanno.hg38_multianno.txt||' LibraryHALO_SNVs_hg38_multianno.txt
