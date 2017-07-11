#!/bin/bash
echo "######################################################"
echo "#                lympSeq_Indels_v1.0                 #"
echo "#    Writen by Eguzkine Ochoa (eguzki8a@gmail.com)   #"
echo "#                                                    #"
echo "#                  March 2nd 2016                    #"
echo "# This Bash script uses the following software under #"
echo "#         GNU Public license v2: vim, Pindel,        #"
echo "#            vcftools, Picard and Annovar.           #"
echo "#                                                    #"
echo "######################################################"
hg38="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.fa"
hg38Dict="/media/eguz/darwin/Resources/hg38.p6/hg38_2MergeAll.dict"
TABLE_ANNOVAR="/media/eguz/darwin/Resources/Software/annovar/table_annovar.pl"
humandb="/media/eguz/darwin/Resources/Software/annovar/humandb"
PICARD="/media/eguz/darwin/Resources/Software/picard-tools-1.141/picard.jar"

#BEFORE TO START: In the folder should be included sorted.bam files and header.txt
#Create list of name samples
ls *.sorted.bam > samplesPindel.txt
vim -c "%s/.sorted.bam//g|wq" samplesPindel.txt
#Create files: pindelinput.txt
for k in `cat samplesPindel.txt`; do
	ls ${k}.sorted.bam > ${k}_pindelinput.txt
	vim -c "%s/\(\S\+\)\.sorted.bam/\1\.sorted.bam\t350\t\1/e|wq" ${k}_pindelinput.txt
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
#        vcftools_0.1.13 --vcf ${i}_SI.vcf --max-missing 1 --recode --out ${i}_SI_F1
#Merging Deletion and Insertion
	vcf-concat ${i}_D.vcf.gz ${i}_SI.vcf.gz > ${i}_D+SI.vcf
	java -jar ${PICARD} SortVcf I=${i}_D+SI.vcf O=${i}_D+SI.sorted.vcf SD=${hg38Dict}
	less ${i}_D+SI.sorted.vcf | grep -v "#" > ${i}_headless.vcf
        cat header.txt ${i}_headless.vcf > ${i}_NewHead.vcf
#Annotation hg38
	perl ${TABLE_ANNOVAR} ${i}_NewHead.vcf ${humandb} -buildver hg38 -out ${i}_D+SI.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
done
echo "lympSeq_Indels analysis was completed."

#Building result files
ls *_D+SI.myanno.hg38_multianno.txt > librarynames_D+SI.txt
for i in `cat librarynames_D+SI.txt`; do
	sed -i "s/^/${i}\t/" ${i}
done
cat *_D+SI.myanno.hg38_multianno.txt > Library_D+SI_hg38_multianno.txt
vim -c "%s/_D+SI.myanno.hg38_multianno.txt//g|wq" Library_D+SI_hg38_multianno.txt