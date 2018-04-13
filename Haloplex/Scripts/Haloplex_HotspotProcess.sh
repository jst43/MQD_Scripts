#!/bin/bash

GATK_v3_8="/mnt/raid/Resources/Software/GenomeAnalysisTK.jar"
GATK_v3_5="/home/joe/GenomeAnalysisTK.jar"
hg38="/mnt/raid/Resources/hg38.p6/hg38_2MergeAll.fa"
All="/mnt/raid/Resources/hg38.p6/All.vcf"
dbsnp="/mnt/raid/Resources/hg38.p6/All.vcf"
TABLE_ANNOVAR="/mnt/raid/Resources/Software/annovar/table_annovar.pl"
humandb="/mnt/raid/Resources/Software/annovar/humandb/"

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

mkdir MuTectBam_In
mkdir MuTectBam_Out
mkdir MuTectVCF

while read pfx; do
	java -Xmx40g -jar ${GATK_v3_5} -nct 20 -T BaseRecalibrator -R $hg38 -I tempfiles/${pfx}.sorted.bam -l info -knownSites $All -o tempfiles/${pfx}.sorted.table
	java -Xmx40g -jar ${GATK_v3_5} -nct 20 -T PrintReads -R $hg38 -I tempfiles/${pfx}.sorted.bam -l INFO -BQSR tempfiles/${pfx}.sorted.table -o MuTectBam_In/${pfx}.sorted.recal.bam
	java -jar ${GATK_v3_8} -T MuTect2 -R $hg38 -I:tumor MuTectBam_In/${pfx}.sorted.recal.bam --dbsnp $dbsnp --maxReadsInRegionPerSample 10000 -A Coverage -A AlleleBalance --max_alternate_alleles 40 -bamout MuTectBam_Out/${pfx}.out.bam -o MuTectVCF/${pfx}_hotspot.vcf -L hotspotsMuTect2.list
	perl $TABLE_ANNOVAR MuTectVCF/${pfx}_hotspot.vcf $humandb -buildver hg38 -out MuTectVCF/${pfx}_hotspot.myanno -remove -protocol refGene,cytoBand,genomicSuperDups,1000g2015aug_eur,avsnp144,cosmic70,clinvar_20160302,ljb26_all -operation g,r,r,f,f,f,f,f -nastring . -vcfinput
done <samples_noLane.txt
