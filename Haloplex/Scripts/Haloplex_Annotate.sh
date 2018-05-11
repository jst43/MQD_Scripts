#!/bin/bash

#CONSTANTS
VEP="../../../Software/ensembl-vep/vep"
VEP_dir="../../../Software/ensembl-vep/VEP_cache"

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
	filepath=`pwd`"/"
fi

if [ ! -d $filepath ]; then
	echo "This script requires a valid filepath to the fastq file directory"
	exit 1
fi

if [ ! -d ${filepath}snv/tempfiles ]; then
        echo "Directory tempfiles not found"
        exit 1
fi

if [ ! -d ${filepath}pindel/tempfiles ]; then
        echo "Directory pindeltemp not found"
        exit 1
fi

if [ ! -d ${filepath}hotspot/tempfiles ]; then
        echo "Directory hotspot_vcf not found"
        exit 1
fi

#SCRIPT

mkdir ${filepath}snv/annotated_vcf
mkdir ${filepath}pindel/annotated_vcf
mkdir ${filepath}hotspot/annotated_vcf
mkdir ${filepath}Output

while read pfx; do
	for directory in snv/tempfiles pindel/tempfiles hotspot/tempfiles; do
		echo $pfx
		if [ $directory == "snv/tempfiles" ]; then
			suffix=".vcf"
			outdir="snv/annotated_vcf"
			new_suffix="_snv.tsv"
		elif [ $directory == "pindel/tempfiles" ]; then
			suffix="_DSI.bedfiltered.sorted.vcf"
			outdir="pindel/annotated_vcf"
			new_suffix="_DSI.tsv"
		elif [ $directory == "hotspot/tempfiles" ]; then
			suffix="_hotspot.vcf"
			outdir="hotspot/annotated_vcf"
			new_suffix="_hotspot.tsv"
		fi
		$VEP -i ${filepath}${directory}/${pfx}${suffix} \
			--dir $VEP_dir \
			--cache \
			--force_overwrite \
			--tab \
			--merged \
			--variant_class \
			--hgvs \
			--symbol \
			--check_existing \
			--af_1kg \
			--af_gnomad \
			--pick \
			--plugin dbNSFP,${VEP_dir}/dbNSFP.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,VEST3_score,CADD_raw,CADD_phred,GERP++_RS,phyloP100way_vertebrate,SiPhy_29way_logOdds,clinvar_clnsig,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred \
			-o ${filepath}${outdir}/${pfx}${new_suffix}
		sed -i 's|#Uploaded_variation|Uploaded_variation|' ${filepath}${outdir}/${pfx}${new_suffix}
		grep -v "#" ${filepath}${outdir}/${pfx}${new_suffix} > ${filepath}${outdir}/temp.tsv
		mv ${filepath}${outdir}/temp.tsv ${filepath}${outdir}/${pfx}${new_suffix}
	done
done <${filepath}samples_noLane.txt

./Haloplex_MergeTSV.R ${filepath}
