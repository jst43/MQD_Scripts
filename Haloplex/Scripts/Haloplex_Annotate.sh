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

if [ ! -d ${filepath}tempfiles ]; then
        echo "Directory tempfiles not found"
        exit 1
fi

if [ ! -d ${filepath}pindeltemp ]; then
        echo "Directory pindeltemp not found"
        exit 1
fi

if [ ! -d ${filepath}hotspot_vcf ]; then
        echo "Directory hotspot_vcf not found"
        exit 1
fi

#SCRIPT

mkdir ${filepath}snv_anno
mkdir ${filepath}pindel_anno
mkdir ${filepath}hotspot_anno

while read pfx; do
	for directory in tempfiles pindeltemp hotspot_vcf; do
		echo $pfx
		if [ $directory == "tempfiles" ]; then
			suffix=".vcf"
			outdir="snv_anno"
			new_suffix="_snv.tsv"
		elif [ $directory == "pindeltemp" ]; then
			suffix="_DSI.bedfiltered.sorted.vcf"
			outdir="pindel_anno"
			new_suffix="_DSI.tsv"
		elif [ $directory == "hotspot_vcf" ]; then
			suffix="_hotspot.vcf"
			outdir="hotspot_anno"
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
	done
done <${filepath}samples_noLane.txt
