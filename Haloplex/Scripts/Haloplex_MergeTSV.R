#!/usr/bin/Rscript

load_packages <- function(){
  library(tidyverse)
  library(reshape2)
}


read_prefixes <- function(filepath){
  prefixes <- scan(paste0(filepath, 'samples_noLane.txt'),
              character())
  return(prefixes)
}


find_tsv_files <- function(prefixes, filepath){
  hotspot_files <- Sys.glob(paste0(filepath,
                                   "hotspot/annotated_vcf/",
                                   prefixes,
                                   "*.tsv"))
  snv_files <- Sys.glob(paste0(filepath,
                               "snv/annotated_vcf/",
                               prefixes,
                               "*.tsv"))
  pindel_files <- Sys.glob(paste0(filepath,
                                  "pindel/annotated_vcf/",
                                  prefixes,
                                  "*.tsv"))
  return(c(hotspot_files,
           snv_files,
           pindel_files))
}


find_vcf_files <- function(prefixes, filepath){
  hotspot_files <- Sys.glob(paste0(filepath,
                                   "hotspot/tempfiles/",
                                   prefixes,
                                   "*.vcf"))
  snv_files <- Sys.glob(paste0(filepath,
                               "snv/tempfiles/",
                               prefixes,
                               "*.vcf"))
  pindel_files <- Sys.glob(paste0(filepath,
                                  "pindel/tempfiles/",
                                  prefixes,
                                  "*sorted.vcf"))
  return(c(hotspot_files,
           snv_files,
           pindel_files))
}


create_file_frame <- function(filepath){
  prefixes <- read_prefixes(filepath)
  vcf_filenames <- find_vcf_files(prefixes, filepath)
  tsv_filenames <- find_tsv_files(prefixes, filepath)
  fileframe <- data.frame(pfx=rep(prefixes, 3),
                          tsv=tsv_filenames,
                          vcf=vcf_filenames,
                          stringsAsFactors = FALSE) %>%
    arrange(pfx)
  return(fileframe)
}


get_unique_vcf_format <- function(vcf_filename){
  unique_format <- unique(system(paste0("grep -v '#' ",
                                        vcf_filename,
                                        " | awk '{ print $9 }'"),
                                 intern=TRUE))
  return(unique_format)
}


get_all_vcf_cols <- function(vcf_filenames){
  columns <- character()
  for(fname in vcf_filenames){
    unique_format <- get_unique_vcf_format(fname)
    columns <- c(columns,
                 unlist(strsplit(unique_format, ":")))
  }
  columns <- unique(columns)
  return(columns)
}


read_annotation <- function(tsv_filename, prefix){
  tsv <- read_tsv(tsv_filename)
  tsv <- cbind(Sample=rep(prefix, nrow(tsv)),
               checkBamFile=rep('',nrow(tsv)),
               tsv,
               Hotspot=rep(0, nrow(tsv)),
               SNV=rep(0, nrow(tsv)),
               Pindel=rep(0, nrow(tsv)))
  if(grepl("hotspot", tsv_filename, fixed=T)){
    tsv$Hotspot <- 1
  }else if(grepl("snv", tsv_filename, fixed=T)){
    tsv$SNV <- 1
  }else if(grepl("pindel", tsv_filename, fixed=T)){
    tsv$Pindel <- 1
  }
  Loc <- colsplit(tsv$Location,
                  ":",
                  names=c('Chr', 'Location'))
  Pos <- colsplit(Loc$Location,
                  "-",
                  names=c('Start', 'End'))
  Pos$End[which(is.na(Pos$End))] <- Pos$Start[which(is.na(Pos$End))]
  variant_id <- paste0('chr',
                       tsv$Location,
                       '_',
                       tsv$USED_REF,
                       '>',
                       tsv$Allele)
  tsv <- cbind(tsv,
               Chr=Loc$Chr,
               Pos,
               variant_id)
  return(tsv)
}


add_vcf_cols <- function(tsv, vcf_filename, unique_format){
  for(cname in unique_format){
    tsv[[cname]] <- rep(NA, nrow(tsv))
  }
  specific_cols <- system(paste0("grep -v '#' ",
                                 vcf_filename,
                                 " | awk '{ print $9 }'"),
                          intern=TRUE)
  specific_vals <- system(paste0("grep -v '#' ",
                                 vcf_filename,
                                 " | awk '{ print $10 }'"),
                          intern=TRUE)
  for(i in 1:nrow(tsv)){
    column_data <- unlist(strsplit(specific_cols[i], ':'))
    values_data <- unlist(strsplit(specific_vals[i], ':'))
    for(c in 1:length(column_data)){
      tsv[i, column_data[c]] <- values_data[c]
    }
  }
  return(tsv)
}


cat_annotations <- function(fileframe, unique_format){
  print(paste("Reading in",
              fileframe$tsv[1]))
  tsv <- read_annotation(fileframe$tsv[1], fileframe$pfx[1])
  tsv <- add_vcf_cols(tsv, fileframe$vcf[1], unique_format)
  for(i in 2:nrow(fileframe)){
    print(paste("Reading in",
                fileframe$tsv[i]))
    temp_tsv <- read_annotation(fileframe$tsv[i], fileframe$pfx[i])
    temp_tsv <- add_vcf_cols(temp_tsv, fileframe$vcf[i], unique_format)
    tsv <- rbind(tsv, temp_tsv)
  }
  return(tsv)
}


merge_tsv <- function(tsv){
  tsv <- tsv %>%
    arrange(Sample, Chr, Start)
  return(tsv)
}


calc_VAF <- function(tsv){
  AD_Info <- colsplit(tsv$AD,
                      ',',
                      names=c('AD.Ref',
                              'AD.Alt1',
                              'AD.Alt2'))
  tsv$DP <- rowSums(AD_Info,
                    na.rm=TRUE)
  VAF <- AD_Info / tsv$DP
  colnames(VAF) <- c('VAF.Ref',
                     'VAF.Alt1',
                     'VAF.Alt2')
  tsv <- cbind(tsv, AD_Info, VAF)
  return(tsv)
}


remove_unnecessary_cols <- function(tsv, unique_format){
  trimmed_vcf_cols <- unique_format[(unique_format != "AD") &
                                      (unique_format != "AF") &
                                      (unique_format != "DP")]
  tsv <- tsv %>%
    select(Sample,
           Location,
           Uploaded_variation,
           Existing_variation,
           USED_REF,
           Allele,
           checkBamFile,
           SYMBOL,
           AD.Ref,
           AD.Alt1,
           AD.Alt2,
           DP,
           VAF.Ref,
           VAF.Alt1,
           VAF.Alt2,
           VARIANT_CLASS,
           Consequence,
           HGVSc,
           HGVSp,
           AFR_AF,
           AMR_AF,
           EAS_AF,
           EUR_AF,
           SAS_AF,
           gnomAD_AF,
           gnomAD_AFR_AF,
           gnomAD_AMR_AF,
           gnomAD_ASJ_AF,
           gnomAD_EAS_AF,
           gnomAD_EAS_AF,
           gnomAD_FIN_AF,
           gnomAD_NFE_AF,
           gnomAD_OTH_AF,
           gnomAD_SAS_AF,
           CLIN_SIG,
           clinvar_clnsig,
           CADD_phred,
           CADD_raw,
           FATHMM_pred,
           FATHMM_score,
           `GERP++_RS`,
           LRT_pred,
           LRT_score,
           MetaLR_pred,
           MetaLR_score,
           MetaSVM_pred,
           MetaSVM_score,
           MutationAssessor_pred,
           MutationAssessor_score,
           MutationTaster_pred,
           MutationTaster_score,
           Polyphen2_HDIV_pred,
           Polyphen2_HDIV_score,
           Polyphen2_HVAR_pred,
           Polyphen2_HVAR_score,
           SIFT_pred,
           SIFT_score,
           SiPhy_29way_logOdds,
           VEST3_score,
           phyloP100way_vertebrate,
           SNV,
           Pindel,
           Hotspot,
           Chr,
           Start,
           End,
           variant_id,
           trimmed_vcf_cols)
  return(tsv)
}


get_output_name <- function(filepath){
  directories <- strsplit(filepath, '/')[[1]]
  name <- directories[length(directories)]
  return(name)
}


main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)!=1){
    stop("Haloplex_MergeTSV.R only accepts one input")
  }
  if(args[1]=="h" | args[1]=="help"){
    stop("Usage: ./Haloplex_MergeTSV.R /path/to/working/directory/")
  }
  filepath <- args[1]
  if(substr(filepath, start=nchar(filepath), stop=nchar(filepath)) != '/'){
    filepath <- paste0(filepath, '/')
  }
  print(paste("filepath is", filepath))
  load_packages()
  fileframe <- create_file_frame(filepath)
  unique_format <- get_all_vcf_cols(fileframe$vcf)
  tsv <- cat_annotations(fileframe, unique_format)
  tsv <- merge_tsv(tsv)
  tsv <- calc_VAF(tsv)
  tsv <- remove_unnecessary_cols(tsv, unique_format)
  name <- get_output_name(filepath)
  write_tsv(tsv, paste0(filepath, 'Output/', name, '.tsv'))
}

main()
