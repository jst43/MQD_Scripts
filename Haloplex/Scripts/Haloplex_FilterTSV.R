#!/usr/bin/Rscript


load_packages <- function(){
  library(readr)
  library(reshape2)
  library(tidyverse)
}


read_in_tsv <- function(tsv_path){
  tsv <- read_delim(tsv_path,
                    "\t",
                    escape_double = FALSE,
                    trim_ws = TRUE,
                    col_types = cols(Chr=col_character(),
                                     DP=col_double()))
  return(tsv)
}


dropped_variants <- function(nrow1, nrow2){
  n_var <- nrow1 - nrow2
  perc_var <- (n_var / nrow1) * 100
  print(paste0(n_var, " variants removed (", perc_var, "%)"))
}


filtering <- function(tsv){
  print("Removing synonymous SNVs")
  row1 <- nrow(tsv)
  tsv <- tsv %>% filter(!grepl('synonymous', Consequence) | VARIANT_CLASS!='SNV')
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print('Filtering 1000G Results')
  row1 <- row2
  tsv <- tsv %>%
    filter(EUR_AF=='-' | as.numeric(EUR_AF)<=0.01 | as.numeric(EUR_AF)>=0.1 | grepl('COSM', Existing_variation))
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print('Filtering ExAC Results')
  row1 <- row2
  tsv <- tsv %>% filter(gnomAD_AF<0.0001 | grepl('COSM', Existing_variation))
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print("Removing intronic variants")
  row1 <- row2
  tsv <- tsv %>% filter(!grepl('intron', Consequence),
                        !grepl('intergenic', Consequence),
                        !grepl('downstream', Consequence))
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print("Filtering Allele Depth")
  row1 <- row2
  tsv <- tsv %>% filter((Pindel != 1) | (Pindel == 1 & (AD.Alt1>=5 | AD.Alt2>=5)))
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print('Removing non-heterogenous SNPs')
  row1 <- row2
  tsv <- tsv %>% filter(!grepl('rs', Existing_variation) | GT=='0/1' | EUR_AF<0.1 | Pindel==1 | VARIANT_CLASS!='SNV')
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  return(tsv)
}


check_for_assoc <- function(tsv, row_index){
  d <- dim(tsv %>% filter(Sample==tsv$Sample[row_index] &
                            Chr==tsv$Chr[row_index] &
                            SNP!=1 &
                            GT=='0/1' &
                            Start>=(tsv$Start[row_index] - 2000) &
                            End<=(tsv$End[row_index] + 2000)))
  if(d[1]==0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}


find_SNP_assoc <- function(tsv, row_index){
  associated_variants <- which(tsv$Sample==tsv$Sample[row_index] &
                                          tsv$Chr==tsv$Chr[row_index] &
                                          tsv$SNP!=1 &
                                          tsv$GT=='0/1' &
                                          tsv$Start>=(tsv$Start[row_index] - 2000) &
                                          tsv$End<=(tsv$End[row_index] + 2000))
  if(!is_empty(associated_variants)){
    assoc_message <- paste0('Associated with SNP ',
                            tsv$Existing_variation[row_index],
                            ' (row ',
                            row_index,
                            ')')
    tsv$Assoc.SNP[associated_variants] <- assoc_message
  }
  return(tsv)
}


add_SNP_annotation <- function(tsv){
  tsv <- data.frame(tsv,
                    SNP=rep(0, nrow(tsv)),
                    Assoc.SNP=rep('-', nrow(tsv)),
                    stringsAsFactors=FALSE)
  snp_indexes <- which(grepl('rs', tsv$Existing_variation) &
                         (tsv$GT=='0/1') &
                         (tsv$EUR_AF>=0.1) &
                         (tsv$VARIANT_CLASS=='SNV'))
  tsv$SNP[snp_indexes] <- 1
  assoc_bool <- rep(TRUE, nrow(tsv))
  for(index in snp_indexes){
    assoc_bool[index] <- check_for_assoc(tsv, index)
  }
  print('Removing SNPs unassociated with variants')
  row1 <- nrow(tsv)
  tsv <- tsv[assoc_bool,]
  dropped_variants(row1, nrow(tsv))
  snp_indexes <- which(tsv$SNP==1)
  for(index in snp_indexes){
    tsv <- find_SNP_assoc(tsv, index)
  }
  return(tsv)
}


get_CNV_filepaths <- function(filepath){
  cnv_paths <- Sys.glob(paste0(filepath,
                               'CNV/*.csv'))
  return(cnv_paths)
}


get_CNV_info <- function(cnv_path, tsv){
  cnv <- read_csv(cnv_path)
  cnv_name <- gsub(paste0(filepath, 'CNV/'),
                   '',
                   cnv_path)
  cnv_name <- gsub('_CNVs.csv',
                   '',
                   cnv_name)
  for(i in 1:nrow(cnv)){
    associated_indexes <- which(grepl(cnv_name, tsv$Sample) &
                                  tsv$Chr==cnv$chromosome[i] &
                                  tsv$Start>=cnv$start[i] &
                                  tsv$End<=cnv$end[i])
    tsv$cnv_id[associated_indexes] <- cnv$id[i]
    tsv$cnv_type[associated_indexes] <- cnv$type[i]
    tsv$cnv_ratio[associated_indexes] <- cnv$reads.ratio[i]
    tsv$cnv_BF[associated_indexes] <- cnv$BF[i]
  }
  return(tsv)
}


add_CNV_annotation <- function(tsv, filepath){
  tsv <- data.frame(tsv,
                    cnv_id=rep('', nrow(tsv)),
                    cnv_type=rep('', nrow(tsv)),
                    cnv_ratio=rep('', nrow(tsv)),
                    cnv_BF=rep(0, nrow(tsv)),
                    stringsAsFactors=FALSE)
  cnv_paths <- get_CNV_filepaths(filepath)
  for(cnv_path in cnv_paths){
    tsv <- get_CNV_info(cnv_path, tsv)
  }
  return(tsv)
}


add_freq_reprod <- function(tsv){
  tsv <- cbind(tsv, Frequency=rep(0,nrow(tsv)), Reproducible=rep(0,nrow(tsv)))
  tsv$Frequency <- unlist(lapply(tsv$variant_id, function(x) sum(tsv$variant_id %in% x)))
  tsvDupVec <- paste0(gsub('.{1}$', '', tsv$Sample), tsv$variant_id)
  tsv$Reproducible <- unlist(lapply(tsvDupVec, function(x) sum(tsvDupVec %in% x)>1))
  tsv$Reproducible[!grepl("dup", tolower(tsv$Sample))] <- NA
  return(tsv)
}


rearrange_tsv <- function(tsv){
  col_loc1 <- which(colnames(tsv)=='VAF.Alt2')
  col_loc2 <- which(colnames(tsv)=='VARIANT_CLASS')
  col_loc3 <- which(colnames(tsv)=='Hotspot')
  col_loc4 <- which(colnames(tsv)=='Chr')
  col_loc5 <- length(colnames(tsv)) - 8
  tsv <- data.frame(tsv[,1:col_loc1],
                    Frequency=tsv$Frequency,
                    Reproducible=tsv$Reproducible,
                    Assoc.SNP=tsv$Assoc.SNP,
                    cnv_type=tsv$cnv_type,
                    cnv_ratio=tsv$cnv_ratio,
                    tsv[,col_loc2:col_loc3],
                    SNP=tsv$SNP,
                    tsv[,col_loc4:col_loc5],
                    cnv_id=tsv$cnv_id,
                    cnv_BF=tsv$cnv_BF,
                    stringsAsFactors=FALSE)
  return(tsv)
}


write_to_tsv <- function(tsv, tsv_path){
  print("Writing tsv")
  write.table(tsv,
              file=tsv_path,
              quote=FALSE, sep="\t", row.names=FALSE)
}


main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)!=1){
    stop("Haloplex_FilterTSV.R only accepts one input")
  }
  if(args[1]=="h" | args[1]=="help"){
    stop("Usage: ./Haloplex_FilterTSV.R path/to/working/directory")
  }
  filepath <- args[1]
  if(substr(filepath, start=nchar(filepath), stop=nchar(filepath)) != '/'){
    filepath <- paste0(filepath, '/')
  }
  print(paste("filepath is", filepath))
  load_packages()
  tsv_path <- Sys.glob(paste0(filepath, 'Output/*tsv'))
  tsv <- read_in_tsv(tsv_path)
  tsv <- filtering(tsv)
  tsv <- add_SNP_annotation(tsv)
  tsv <- add_CNV_annotation(tsv, filepath)
  tsv <- add_freq_reprod(tsv)
  tsv <- rearrange_tsv(tsv)
  write_to_tsv(tsv, tsv_path)
}


main()
