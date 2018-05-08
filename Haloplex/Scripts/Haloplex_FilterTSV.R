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
  synonymous_variants <- grepl('synonymous', tsv$Consequence)
  snvs <- tsv$VARIANT_CLASS == 'SNV'
  keep_rows <- (!synonymous_variants) | (!snvs)
  tsv <- tsv[keep_rows,]
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print('Filtering 1000G Results')
  row1 <- row2
  cosmic <- grepl('COSM', tsv$Existing_variation)
  tsv <- subset(tsv, EUR_AF=="-" | as.numeric(EUR_AF)<=0.01 | (as.numeric(EUR_AF)>0.01 & cosmic))
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print("Removing intronic variants")
  row1 <- row2
  tsv <- tsv %>% filter(HGVSc != "-")
  introns <- !grepl('intron', tsv$Consequence)
  tsv <- tsv[introns,]
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  print("Filtering Allele Depth")
  row1 <- row2
  tsv <- tsv %>% filter((Pindel != 1) | (Pindel == 1 & (AD.Alt1>=5 | AD.Alt2>=5)))
  row2 <- nrow(tsv)
  dropped_variants(row1, row2)
  return(tsv)
}


add_freq_reprod <- function(tsv){
  tsv <- cbind(tsv, Frequency=rep(0,nrow(tsv)), Reproducible=rep(0,nrow(tsv)))
  tsv$Frequency <- unlist(lapply(tsv$HGVSc, function(x) sum(tsv$HGVSc %in% x)))
  tsvDupVec <- paste0(gsub('.{1}$', '', tsv$Sample), tsv$HGVSc)
  tsv$Reproducible <- unlist(lapply(tsvDupVec, function(x) sum(tsvDupVec %in% x)>1))
  tsv$Reproducible[!grepl("dup", tolower(tsv$Sample))] <- NA
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
    stop("Usage: ./Haloplex_FilterTSV.R path/to/input.tsv")
  }
  tsv_path <- args[1]
  print(paste("TSV Path is", tsv_path))
  load_packages()
  tsv <- read_in_tsv(tsv_path)
  tsv <- filtering(tsv)
  tsv <- add_freq_reprod(tsv)
  write_to_tsv(tsv, tsv_path)
}


main()