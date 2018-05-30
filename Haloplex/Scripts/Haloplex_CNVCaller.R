#!/usr/bin/Rscript

load_packages <- function(){
  library(readr)
  library(methods)
  library(tidyverse)
  library(devtools)
  load_all('/mnt/raid/Resources/Software/ExomeDepth/')
  options(warn=-1)
}


load_amplicons <- function(){
  amplicons <- read_delim("/mnt/raid/Resources/MQD_Scripts/Haloplex/Dependent_Files/coverage_amplicons38.bed",
                        "\t",
                        escape_double = FALSE,
                        col_names = FALSE,
                        col_types = cols(X5 = col_skip()),
                        trim_ws = TRUE)
  colnames(amplicons) <- c("chromosome", "start", "end", "name", "strand")
  ampliconID <- paste(amplicons$chromosome,amplicons$start,amplicons$end, sep="-")
  amplicons <- cbind(amplicons, ID=ampliconID)
  amplicons$ID <- as.character(amplicons$ID)
  nUnique <- length(unique(amplicons$ID))
  keepRows <- numeric(length=nUnique)
  for(i in 1:nUnique){
    dupRows <- which(amplicons$ID %in% unique(amplicons$ID)[i])
    keepRows[i] <- dupRows[1]
    }
  amplicons <- amplicons[keepRows,]
  rownames(amplicons) <- seq(1,nUnique)
  amplicons <- amplicons[,1:4]
  #Remove Amplicons in ChrX as sex of patients unknown
  amplicons <- amplicons[which(amplicons$chromosome!="X"),]
  return(amplicons)
}


get_bam_paths <- function(filepath){
  bampaths <- Sys.glob(paste0(filepath, "realigned_recal_bam/*bam"))
  baipaths <- Sys.glob(paste0(filepath, "realigned_recal_bam/*bai"))
  bamnames <- gsub(paste0(filepath, "realigned_recal_bam/"), "", bampaths, fixed=TRUE)
  bainames <- gsub(paste0(filepath, "realigned_recal_bam/"), "", baipaths, fixed=TRUE)
  fileframe <- data.frame(bampath=bampaths,
                          baipath=baipaths,
                          bam=bamnames,
                          bai=bainames,
                          stringsAsFactors=FALSE)
  return(fileframe)
}


get_bam_count_frame <- function(amplicons, fileframe){
  bamcount <- getBamCounts(bed.frame=amplicons,
                           bam.files=fileframe$bampath,
                           index.files=fileframe$baipath,
                           include.chr=FALSE)
  bamcount.dafr <- as.data.frame(bamcount[,colnames(bamcount)])
  return(bamcount.dafr)
}


get_sample_names <- function(filepath){
  sample_names <- scan(paste0(filepath, 'samples_noLane.txt'),
                       character())
  sample_names <- gsub('_Dup[0-9]', '', sample_names)
  sample_names <- unique(sample_names)
  return(sample_names)
}


runTest <- function(bamsample, fileframe, bamcount, filepath){
  #Find Columns in BamCount.dafr that have belong to BamFile sample
  matchbams <- grepl(bamsample, fileframe$bam)
  mismatchbams <- !matchbams
  if(any(mismatchbams==TRUE)){
  	#Generate test coverage from average of MatchBams
  	my_test <- ceiling(rowMeans(as.data.frame(bamcount[, c(rep(FALSE, 5), matchbams)])))
  	#Generate list of reference samples
  	my_ref_samples <- fileframe$bam[mismatchbams]
  	my_reference_set <- as.matrix(bamcount[, my_ref_samples])
  	#Find appropriate reference bams
  	my_choice <- select.reference.set(test.counts = my_test,
  	                                  reference.counts = my_reference_set,
        	                            bin.length = (bamcount$end - bamcount$start)/ 1000,
                	                    n.bins.reduced = 10000)
  	#Create reference set to compare to test set
  	my_matrix <- as.matrix(bamcount[, my_choice$reference.choice, drop=FALSE])
  	my_reference_selected <- apply(my_matrix, MARGIN=1, FUN=sum)
  	#Fit the beta-binomial model to the data
  	all_regions <- new('ExomeDepth',
        	             test=my_test,
                	     reference=my_reference_selected,
        	             formula='cbind(test,reference) ~ 1')
  	#Call CNVs
  	all_exons <- CallCNVs(x=all_regions,
        	                transition.probability = 10^-4,
                	        chromosome=bamcount$space,
                	        start=bamcount$start,
                	        end=bamcount$end,
        	                name=bamcount$names)
  	if(all_exons[[1]] >= 0.97){
  	  all_exons <- all_exons[[2]]
  	}else{
  	  print('Correlation is below threshold')
  	  return()
  	}
  	my_cnvs <- all_exons@CNV.calls[order(all_exons@CNV.calls$BF, decreasing = TRUE),]
  	my_cnvs <- my_cnvs %>% filter(BF >= 55)
  	if(dim(my_cnvs)[1]!=0){
  	  outputFile <- paste0(filepath,
  	                       "CNV/",
  	                       bamsample,
  	                       "_CNVs.csv")
  	  write.csv(my_cnvs,
  	            outputFile,
  	            quote=FALSE,
  	            row.names=FALSE)
  	}
   }
}


main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)!=1){
    stop("Haloplex_CNVCaller.R only accepts one input")
  }
  if(args[1]=="h" | args[1]=="help"){
    stop("Usage: ./Haloplex_CNVCaller.R /path/to/working/directory/")
  }
  filepath <- args[1]
  if(substr(filepath, start=nchar(filepath), stop=nchar(filepath)) != '/'){
    filepath <- paste0(filepath, '/')
  }
  print(paste("filepath is", filepath))
  load_packages()
  amplicons <- load_amplicons()
  fileframe <- get_bam_paths(filepath)
  sample_names <- get_sample_names(filepath)
  bamcount <- get_bam_count_frame(amplicons, fileframe)
  dir.create(paste0(filepath,'CNV'),
             showWarnings = FALSE)
  for(bam in sample_names){
    print(bam)
    runTest(bam, fileframe, bamcount, filepath)
  }
}


main()
