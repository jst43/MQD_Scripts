#!/usr/bin/Rscript


read_coverage <- function(){
  print("Reading in Coverage File")
  coverage <- read.table("../Dependent_Files/Coverage_IntronsExons.csv",
                         header=TRUE,
                         quote="\"",
                         sep=",",
                         stringsAsFactors=FALSE)
  return(coverage)
}


get_coverage_files <- function(filepath){
  filenames <- Sys.glob(paste0(filepath,
                               "coverage/*coverage"))
  return(filenames)
}


make_prefixes <- function(filenames, filepath){
  prefixes <- gsub('.coverage',
                 '',
                 filenames,
                 fixed=TRUE)
  prefixes <- gsub(paste0(filepath, 'coverage/'),
                 '',
                 prefixes)
  return(prefixes)
}


import_coverage_data <- function(coverage_files, coverage, prefixes){
  print("Importing coverage data")
  for(i in 1:length(coverage_files)){
    coverage <- cbind(coverage, read.csv(file=coverage_files[i],
                                         sep="\t",
                                         stringsAsFactors=FALSE,
                                         colClasses=c("NULL", "NULL", "NULL", NA)))
    colnames(coverage) <- c(colnames(coverage)[1:(ncol(coverage)-1)], prefixes[i])
  }
  return(coverage)
}


remove_introns <- function(coverage){
  coverage <- coverage[grepl("Exon", coverage[,1]),]
  rownames(coverage) <- NULL
  return(coverage)
}


depth_at_locus <- function(coverage, filenames, filepath, prefixes){
  if(length(filenames)>1){
    coverage <- cbind(coverage[, 1:2],
                      Total_Depth=rowSums(coverage[, 3:ncol(coverage)]),
                      Ave_Depth=round(rowMeans(coverage[, 3:ncol(coverage)]), digits=2),
                      coverage[, 3:ncol(coverage)])
  } else {
    coverage <- cbind(coverage[, 1:2],
                      Total_Depth=coverage[, 3],
                      Ave_Depth=coverage[, 3],
                      coverage[, 3])
    colnames(coverage) <- c(colnames(coverage)[1:(ncol(coverage)-1)],
                            prefixes)
  }
  return(coverage)
}


n_loci_over_threshold <- function(threshold, coverage){
  sum_threshold <- numeric(length=(ncol(coverage)-3))
  for(i in 1:length(sum_threshold)){
    sum_threshold[i] <- sum(coverage[, (i+3)]>=threshold)
  }
  return(sum_threshold)
}


make_stats_matrix <- function(coverage){
  print("Generating Statistics")
  n_loci_10 <- n_loci_over_threshold(10, coverage)
  n_loci_20 <- n_loci_over_threshold(20, coverage)
  n_loci_50 <- n_loci_over_threshold(50, coverage)
  percentLoci_10 <- round((n_loci_10 / nrow(coverage)) * 100, digits=2)
  percentLoci_20 <- round((n_loci_20 / nrow(coverage)) * 100, digits=2)
  percentLoci_50 <- round((n_loci_50 / nrow(coverage)) * 100, digits=2)
  sample_ave_depth <- round(unname(colMeans(coverage[, 4:ncol(coverage)])), digits=2)
  sample_max_depth <- numeric(length=(ncol(coverage) - 3))
  for(i in 1:length(sample_max_depth)){
    sample_max_depth[i] <- max(coverage[, (i+3)])
  }
  stats_matrix <- rbind(n_loci_10,
                       percentLoci_10,
                       n_loci_20,
                       percentLoci_20,
                       n_loci_50,
                       percentLoci_50,
                       sample_ave_depth,
                       sample_max_depth)
  stats_matrix <- cbind(c("# > 10",
                         "% > 10",
                         "# > 20",
                         "% > 20",
                         "# > 50",
                         "% > 50",
                         "Ave. Depth",
                         "Max Depth"),
                       stats_matrix)
  return(stats_matrix)
}


get_samples_below_20 <- function(coverage){
  samples_below_20 <- numeric(length=nrow(coverage))
  for(i in 1:length(samples_below_20)){
    samples_below_20[i] <- sum(coverage[i, 5:ncol(coverage)]<20)
  }
  coverage <- cbind(coverage, samples_below_20)
  return(coverage)
}


get_dup_prefixes <- function(prefixes){
  prefixes <- prefixes[grepl("dup", tolower(prefixes))]
  if(length(prefixes)>1){
    prefixes <- unique(unlist(strsplit(prefixes, "_"))[seq(1, (2*length(prefixes)), 2)])
  }
  return(prefixes)
}


generate_suffix <- function(prefix, coverage){
  indices <- grep(prefix, colnames(coverage))
  full_names <- colnames(coverage)[indices]
  suffixes <- unlist(strsplit(full_names, paste0(prefix,"_")))
  suffixes <- suffixes[seq(2, length(suffixes), 2)]
  return(suffixes)
}


create_sample_fail_vec <- function(prefix, coverage){
  suffixes <- generate_suffix(prefix, coverage)
  indices <- grep(prefix, colnames(coverage))
  temp_frame <- coverage[, c(1, indices)]
  fail_vec <- character(length=nrow(coverage))
  for(i in 1:nrow(temp_frame)){
    binary_vec <- temp_frame[i, 2:ncol(temp_frame)]>=20
    if(sum(binary_vec)==length(suffixes)){
      fail_vec[i] <- "Covered"
    }
    else if(sum(binary_vec)==0){
      fail_vec[i] <- "All Failed"
    }
    else{
      fail_vec[i] <- paste0(paste(suffixes[!binary_vec], collapse=" and "), " Failed")
    }
  }
  return(fail_vec)
}


add_failure_at_locus <- function(coverage, prefixes){
  prefixes <- get_dup_prefixes(prefixes)
  print("Adding failure at Locus data")
  for(prefix in prefixes){
    print(prefix)
    fail_vec <- create_sample_fail_vec(prefix, coverage)
    coverage_names <- c(colnames(coverage), prefix)
    coverage <- cbind(coverage, fail_vec)
    colnames(coverage) <- coverage_names
  }
  return(coverage)
}


write_coverage <- function(coverage, stats_matrix, filepath){
  print("Writing to file")
  write.csv(coverage,
            file=paste0(filepath,"coverage/Coverage.csv"),
            row.names=FALSE,
            quote=FALSE)
  for(i in 1:nrow(stats_matrix)){
    line <- paste0(",,",
                   paste(stats_matrix[i,],
                         collapse=","))
    write(line,
          file=paste0(filepath,"coverage/Coverage.csv"),
          append=TRUE)
  }
  print("Done")
}


main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)!=1){
    stop("Haloplex_coverageSummary.R only accepts one input")
  }
  if(args[1]=="h" | args[1]=="help"){
    stop("Usage: ./Haloplex_coverageSummary.R /path/to/working/directory/")
  }
  filepath <- args[1]
  if(substr(filepath, start=nchar(filepath), stop=nchar(filepath)) != '/'){
    filepath <- paste0(filepath, '/')
  }
  print(paste("filepath is", filepath))
  coverage <- read_coverage()
  coverage_files <- get_coverage_files(filepath)
  prefixes <- make_prefixes(coverage_files, filepath)
  coverage <- import_coverage_data(coverage_files, coverage, prefixes)
  coverage <- remove_introns(coverage)
  coverage <- depth_at_locus(coverage, coverage_files, filepath, prefixes)
  stats_matrix <- make_stats_matrix(coverage)
  coverage <- get_samples_below_20(coverage)
  coverage <- add_failure_at_locus(coverage, prefixes)
  write_coverage(coverage, stats_matrix, filepath)
}


#main()
