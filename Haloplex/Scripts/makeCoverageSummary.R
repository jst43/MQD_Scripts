sample_Line <- commandArgs(trailingOnly=TRUE)
print(paste0("Biological Sample Prefix is ", sample_Line))

#Read in Coverage table
print("Reading in Coverage File")
Coverage <- read.table("/home/joe/MQD_Scripts/Haloplex/Dependent_Files/Coverage_IntronsExons.csv",
                       header=TRUE,
                       quote="\"",
                       sep=",",
                       stringsAsFactors=FALSE)

filenames <- Sys.glob(paste0(sample_Line,"*coverage"))

#For each file, read in the total depth column and add to Coverage
print("Importing coverage data")
for(Coveragefile in filenames){
  Prefix <- strsplit(Coveragefile, split=".c")[[1]][1]
  Coverage <- cbind(Coverage, read.csv(file=Coveragefile,
                                       sep="\t",
                                       stringsAsFactors=FALSE,
                                       colClasses=c("NULL", "NULL", "NULL", NA)))
  colnames(Coverage) <- c(colnames(Coverage)[1:(ncol(Coverage)-1)], Prefix)
}
rm(Coveragefile, Prefix)

#Remove all rows in Coverage which do not correspond to Exons
Coverage <- Coverage[which(grepl("Exon", Coverage[,1])),]
rownames(Coverage) <- NULL

#Add columns for Average Depth at locus, and Total Depth at locus
if(length(filenames)>1){
  Coverage <- cbind(Coverage[,1:2],
                    Total_Depth=rowSums(Coverage[,3:ncol(Coverage)]),
                    Ave_Depth=round(rowMeans(Coverage[,3:ncol(Coverage)]), digits=2),
                    Coverage[,3:ncol(Coverage)])
} else {
  Coverage <- cbind(Coverage[,1:2],
                    Total_Depth=Coverage[,3],
                    Ave_Depth=Coverage[,3],
                    Coverage[,3])
  colnames(Coverage) <- c(colnames(Coverage)[1:(ncol(Coverage)-1)], strsplit(filenames, split=".c")[[1]][1])
}


#Function to determine the number of loci for each sample with depth >= threshold
nLociOverThreshold <- function(threshold){
  sumThreshold <- numeric(length=(ncol(Coverage)-3))
  for(i in 1:length(sumThreshold)){
    sumThreshold[i] <- sum(Coverage[,(i+3)]>=threshold)
  }
  return(sumThreshold)
}

#Generate matrix to fit underneath Coverage in the csv file, giving statistics for each
#sample
print("Generating Statistics")
nLoci_10 <- nLociOverThreshold(10)
nLoci_20 <- nLociOverThreshold(20)
nLoci_50 <- nLociOverThreshold(50)
percentLoci_10 <- round((nLoci_10/nrow(Coverage))*100, digits=2)
percentLoci_20 <- round((nLoci_20/nrow(Coverage))*100, digits=2)
percentLoci_50 <- round((nLoci_50/nrow(Coverage))*100, digits=2)
sample_ave_depth <- round(unname(colMeans(Coverage[,4:ncol(Coverage)])), digits=2)
sample_max_depth <- numeric(length=(ncol(Coverage)-3))
for(i in 1:length(sample_max_depth)){
  sample_max_depth[i] <- max(Coverage[,(i+3)])
}
statsMatrix <- rbind(nLoci_10, percentLoci_10,
                     nLoci_20, percentLoci_20,
                     nLoci_50, percentLoci_50,
                     sample_ave_depth,
                     sample_max_depth)
rm(nLoci_10,
   nLoci_20, nLoci_50,
   percentLoci_10, percentLoci_20,
   percentLoci_50, sample_ave_depth,
   sample_max_depth)
statsMatrix <- cbind(c("# > 10",
                           "% > 10",
                           "# > 20",
                           "% > 20",
                           "# > 50",
                           "% > 50",
                           "Ave. Depth",
                           "Max Depth"), statsMatrix)

#Add column to Coverage which gives the number of samples with depth < 20 at each locus
samples_below_20 <- numeric(length=nrow(Coverage))
for(i in 1:length(samples_below_20)){
  samples_below_20[i] <- sum(Coverage[i, 5:ncol(Coverage)]<20)
}
Coverage <- cbind(Coverage, samples_below_20)
rm(samples_below_20)

#Generate columns for Coverage which show how many duplicates of a sample cover that locus
filenames <- filenames[grepl("dup", tolower(filenames))]
filenames <- unique(unlist(strsplit(filenames, "_"))[seq(1,(2*length(filenames)),2)])
generateSuffix <- function(filePrefix){
  indices <- grep(filePrefix, colnames(Coverage))
  fullNames <- colnames(Coverage)[indices]
  suffixes <- unlist(strsplit(fullNames, paste0(filePrefix,"_")))
  suffixes <- suffixes[seq(2, length(suffixes), 2)]
  return(suffixes)
}
createSampleFailvec <- function(filePrefix){
  suffixes <- generateSuffix(filePrefix)
  indices <- grep(filePrefix, colnames(Coverage))
  tempFrame <- Coverage[,c(1,indices)]
  failVec <- character(length=nrow(Coverage))
  for(i in 1:nrow(tempFrame)){
    binaryVec <- tempFrame[i,2:ncol(tempFrame)]>=20
    if(sum(binaryVec)==length(suffixes)){
      failVec[i] <- "Covered"
    }
    else if(sum(binaryVec)==0){
      failVec[i] <- "All Failed"
    }
    else{
      failVec[i] <- paste0(paste(suffixes[!binaryVec], collapse=" and "), " Failed")
    }
  }
  return(failVec)
}

print("Adding failure at Locus data")
for(filePrefix in filenames){
  print(filePrefix)
  failVec <- createSampleFailvec(filePrefix)
  CoverageNames <- c(colnames(Coverage), filePrefix)
  Coverage <- cbind(Coverage, failVec)
  colnames(Coverage) <- CoverageNames
}

#Write Coverage to file
print("Writing to file")
write.csv(Coverage, file=paste0(sample_Line,"_Coverage.csv"), row.names=FALSE, quote=FALSE)

#Append statsMatrix
for(i in 1:nrow(statsMatrix)){
  line <- paste0(",,",paste(statsMatrix[i,], collapse=","))
  write(line, file=paste0(sample_Line,"_Coverage.csv"), append=TRUE)
}
print("Done")
