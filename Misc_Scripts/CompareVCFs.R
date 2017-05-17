BasicFilepath <- "/home/joe/Analysis/Ongoing/ValidationLibraryTesting/"
setwd(BasicFilepath)
GenePanel <- "10"
FilepathNew <- paste0("New/", GenePanel, "GenePanel/", GenePanel, "Panel_out/")
FilepathOld <- paste0("Old/", GenePanel, "GenePanel/")

#Get new filenames
NewNames <- Sys.glob(paste0(FilepathNew, "*multianno.vcf"))
OldNames <- Sys.glob(paste0(FilepathOld, "*multianno.vcf"))
Prefixes <- unlist(strsplit(NewNames, "/"))[seq(from=4,
                                                   by=4,
                                                   to=(4*length(NewNames)))]
Prefixes <- unlist(strsplit(Prefixes, ".myanno.hg38_multianno.vcf"))

nonMatching <- function(Old, New){
  OldPos <- apply(Old, 1, function(x) paste0(x[1], ",", x[2]))
  NewPos <- apply(New, 1, function(x) paste0(x[1], ",", x[2]))
  OldOnly <- Old[which(!(OldPos %in% NewPos)),]
  NewOnly <- New[which(!(NewPos %in% OldPos)),]
  return(list(OldOnly, NewOnly))
}

for(sample in 1:length(NewNames)){
  print(NewNames[sample])
  print(OldNames[sample])
  New <- read.table(NewNames[sample], quote="\"", stringsAsFactors=FALSE)
  Old <- read.table(OldNames[sample], quote="\"", stringsAsFactors=FALSE)
  Results <- nonMatching(Old, New)
  OldResults <- Results[[1]]
  NewResults <- Results[[2]]
  write.table(OldResults,
              file=paste0(Prefixes[sample], "OldOnly.vcf"),
              row.names=FALSE,
              quote=FALSE)
  write.table(NewResults,
              file=paste0(Prefixes[sample], "NewOnly.vcf"),
              row.names=FALSE,
              quote=FALSE)
}
