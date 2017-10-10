#Script to process VAF files
#Packages
require(readr)
require(reshape2)

#Constants
TargetGene_Path <- "/home/joe/MQD_Scripts/Haloplex/Dependent_Files/Haloplex_targetGenes.csv"

#Commandline variables
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=1){
  stop("FilterVCF.R only accepts one input")
}
if(args[1]=="h" | args[1]=="help"){
  stop("Usage: Rscript FilterVCF.R inputVCF.txt")
}
VCF_Path <- args[1]
print(paste("VCF Path is", VCF_Path))

#Code
VCF_Name <- tail(strsplit(VCF_Path, "/")[[1]], n=1)
Pipeline <- tail(strsplit(VCF_Name, "_")[[1]], n=2)[1]

#Read in file
VCF <- read_delim(VCF_Path,
                  "\t", escape_double = FALSE,
                  col_types = cols(Chr = col_character(),
                                   A3 = col_character()),
                  trim_ws = TRUE)


#Generate Data
print("Generating Column Data")
alt_allele <- colsplit(VCF$A7, ",", names=c("Alt1", "Alt2"))
VCF <- VCF[,c(seq(1,49), 59)]
if(Pipeline=="SNV"){
  allele_info <- colsplit(VCF$A12, ":", names=c("GT", "AD", "DP", "GQ", "PL"))[,seq(1,3)]
  AD_Info <- colsplit(allele_info$AD, ",", names=c("AD.Ref", "AD.Alt1", "AD.Alt2"))
}else if(Pipeline=="Pindel"){
  allele_info <- colsplit(VCF$A12, ":", names=c("GT", "AD"))
  AD_Info <- colsplit(allele_info$AD, ",", names=c("AD.Ref", "AD.Alt1", "AD.Alt2"))
  allele_info <- data.frame(allele_info, DP=rowSums(AD_Info, na.rm=TRUE))
}
VAF <- data.frame(VAF.Ref=(AD_Info$AD.Ref/allele_info$DP), VAF.Alt1=(AD_Info$AD.Alt1/allele_info$DP), VAF.Alt2=(AD_Info$AD.Alt2/allele_info$DP))
Variant_ID <- paste0("Chr", VCF$Chr, ":", VCF$Start, "_", VCF$Ref, ">", VCF$Alt)

#Create Data Frame
print("Creating VCF")
VCF <- data.frame(Sample=VCF$Sample,
           Variant.ID=Variant_ID,
           Frequency=character(length=nrow(VCF)),
           Check.Bam.File=character(length=nrow(VCF)),
           VCF[,8],
           AD_Info,
           DP=allele_info$DP,
           VAF,
           VCF[,c(7,10,11)],
           Reproducible=character(length=nrow(VCF)),
           VCF[,seq(14,46)],
           VCF[,c(2,5)],
           alt_allele,
           VCF[,c(12,3,4)])

rm(AD_Info, allele_info, alt_allele, VAF, Variant_ID)

#Filtering
print("Filtering VCF")
Haloplex_targetGenes <- read_csv(TargetGene_Path, col_names = FALSE)
VCF <- VCF[which(VCF$Gene.refGene %in% Haloplex_targetGenes$X1),]
VCF <- subset(VCF, ExonicFunc.refGene!="synonymous SNV")
VCF <- subset(VCF, X1000g2015aug_eur=="." | as.numeric(X1000g2015aug_eur)<=0.01 | (as.numeric(X1000g2015aug_eur)>0.01 & cosmic70!="."))
VCF <- subset(VCF, Func.refGene=="exonic" | Func.refGene=="exonic;splicing" | Func.refGene=="splicing")
if(Pipeline=="Pindel"){
  VCF <- subset(VCF, AD.Alt1>=5 | AD.Alt2>=5)
  VCF <- subset(VCF, VAF.Alt1>=0.015 | VAF.Alt2>=0.015)
}

#Fill Frequency and Reproducible feeds
VCF$Frequency <- unlist(lapply(VCF$Variant.ID, function(x) sum(VCF$Variant.ID %in% x)))
#Generate reproducible counts
VCFDupVec <- paste0(gsub('.{1}$', '', VCF$Sample), VCF$Variant.ID)
VCF$Reproducible <- unlist(lapply(VCFDupVec, function(x) sum(VCFDupVec %in% x)>1))
VCF$Reproducible[!grepl("dup", tolower(VCF$Sample))] <- NA

#Write to file
print("Writing VCF")
write.table(VCF, file=VCF_Name, quote=FALSE, sep="\t", row.names=FALSE)
