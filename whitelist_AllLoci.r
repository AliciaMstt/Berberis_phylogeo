rm(list = ls())

# Define WD and other directories
#WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
WD<-"~/BerL_1_2_3/3Berberis_phylogeo"
setwd(WD) 
list.files()

outfolder = "/data.out/PopSamples_m3/" # tvsdirectory, where the tsv SNP matrices are


### load samples meta information (lane, barcode, pop, etc) and list of potetial paralog loci as blacklist
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

## Load custom functions to create white list of loci and population maps:
source(paste0(WD,"/bin/whiteRADlist.R"))
source(paste0(WD,"/bin/whitePopMap.R"))

#Define desired population order
PopOrder <- c(1,2,3,4,9,5,6,7,8)

#### BerAll
# Run whiteRAD list for m3
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"), 
  writedirectory= paste0(WD,outfolder, "AllLoci/"))
  
# Create population map ALL samples, EXCLUDING replicates
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"), 
  writedirectory= paste0(WD,outfolder, "AllLoci/"), matinfo=matinfo, PopOrder=PopOrder)
# Rename output file to include BerAll name
system("mv data.out/PopSamples_m3/AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep.tsv data.out/PopSamples_m3/AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_BerAll.tsv")



# Create population map ALL samples, INCLUDING replicates
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"), 
  writedirectory= paste0(WD,outfolder, "AllLoci/"), matinfo=matinfo, PopOrder=PopOrder)
# Rename output file to include BerAll name
system("mv data.out/PopSamples_m3/AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep.tsv data.out/PopSamples_m3/AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_BerAll.tsv")

#### BerwoOut
# Create population map w/o and Out samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Out
s<-grep("Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_BerwoOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o and Out samples, INCLUDING replicates
x<-read.delim(file= paste0(WD, outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD, outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_BerwoOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)


#### woZaOut
# Create population map w/o Za and Out samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_woZaOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o Za and Out samples, INCLUDING replicates
x<-read.delim(file= paste0(WD, outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD, outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_woZaOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

#### BerSS
# Create population map w/o Za, Out and An samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out|An", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_BerSS.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o Za and Out and An samples, INCLUDING replicates
x<-read.delim(file= paste0(WD, outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out|An", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD, outfolder, "AllLoci/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_BerSS.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)
