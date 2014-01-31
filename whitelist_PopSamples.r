rm(list = ls())

# Define WD and other directories
WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
setwd(WD) 
list.files()

outfolder = "/data.out/" # tvsdirectory, where the SNP matrices are

### load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

## Load custom functions to create white list of loci and population maps:
source(paste0(WD,"/bin/whiteRADlist.R"))
source(paste0(WD,"/bin/whitePopMap.R"))

#Define desired population order
PopOrder <- c(1,2,3,4,9,5,6,7,8)

# Run whiteRAD list for m3
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))

# Create population map ALL samples, EXCLUDING replicates
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"), matinfo=matinfo, PopOrder=PopOrder)

# Create population map ALL samples, INCLUDING replicates
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"), matinfo=matinfo, PopOrder=PopOrder)



# Create population map w/o Za and Out samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,"/docs/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,"/docs/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep_woZaOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o Za and Out ALL samples, INCLUDING replicates
x<-read.delim(x, file= paste0(WD,"/docs/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,"/docs/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep_woZaOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)
