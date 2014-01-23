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

# Run whiteRAD lsit for m3
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))

# Create population map EXCLUDING replicates
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"), matinfo=matinfo, PopOrder=PopOrder)

# Create population map INCLUDING replicates
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"), matinfo=matinfo, PopOrder=PopOrder)


