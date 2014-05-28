rm(list = ls())

# Define WD and other directories
#WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
WD<-"~/BerL_1_2_3/3Berberis_phylogeo"
setwd(WD) 
list.files()

outfolder = "/data.out/PopSamples_m3/" # tvsdirectory, where the tsv SNP matrices are


### load samples meta information (lane, barcode, pop, etc) and list of potetial paralog loci as blacklist
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

### Soruce functions
source(paste0(WD,"/bin/PopMapEQsz.R")) # to create pop map with homogenized sampling size
source(paste0(WD,"/bin/whiteRADlist.R")) # to create whitelist of RAD loci

### Run whiteRAD list for m3
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"), 
  writedirectory= paste0(WD,outfolder, "AllLoci/"))

### Create population map
# Define custom population order
PopOrder <- c(1,2,3,4,9,5,6,7,8)
PopMapEQsz(tsv=paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"),
           writedirectory= paste0(WD,outfolder, "AllLoci/"),
           matinfo=matinfo, homogp=c(1:8), dsz=4, PopOrder=PopOrder)


