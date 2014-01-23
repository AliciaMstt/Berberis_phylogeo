rm(list = ls())

# Define WD and other directories
WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
setwd(WD) 
list.files()

outfolder = "/data.out/" # tvsdirectory, where the SNP matrices are

### load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

##### Read Stacks population output file batch_1.sumstats_summary.tsv

source(paste0(WD,"/bin/read.sumstats_summary.R"))

read.sumstats_summary(file=paste0(WD,"/data.out/PopSamples_m3/Popsouts_RselecStructure/out.noreplicates/batch_1.sumstats_summary.tsv"),
                      allP=TRUE, npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  

