
# Generate fasta file from .SNPs.SNPs matrix produced by postcleaning
rm(list = ls())
WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
setwd(WD) 
list.files()


# Call function
source("bin/tsv2Fasta.R")

# Run funtion with data
tsv2Fasta(tsvfile=paste0(WD,"/data.out/PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), 
          outputfile=paste0(WD,"/data.out/PopSamples_m3/PopSamples_BeralpBt_m3.FASTA"), len=83)
