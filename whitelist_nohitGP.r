rm(list = ls())

# Define WD and other directories
#WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
WD<-"~/BerL_1_2_3/3Berberis_phylogeo"
setwd(WD) 
list.files()

outfolder = "/data.out/PopSamples_m3/" # tvsdirectory, where the tsv SNP matrices are

# blast outputfile on which results are based
blastoutput = paste0(WD,"/blast/PopSamples_BeralpBt_m3_blastGP.out") 

#name
bloutname<- basename(blastoutput)
bloutname<- sub(".out", "", bloutname)

### load samples meta information (lane, barcode, pop, etc) and list of potetial paralog loci as blacklist
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)


## Load custom functions to create population maps:
source(paste0(WD,"/bin/whitePopMap.R"))

#Define desired population order
PopOrder <- c(1,2,3,4,9,5,6,7,8)

################ Generate Whitelist_HitPlants with blasted loci

#### Get blast output data
hitGP<-read.delim(file=blastoutput,
  header=FALSE, stringsAsFactors=FALSE)

## Add header
colnames(hitGP)<-c("qacc","sacc", "evalue", "bitscore","qcovs", "length", "pident", "staxids", "sscinames", "stitle")
# qacc means Query accesion
# evalue means Expect value
# bitscore means Bit score (normalized raw score that considers the search space)
# qcovs means Query Coverage Per Subject
# length means Alignment length
# pident means Percentage of identical matches
# staxids means unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
# sscinames means unique Subject Scientific Name(s), separated by a ';'
# stitle means Subject Title


## Generate whitelist of NohitGP loci
# get list of all loci used for blast
final = read.delim(paste(file=paste0(WD, outfolder, "PopSamples_BeralpBt_m3.SNP.SNPs"), sep = ""), header = T) 
# substract the loci that blasted from the loci used
"%w/o%" <- function(x, y) x[!x %in% y]
NohitGP <- "%w/o%"(final$CatalogID, hitGP$qacc)


## How many loci did not blast against green plants?
length(NohitGP)

#### Save whitelist of loci 
# NohitGP
write(NohitGP, file= paste0(WD,outfolder,"NohitGP/",  bloutname, "NohitGP_whitelist.tsv"), ncolumns = 1)



############## Generate Population Maps

#### BerAll

# Create population map ALL samples, EXCLUDING replicates
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"), 
  writedirectory= paste0(WD,outfolder, "NohitGP/"), matinfo=matinfo, PopOrder=PopOrder)
# Rename output file to include BerAll name
system(paste0("mv data.out/PopSamples_m3/NohitGP/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_norep.tsv data.out/PopSamples_m3/NohitGP/",bloutname, "_PopMap_norep_BerAll.tsv"))

# Create population map ALL samples, INCLUDING replicates
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_BeralpBt_m3.SNP.SNPs"), 
  writedirectory= paste0(WD,outfolder, "NohitGP/"), matinfo=matinfo, PopOrder=PopOrder)
# Rename output file to include BerAll name
system(paste0("mv data.out/PopSamples_m3/NohitGP/PopSamples_BeralpBt_m3.SNP.SNPs_PopMap_withrep.tsv data.out/PopSamples_m3/NohitGP/",bloutname, "_PopMap_withrep_BerAll.tsv"))

#### BerwoOut
# Create population map w/o Za and Out samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_norep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_norep_BerwoOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o Za and Out samples, INCLUDING replicates
x<-read.delim(file= paste0(WD, outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_withrep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD, outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_withrep_BerwoOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)


#### woZaOut
# Create population map w/o Za and Out samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_norep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_norep_woZaOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o Za and Out samples, INCLUDING replicates
x<-read.delim(file= paste0(WD, outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_withrep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD, outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_withrep_woZaOut.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

#### BerSS
# Create population map w/o Za, Out and An samples, EXCLUDING replicates
x<-read.delim(file= paste0(WD,outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_norep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out|An", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD,outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_norep_BerSS.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)

# Create population map w/o Za and Out and An samples, INCLUDING replicates
x<-read.delim(file= paste0(WD, outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_withrep_BerAll.tsv"), header=FALSE)
# look for sample names that do not include Za or Out
s<-grep("Za|Out|An", x[,1], invert= TRUE)
#save those samples
x<-x[s,]
write.table(x, file= paste0(WD, outfolder, "NohitGP/PopSamples_BeralpBt_m3_blastGP_PopMap_withrep_BerSS.tsv"),
  sep = "\t",
  row.names =FALSE, col.names=FALSE,
  quote=FALSE)
