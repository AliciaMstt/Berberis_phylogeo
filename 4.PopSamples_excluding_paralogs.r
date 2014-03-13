# This script uses Stacks populations summary stats output to identify potential paralog loci and create a whitelist excluding them
rm(list = ls())
WD<-"/gpfs/bio/bxh10nyu/BerL_1_2_3/PopSamples/" # to work in the cluster
setwd(WD) 
list.files()

##### Extract paralogs 
#     defined as those loci are P=0.5 in a given minimal number of populations

# Get functions
source(paste0(WD, "/bin/read.sumstats.R"))
source(paste0(WD, "/bin/getpparal.R"))
library(plyr)

## Loop for every dataset
for (i in c("m3", "m4", "m10", "def")){
  popsumstats <-read.sumstats(file=paste0(WD,"/data.out/PopSamples_",i,"/Popsouts_Rselec/out.noreplicates/batch_1.sumstats.tsv"),
    npop=9, popNames=c("Aj","An","Iz","Ma","Out","Pe","Tl","To","Za"))                  
  #get all loci names  
  allloci<-levels(as.factor(popsumstats$Locus.ID))
  
  # get paralogs
  pl<-getpparal(popsumstats, mp=6)
  
  # get no paralogs
  "%w/o%" <- function(x, y) x[!x %in% y]
  noparalogs<-allloci %w/o% pl
  
  # Save white list
  writedirectory= paste0(WD,"/docs")
  x<-grep(pattern=i, x=list.files(writedirectory), value=TRUE)
  x<-grep(pattern="whitelist.tsv", x=x, value=TRUE)  
  name<-gsub(pattern="whitelist.tsv", replacement="whitelistnoparlgs.tsv", x)
  write(noparalogs, file= paste0(writedirectory, "/", name), ncolumns = 1)
}