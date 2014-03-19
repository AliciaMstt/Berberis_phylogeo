# This script uses Stacks populations summary stats output to identify potential paralog loci and create a whitelist excluding them
rm(list = ls())
WD <- "/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo/"
#WD<-"/gpfs/bio/bxh10nyu/BerL_1_2_3/PopSamples/" # to work in the cluster
setwd(WD) 
list.files()

datadirectory<-paste0(WD,"/data.out/PopSamples_m3")


##### Extract paralogs 
#     defined as those loci are P=0.5 in a given minimal number of populations


  ## Get data
source(paste0(WD, "/bin/read.sumstats.R"))
popsumstats <-read.sumstats(file=paste0(datadirectory,"/Popsouts_Rselec/out.noreplicates/batch_1.sumstats.tsv"),
    npop=9, popNames=c("Aj","An","Iz","Ma","Out","Pe","Tl","To","Za"))                  
  #get all loci names  
  allloci<-levels(as.factor(popsumstats$Locus.ID))
  
  ## Get potential paralogs
  ## Function to subset the potential paralog loci from the batch_1.sumstats.tsv output of Stacks 
  ## paralogs are defined as those loci are P=0.5 in at least mp populations. 
  source(paste0(WD, "/bin/getP0.5.R"))
  P0.5<-getP0.5(popsumstats=popsumstats, mp=2)
  pl<-levels(as.factor(P0.5$Locus.ID))

  # Save list of potential paralots
  writedirectory= paste0(WD,"/docs")
  write(pl, file= paste0(writedirectory, "/", "potentialparalogs"), ncolumns = 1)

###### Plot only potential paralogs to see if they are all heterozygous
# Inport using the plink file exported from populations Stacks program
# using whilelist file (i.e. loci in final matrix) and pop map excluding replicates
require(adegenet)
plinkfile=paste0(datadirectory, "/Popsouts_Rselec/out.noreplicates/plink.raw")
liSNPs<- read.PLINK(file=plinkfile)

## Check number of recovered samples and loci
print("number of recovered samples and total number of SNPs")
print(nInd(liSNPs)) 
print(nLoc(liSNPs))

## Visualize the Matrix as a plot
par(mfrow = c(1, 1))
glPlot(liSNPs)
title("All loci")

# Subset paralog loci
lociall<-locNames(liSNPs, withAlleles=FALSE) #get loci 
lociall<-gsub("_._.|_.._.", replacement="", lociall) #take out the BP position and base to keep only loci names
plplink<- lociall %in% pl # get the ones that match against the paralog list
paliSNPs<-liSNPs[, plplink] # keep all individuals, and loci in plplink
  
#Check number of recovered samples and loci
print("number of recovered samples and total number of SNPs")
print(nInd(paliSNPs)) 
print(nLoc(paliSNPs))

# Visualize the Matrix as a plot
par(mfrow = c(1, 1))
glPlot(paliSNPs)
title("Putatively paralogous SNP loci")


# Plot both
par(mfrow = c(2, 1))
glPlot(liSNPs)
title("All SNP loci")
glPlot(paliSNPs)
title("Putatively paralogous SNP loci")





