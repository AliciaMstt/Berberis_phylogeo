# This script uses Stacks populations summary stats output to identify potential paralog loci and create a whitelist excluding them
rm(list = ls())
WD <- "/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo/"
#WD<-"/gpfs/bio/bxh10nyu/BerL_1_2_3/PopSamples/" # to work in the cluster
setwd(WD) 
list.files()

datadirectory<-paste0(WD,"/data.out/PopSamples_m3")


##### Extract paralogs 
#defined as those loci are P=0.5 in a given minimal number of populations or species 
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
  P0.5<-getP0.5(popsumstats=popsumstats, minp=2, maxp=9)
  pl<-levels(as.factor(P0.5$Locus.ID))

  # Save list of potential paralogs
  writedirectory= paste0(WD,"/docs")
  write(pl, file= paste0(writedirectory, "/", "potentialparalogs"), ncolumns = 1)

### Examine Hobs and FIS of extracted loci
# Get subset of loci
df<-popsumstats[popsumstats$Locus.ID %in% pl,]
  # in B. alpina ingroup pop.
  ingroup=c("Aj","Iz","Ma","Pe","Tl","To")
  # keep only ingroup
  df<-popsumstats$Pop.Name %in% ingroup
  df<-popsumstats[df,]
  # Split data to have managable sizes
  dfs<-split(df,df$Locus.ID)
  # Create column id with locus+BP info 
  dfs<-lapply(dfs, function(x) transform(x,id =interaction(Locus.ID,BP)))
  # select locus+bp (id) that are P=0.5 in at least mixp populations and no more than maxp populations
  desired<- lapply(dfs, function(x) ddply(x, .(id), function(x)
    if(sum(x$P==0.5)>=2)x) )
  # put in a single df again
  x<- do.call("rbind", desired)

# observed Het
obsHet1<-x$Obs.Het==1 #get subset where observed Heterosigosity = 1.
# % p=0.5 loci where all individuals are heterozygous (obsHe=1)
sum(obsHet1) * 100 / nrow(x)
# FIS
negFIS<-x$Fis<1 #get subset where FIS is negative
# % p=0.5 loci where FIS is negative
sum(negFIS) * 100 / nrow(x)
# % of negFIS loci with < 0.5 
negFIS<-x[negFIS,]
negFIS.5<-negFIS$Fis<=-0.5
sum(negFIS.5) * 100 / nrow(negFIS)


##### Extract loci were P=0.5 in any popupation
P0.5<-getP0.5(popsumstats=popsumstats, minp=1, maxp=9)
pl<-levels(as.factor(P0.5$Locus.ID))

# Save list of loci were p=0.5 
writedirectory= paste0(WD,"/docs")
write(pl, file= paste0(writedirectory, "/", "lociP05"), ncolumns = 1)


### Examine Hobs and FIS of extracted loci
# Get subset of loci
  df<-popsumstats[popsumstats$Locus.ID %in% pl,]
  # get SNP-loci where p=0.5
  P0.5<-df$P==0.5 
  # how many?
  sum(P0.5)
  x<-df[P0.5,] #subset
# observed Het
obsHet1<-x$Obs.Het==1 #get subset where observed Heterosigosity = 1.
# % p=0.5 loci where all individuals are heterozygous (obsHe=1)
sum(obsHet1) * 100 / nrow(x)
# FIS
negFIS<-x$Fis<1 #get subset where FIS is negative
# % p=0.5 loci where FIS is negative
sum(negFIS) * 100 / nrow(x)
# % of negFIS loci with < 0.5 
negFIS<-x[negFIS,]
negFIS.5<-negFIS$Fis<=-0.5
sum(negFIS.5) * 100 / nrow(negFIS)


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
liSNPs$ind.names

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










