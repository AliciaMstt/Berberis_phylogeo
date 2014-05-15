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
  P0.5<-getP0.5(popsumstats=popsumstats, minp=2, maxp=9)
  pl<-levels(as.factor(P0.5$Locus.ID))

  # Save list of potential paralogs
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


######## See distribution of loci P=0.5 among populations and species
library(ggplot2)

##### Loci at P=0.5 by Population
# P=0.5 in at least one pop
P0.5atleast1<-getP0.5(popsumstats=popsumstats, minp=1, maxp=9) 
subdf<-P0.5atleast1[P0.5atleast1$P==0.5, ] #subset dataframe to loci P=0.5
subdf<-subdf[!duplicated(subdf[, c("Pop.Name", "Locus.ID")]), ] # there could be >1 SNPs at P=.5 per RAD-locus, keep only the first
plt<- ggplot(subdf, aes(x=Pop.Name)) + theme_bw() # plot by pop
plt<- plt + geom_bar(fill="white", colour="black")
plt<- plt + scale_x_discrete(limits=c("Aj","Iz","Ma","Pe","Tl","To","An", "Za", "Out")) 
plt<- plt + ylab("Number of RAD-loci") + xlab("Population")
plt 
plt1<-plt

#### Loci at P=0.5 UNIQUE by population
# P=0.5 in only one population
P0.5only1<-getP0.5(popsumstats=popsumstats, minp=1, maxp=1) 
# Plot number of loci with at least one SNP at P=0.5 in only one population
subdf<-P0.5only1[P0.5only1$P==0.5, ] #subset dataframe to loci P=0.5
subdf<-subdf[!duplicated(subdf[, c("Pop.Name", "Locus.ID")]), ] # there could be >1 SNPs at P=.5 per RAD-locus, keep only the first
plt<- ggplot(subdf, aes(x=Pop.Name)) + theme_bw() # plot by pop
plt<- plt + geom_bar(fill="white", colour="black")
plt<- plt + scale_x_discrete(limits=c("Aj","Iz","Ma","Pe","Tl","To","An", "Za", "Out")) 
plt<- plt + ylab("Number of RAD-loci") + xlab("Population")
plt 


#### Loci at P=0.5 shared with B.alpina sensu stricto
# defina B. alpina sensu stricto
ingroup=c("Aj","Iz","Ma","Pe","Tl","To")
# keep only P=0.5 of the P=0.5 loci in at least one Poo
P0.5<-P0.5atleast1$P==0.5
P0.5<-P0.5atleast1[P0.5, ]
# Split data to have managable sizes
dfs<-split(P0.5,P0.5$Locus.ID)
# select locus+bp (id) that are P=0.5 in at least 2 populations
desired<- lapply(dfs, function(x) ddply(x, .(id), function(x)
  if(sum(x$P==0.5)>=2)x) )
# select locus+bp (id) that are P=0.5 in at least one pop of the ingroup
desired<- lapply(desired, function(x) ddply(x, .(id), function(x)
  if(any(ingroup %in% x$Pop.Name))x) )
#put in a single df again
P0.5shared<- do.call("rbind", desired)
# Plot number of loci with at least one SNP at P=0.5 that are shared with B. alpina ingroup 
subdf<-P0.5shared
subdf<-subdf[!duplicated(subdf[, c("Pop.Name", "Locus.ID")]), ] # there could be >1 SNPs at P=.5 per RAD-locus, keep only the first
plt<- ggplot(subdf, aes(x=Pop.Name)) + theme_bw() # plot by pop
plt<- plt + geom_bar(fill="white", colour="black")
plt<- plt + scale_x_discrete(limits=c("Aj","Iz","Ma","Pe","Tl","To","An", "Za", "Out")) 
plt<- plt + ylab("Number of RAD-loci") + xlab("Population")
plt 


##### Stack in a single plot:
## Generate datasets
# P=0.5 in at least one pop
P0.5atleast1<-getP0.5(popsumstats=popsumstats, minp=1, maxp=9) 
P0.5atleast1<-P0.5atleast1[P0.5atleast1$P==0.5, ] #subset dataframe to loci P=0.5
P0.5atleast1<-P0.5atleast1[!duplicated(P0.5atleast1[, c("Pop.Name", "Locus.ID")]), ] # there could be >1 SNPs at P=.5 per RAD-locus, keep only the first

# P=0.5 UNIQUE to one population
P0.5only1<-getP0.5(popsumstats=popsumstats, minp=1, maxp=1) 
P0.5only1<-P0.5only1[P0.5only1$P==0.5, ] #subset dataframe to loci P=0.5
P0.5only1<-P0.5only1[!duplicated(P0.5only1[, c("Pop.Name", "Locus.ID")]), ] # there could 

## Plot them with doging bars
plt<- ggplot(P0.5atleast1, aes(x=Pop.Name))
plt<- geom_bar()
plt<- plt + geom_bar(data=P0.5only1)
plt


##### See distribution of potential paralogs among spp and pops
require(reshape2)

# Keep only pop and Loci data
df<-P0.5[,c(1,3)]

## Transform to presence/absence matrix of Pop x loci
# cast to wide format
mt<-dcast(df, Pop.Name ~ Locus.ID)
rownames(mt) <- mt[,1] #change first col to rownames
mt<- mt[,2:ncol(mt)]
# Some numbers are >1 because there was more than one SNP by loci, change them to 1
mt[mt > 1] = 1
mt<-as.matrix(mt)


## Build distance matrix 
par(mfrow = c(1, 1))
dist_mt<-dist(mt, method="binary") 


# plot distance matrix
table.dist(dist_mt, csize=.7, labels=labels(dist_mt))
table.dist(lower.tri(dist_mt), csize=.7, labels=labels(dist_mt))

## NJ tree
# Build tree from distance matrix
tree<-nj(dist_mt)
# root
tree<-root(tree, outgroup="Out")
#plot
plot(tree)


#### What happens with sets of 831 random loci?
# select random loci
for(i in 1:10){
set.seed(i)
locinames<-levels(as.factor(popsumstats$Locus.ID)) # extract loci names
randloci <- sample(locinames, 831) # select 831 loci randombly
df<-popsumstats[popsumstats$Locus.ID %in% randloci,]

# Keep only pop and Loci data
df<-df[, c(1,3)]

## Transform to presence/absence matrix of Pop x loci
# cast to wide format
mt<-dcast(df, Pop.Name ~ Locus.ID)
rownames(mt) <- mt[,1] #change first col to rownames
mt<- mt[,2:ncol(mt)]
# Some numbers are >1 because there was more than one SNP by loci, change them to 1
mt[mt > 1] = 1
mt<-as.matrix(mt)

## Build distance matrix 
par(mfrow = c(1, 1))
dist_mt<-dist(mt, method="binary") 

# plot distance matrix
table.dist(dist_mt, csize=.7, labels=labels(dist_mt))

## NJ tree
# Build tree from distance matrix
tree<-nj(dist_mt)
# root
tree<-root(tree, outgroup="Out")
#plot
plot(tree)
}




