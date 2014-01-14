# This script uses the PostCleaning.r outputs to evaluate number of SNPs, RADS and perform a PCoA and other Pop.Structure analyses 

rm(list = ls())
WD<-"/Volumes/TO_GO/BerL_1_2_3/2R/PopSamples" 
setwd(WD) 
list.files()


#### --- Load necessary directories

workfolder = "/savedanalyses" #leave as it is
outfolder = "/data.out/" #leave as it is
infolder = "/data.in/" # specify data.in folder in script directory

### load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)


######## GENERAL TOOLS ----

### Get number of alleles per specimen
nall_cell = function(val){
  nall = NA
  if(is.na(val) == F){
    if(val != ""){
      nall = length(as.numeric(strsplit(val, "/")[[1]]))
    }
  } else {
    nall = NA
  }
  nall 
}
nall_row = function(vec) sapply(vec, nall_cell)

### Get average coverage per rad and per specimen
mcov_cell2 = function(val){
  nall = NA
  if(is.na(val) == F){
    if(val != ""){
      nall = mean(as.numeric(strsplit(as.character(val), "/")[[1]]))
    }
  } else {
    nall = NA
  }
  nall 
}
mcov_row2 = function(vec) sapply(vec, mcov_cell2)


#### Colors functions that will be used for plots
## to color a plot according to LANE group
# define list of nice colors for lanes
colors <-c("red", "green", "blue")
lane.names<-c("1","2","3")
# bind pops to colors to make a stable color key
cols.lane.key =cbind(lane.names, colors)
cols.lane.key
# function
lane.col<- function(x){ # x is the dataframe to plot
  plot.samples = names(x) #takes the samples names from the samples in the data we want to plot
  idx = match(plot.samples, matinfo$sample) # matches the samples names with the samples names of a matrix info
  info.plot = matinfo[idx,] #to add the information of the matinfo to the plot
  lanes = info.plot$Lane # the lane are now in the colum Lane, 
  lanes = as.factor(lanes)
  levels(lanes) = cols.lane.key[match(levels(lanes), cols.lane.key[, 1]), 2] # reassing the names of the lanes to its corresponding color key defined before
  lanes = as.character(lanes) #this is what is read as a color list for col=
}

## to color a plot according to POPULATION

#load RColorBrewer to have 12 nice contrasting colors per pop + grey for outgroups
library(RColorBrewer)
Mycols <-c(brewer.pal(5,name="Set1"), "brown", brewer.pal(7,name="Set2"), "grey")

#all sampling pops (more than Berberis)
all.pops<-c("Aj", "An", "Iz", "Pp", "Ma", "Pe", "Ne", "Ta", "Co", "Tl", "To", "Za", "Pt", "Out")

# bind pops to colors to make a stable color key
cols.pop.key =cbind(all.pops, Mycols)

#function
pops.col<- function(x){ # x is the dataframe to plot
  plot.samples = names(x) #takes the samples names from the samples in the data we want to plot
  idx = match(plot.samples, matinfo$sample) # matches the samples names with the samples names of a matrix info
  info.plot = matinfo[idx,] #to add the information of the matinfo to the plot
  pops = info.plot$Pop # the lane are now in the colum Lane, 
  pops = as.factor(pops)
  levels(pops) = cols.pop.key[match(levels(pops), cols.pop.key[, 1]), 2] # reassing the names of the lanes to its corresponding color key defined before
  pops = as.character(pops) #this is what is read as a color list for col=
}



#### Define functions to display PCoA

## PCoA display by POPULATION and defined colors
source(paste0(WD,"/bin/PCoA_pop.r"))

################################------ ANALYSES ------ #####################

######## Coverage, Information content, error rates and efficacy to detect structuring of genetic variation   ----

#### Function to estimate the no. of loci and SNPs as well 
# as the mean coverage, the loci, allele and SNP error rates

Summary_performance <- function(id, final, final.cov, plink){
#
  
# id is a chracter string to identify the dataset 
# final and final.cov are the paths to the SNP.SNPs and COV.COVs matrices procuded by PostCleaning.r
# plink is the path to the plink.raw file to be used to estimate SNP error rate
# final, final.cov and plink paths should be given from outfolder onwards
  
# open files
final = read.delim(paste(WD, outfolder, final, sep = ""), header = T) 
final.cov= read.delim(paste(WD, outfolder, final.cov, sep = ""), header = T) 

#### See number of loci obtained
  n.loci<-nrow(final)


#####check coverage in final selected samples per RADloci and specimen
  #use funtion  mcov_row2 previously defined (General tools). 
  ### estimate mean coverage (of the two alleles) or the coverage (if only one allele) for each RADloci
  m.finalcov <- lapply(final.cov[4:length(final.cov)], mcov_row2)
  m.finalcov <- data.frame(m.finalcov)
  

###estimate the mean and SD coverage per sample
  m.cov.sample <-sapply(m.finalcov, mean, na.rm=TRUE)
  m.cov.sample
  m.cov<-mean(m.cov.sample)
  sd.cov<-sd(m.cov.sample)
  
  ## merge with matinfo
  idx1 <- match(names(m.cov.sample), matinfo$sample) # matches the samples names with the samples names of a matrix info
  info1 <- matinfo[idx1,] #creates a subset witht he samples that matched
  m.cov.info <-cbind(info1, m.cov.sample) #cbind now that order and size is the same
  names(m.cov.info)[5]<-paste("mean.cov") #sensible name to m.cov.sample, "mean.cov" is the mean coverage of all RADloci of a sample
  m.cov.info
  
  
  #plot
  par(mfrow = c(1, 1))
  barplot(m.cov.info$mean.cov,
    col=lane.col(m.finalcov), #pops.col is the function defined before to match samples to its population and a specified color
    xlab="sample",
    ylab="mean coverage"
  ) 
  title(sub= "bars color according to lane")
  
  
  barplot(m.cov.info$mean.cov,
    col=pops.col(m.finalcov), #pops.col is the function defined before to match samples to its population and a specified color
    xlab="sample",
    ylab="mean coverage"
  )
  title(sub= "bars color according to population")
  
  
  ###Boxplot grouping by lane
  boxplot(m.cov.info$mean.cov ~ m.cov.info$Lane,
    border = cols.lane.key[, 2],
    xlab = "Lane",
    ylab = "Mean coverage")

## Estimate allele and loci error rates:
  # Load custom function
  source(paste0(WD,"/bin/LociAllele_error.R"))
  
  # Get the name of the parameter used 
  param<-id
  
  # Run funtion to get allele and loci error rates
  repliDiff<-LociAllele_error(mat=final, param=param)  
  repliDiff 
  
  # Estimate mean allele and loci error rates
  m.allele.error<-mean(as.numeric(repliDiff$allele.error.rate))
  m.locus.error<-mean(as.numeric(repliDiff$loci.error.rate))
  
  # Estimate SD from allele and loci error rates
  sd.allele.error<-sd(as.numeric(repliDiff$allele.error.rate))
  sd.locus.error<-sd(as.numeric(repliDiff$loci.error.rate))

## Estimate SNP error rate
  # load plink file
  liSNPs<- read.PLINK(file = paste0(WD,outfolder, plink))
  
  # Count number of SNPs (in plink the number of loci = number SNPs)
  n.SNPs<-liSNPs$n.loc
  
  # source SNP:error script
  source(paste0(WD,"/bin/SNPs_error.R"))
  
  # Run funtion to get SNP error rate
  y<-SNP_error(liSNPs=liSNPs, param=id)
  
  # Add SNP.error.rate to the matrix of results
  repliDiff <-merge(repliDiff, y)
  
  # Estimate mean SNP error rate
  m.SNP.error<-mean(as.numeric(paste(repliDiff$SNP.error.rate)))
  
  # Estimate SD of SNP error rate
  sd.SNP.error<-sd(as.numeric(paste(repliDiff$SNP.error.rate)))
  
## Boxplot error rates
  par(mfrow=c(1,3))
  boxplot(as.numeric(paste(repliDiff$loci.error.rate)) ~ repliDiff$parameter,
    las=2, main="Locus error rate")
  boxplot(as.numeric(paste(repliDiff$allele.error.rate)) ~ repliDiff$parameter,
    las=2, main="Allele error rate")
  boxplot(as.numeric(paste(repliDiff$SNP.error.rate)) ~ repliDiff$parameter,
    las=2, main="SNP error rate")
  par(mfrow=c(1,1))

#### Print replicates diferences
print("replicate differences for each pair")
print(repliDiff)  
  
#### Estimate mean of Fst distance matrix
  ## Import Fst pairwise distance matrix created by Stacks population algoritm in the full dataset excluding the replicates
  Fstmat<-data.matrix(read.delim(file = paste0(WD,outfolder,"PopSamples_", id,"/Popsouts_Rselec/out.noreplicates/batch_1.fst_summary.tsv"),
    row.names=1, fill=TRUE))  
  colnames(Fstmat)<- c("Aj","An","Iz","Ma","Out","Pe","Tl","To","Za") # add col names
  
  ## Estimate mean and SD Fst  
  meanFst.all<-mean(Fstmat[upper.tri(Fstmat)])
  sdFst.all<-sd(Fstmat[upper.tri(Fstmat)])
    
  #Estimate mean and SD Fst excluding Out and Za
  Fstmat.ingrp<-Fstmat[c(1:4,6:8),c(1:4,6:8)]
  meanFst.ingrp<-mean(Fstmat.ingrp[upper.tri(Fstmat.ingrp)])
  sdFst.ingrp<-sd(Fstmat.ingrp[upper.tri(Fstmat.ingrp)])
  
#### Return summary of results
summ.performance<- cbind(id, n.loci, n.SNPs, m.cov, sd.cov,
  m.locus.error, sd.locus.error,
  m.allele.error, sd.allele.error,
  m.SNP.error, sd.SNP.error, 
  meanFst.all, sdFst.all,
  meanFst.ingrp, sdFst.ingrp)
  
  
return(summ.performance)
  
}

## For m3
# run funtion
xm3<-Summary_performance(id= "m3",
                    final = "PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs",
                    final.cov = "PopSamples_m3/PopSamples_BeralpBt_m3.COV.COVs",
                    plink="PopSamples_m3/Popsouts_Rselec/out.replicates/plink.raw")

# Show results
xm3


## For m4
# run funtion
xm4<-Summary_performance(id= "m4",
  final = "PopSamples_m4/PopSamples_BeralpBt_m4.SNP.SNPs",
  final.cov = "PopSamples_m4/PopSamples_BeralpBt_m4.COV.COVs",
  plink="PopSamples_m4/Popsouts_Rselec/out.replicates/plink.raw")

# Show results
xm4

## For m10
# run funtion
xm10<-Summary_performance(id= "m10",
  final = "PopSamples_m10/PopSamples_BeralpBt_m10.SNP.SNPs",
  final.cov = "PopSamples_m10/PopSamples_BeralpBt_m10.COV.COVs",
  plink="PopSamples_m10/Popsouts_Rselec/out.replicates/plink.raw")

# Show results
xm10

## For def
# run funtion
xdef<-Summary_performance(id= "def",
  final = "PopSamples_def/PopSamples_BeralpBt_def.SNP.SNPs",
  final.cov = "PopSamples_def/PopSamples_BeralpBt_def.COV.COVs",
  plink<-"PopSamples_def/Popsouts_Rselec/out.replicates/plink.raw")

# Show results
xdef

### Show table with all resutls
summary_performance<-rbind(xm3, xm4, xm10, xdef)
summary_performance

######## SNPs frequency spectrum ----

#### Estimate SNP frequencies for each dataset
# Read plink files to genlight
require(adegenet)
liSNPs.m3<- read.PLINK(file = paste0(WD,outfolder,"PopSamples_m3/Popsouts_Rselec/out.noreplicates/plink.raw"))
liSNPs.m4<- read.PLINK(file = paste0(WD,outfolder,"PopSamples_m4/Popsouts_Rselec/out.noreplicates/plink.raw"))
liSNPs.m10<- read.PLINK(file = paste0(WD,outfolder,"PopSamples_m10/Popsouts_Rselec/out.noreplicates/plink.raw"))
liSNPs.def<- read.PLINK(file = paste0(WD,outfolder,"PopSamples_def/Popsouts_Rselec/out.noreplicates/plink.raw"))


### Estimate SNP frequency distribution and plot
library(ggplot2)

# Estimate SNPs frequencies for each dataset 
SNPF.m3 <- glMean(liSNPs.m3)
SNPF.m4 <- glMean(liSNPs.m4)
SNPF.m10 <- glMean(liSNPs.m10)
SNPF.def <- glMean(liSNPs.def)


# Create a hist for each dataset 
hm3<-hist(SNPF.m3, proba= TRUE, breaks=30)
hm4<-hist(SNPF.m4, proba= TRUE, breaks=30)
hm10<-hist(SNPF.m10, proba= TRUE, breaks=30)
hdef<-hist(SNPF.def, proba= TRUE, breaks=30)


# plot them all togheter 
plot(hm3, col=rgb(0,0,1,1/2),xlab="SNP Frequency", ylab="Density")
plot(hdef, col=rgb(1,1,0,1/4), add=TRUE)
plot(hm4, col=rgb(0,1,0,1/4), add=TRUE)
plot(hm10, col=rgb(1,0,0,1/2), add=TRUE)

## Rescale all hists to make them comparable
# copy hists
rhm3 <- hm3 
rhm4 <- hm4
rhm10 <- hm10
rhdef <- hdef

# rescale them
rhm3$counts<-hm3$counts/max(hm3$counts) 
rhm4$counts<-hm4$counts/max(hm4$counts)
rhm10$counts<-hm10$counts/max(hm10$counts)
rhdef$counts<-hdef$counts/max(hdef$counts)


# plot them all in pannels
par(mfrow=c(2,2))
plot(rhm3, xlab="SNP Frequency", ylab="Density",  main="a) m3")
plot(rhdef, xlab="SNP Frequency", ylab="Density", main="b) m4")
plot(rhm4, xlab="SNP Frequency", ylab="Density",  main="c) m10")
plot(rhm10, xlab="SNP Frequency", ylab="Density", main="d) def")
par(mfrow=c(1,1))


######## Compare distance between individulas of the same pop
# Estimate distance between individuals of the same pop 

dists<-data.frame(stringsAsFactors=FALSE)
source(paste0(WD,"/bin/dist.pop.R"), echo=TRUE)


# Plot to evaluate differences in distance between indvs by pop
require(ggplot2)
plt<-ggplot(data= dists) + theme_bw()
plt + geom_point(aes(x=Pop, y=dist, colour=factor(param), alpha=0.3))

plt<-ggplot(data= dists) + theme_bw()
plt + geom_boxplot(aes(x=Pop, y=dist, colour=factor(param))) +
  ylab("Genetic distance") + xlab("Population") +
  scale_colour_discrete(name="Parameter", labels=c("def", "m=10", "m=4", "m=3"))
 

######## Pop. differentiation and structure ----


### Distance based methods 
for (i in c("m3", "m4", "m10", "def")){
  
  param <- i
  
  ### For log, Tell wich element of the loop is being analysed
  id <- paste("Results for", param)
  print(id)
  

  ####---1) Create genlight object and summarize basic data
  # Inport using the plink file exported from populations Stacks program
  # using whilelist file (i.e. loci in final matrix) and pop map excluding replicates
  require(adegenet)
  liSNPs<- read.PLINK(file = paste0(WD,outfolder,"PopSamples_", param,"/Popsouts_Rselec/out.noreplicates/plink.raw"))
  
  # Change pop labels to Pop names
  levels(pop(liSNPs))<-c("Aj","An","Iz","Ma","Out", "Pe","Tl","To","Za")
  
  
  ### Check number of recovered samples and loci
  print("number of recovered samples and total number of SNPs")
  print(nInd(liSNPs)) 
  print(nLoc(liSNPs))
  
  #check samples names
  print("Sample names")
  print(liSNPs$ind.names)
  
  # Check missing loci per sample
  liSNPs$gen
  x<- NA.posi(liSNPs)
  
  
  ### Visualize the Matrix as a plot
  par(mfrow = c(1, 1))
  glPlot(liSNPs)
  title(paste("Distribution of genotypes by individual using", param))
  
 
  ####### 2) Compute pairwise distances between individuals 
  
  ## separate the whole genlight object in smaller ones by creating blocks of loci. This facilitates computing
  blocks<- seploc(liSNPs, n.block=5) 
  class(blocks) #check if it is a list
  ## estimate distance matrix between individuals of each block
  D<- lapply(blocks, function(e) dist(as.matrix(e))) 
  names(D) #check names correspond to blocks
  # generate final general distance matrix by summing the distance matrixes
  Df<- Reduce("+", D) 
  

  ###  Compute distance matrix (euclidean)
  dat.d = dist(Df)
  
  
  ## add matinfo to dat.d matrix
  x <- match(labels(dat.d), matinfo$sample) #match labels with samples names
  x <- matinfo[x,] #create a new matrix with the samples from dat.d
  pop <- x$Pop #extract Pop info
  lane <-x$Lane #extract Lane info
  
  ####### 3) Plot NJ dendogram
  require(ape)
  ## Function to color the populations according to a specific color
  # to use to color a tree regardless of the subset of samples used
  pops.col.tree<- function(tree.samples){ #samples present in the tree, ie the tips
    tree.samples = tree$tip.label #takes the samples names from the tips of a given tree (ape)
    idx = match(tree.samples, matinfo$sample) # matches the tips samples names with the samples names of a matrix info
    info.tree = matinfo[idx,] #to add population info to the tree
    pops = info.tree$Pop # the populations are now in the colum Pop, these population names may be a subset of the matrix info, if the tree was built using a subset of samples
    pops = as.factor(pops) 
    levels(pops) = cols.pop.key[ match(levels(pops), cols.pop.key[, 1]), 2] # reassing the names of the populations to its corresponding color key defined before
    pops = as.character(pops) #this is what is read as a color list for col=
  }
  
  ## construct, root and plot tree coloring by population
  
  # Build tree
  tree<-nj(Df)
  
  # Re-scale tree so it will be comparable to others with the same height scale = 100
  # Function from http://blog.phytools.org/2012/02/quicker-way-to-rescale-total-length-of.html
  library(phytools)
  rescaleTree<-function(tree,scale){
    tree$edge.length<-
      tree$edge.length/max(nodeHeights(tree)[,2])*scale
    return(tree)
  }
  tree<-rescaleTree(tree,100)
  
  # root
  tree<-root(tree, outgroup="OutBtAl214")
  
  # plot
  plot(tree, 
    type="phylogram",  #type of tree to draw, try fan, looks cool
    show.tip=TRUE, #take out full sample name
    cex=0.5, underscore=TRUE, # use underscore=TRUE to keep the underscore as it and not a space)                         
    tip.col= pops.col.tree(tip.label) #pops.col is the function defined before to match samples to its population and a specified color
  ) 
  axisPhylo()
  title(paste0("rooted NJ dendogram scaled 0-100 using ", param))
  
  
  ###### 4) Perform and Plot PCoA
  par(mfrow = c(1, 2))
  PCoA_pop(dat.d, cols.pop.key=cols.pop.key, pop, pile = TRUE) #labels = pops
 
  
  #### Take out the Za and Out samples and repeat PCoA
  ## Take out Za from the matrix 
  #look for sample names (and other data columns) that do NOT include "Za" 
  s<-grep("Za|Out", labels(dat.d), value= TRUE, invert= TRUE)
  
  #Transform dist matrix to matrix
  dat.d<-as.matrix(dat.d)
  #select the samples to keep
  dat.d<-dat.d[s,s]
  #transform to dist martrix again
  dat.d<-as.dist(dat.d)
  
  ## Repeat PCoA
  # add matinfo to dat.d matrix
  x <- match(labels(dat.d), matinfo$sample) #m atch labels with samples names
  x <- matinfo[x,] # create a new matrix with the samples from dat.d
  pop <- x$Pop # extract Pop info
  lane <-x$Lane # extract Lane info
  
  ## Perform and Plot PCoA
  PCoA_pop(dat.d, cols.pop.key=cols.pop.key, pop, pile = TRUE) #labels = pops
  
}
  

### Fst
for (i in c("m3", "m4", "m10", "def")){
  
  ### For log, Tell wich element of the loop is being analysed
  id <- paste("Results for", i)
  print(id)
  
##### Estimate mean of Fst distance matrix
## Import Fst pairwise distance matrix created by Stacks population algoritm in the full dataset excluding the replicates
Fstmat<-data.matrix(read.delim(file = paste0(WD,outfolder,"PopSamples_", i,"/Popsouts_Rselec/out.noreplicates/batch_1.fst_summary.tsv"),
  row.names=1, fill=TRUE))  
colnames(Fstmat)<- c("Aj","An","Iz","Ma","Out","Pe","Tl","To","Za") # add col names

print("Pairwise Fst mat")
print(as.table(Fstmat))


}  
  


######## Distribution of private alleles and Ho -----

## Read Stacks population results # All positions (variant and fixed) in file batch_1.sumstats_summary.tsv
# Create a data frame to store the data
pop.gen.variant<-data.frame() # for only variant positions
pop.gen.varfix<-data.frame()  # for variant and fixed positions

# read variable sites summary for each file
for (i in c("m3", "m4", "m10", "def")){
  
  ## For only variant postions
  # read part of interest batch_1.sumstats_summary.tsv
  data<-read.delim(file = paste0(WD,outfolder,"PopSamples_", i,"/Popsouts_Rselec/out.noreplicates/batch_1.sumstats_summary.tsv"),
                   skip= 1, # to skip the first line of the file which only indicates that the following are the variant positions
                   nrows= 9 # to do not read the info corresponding to both fixed and variant positions)  
                   )
  
  # Create vector with pop names in same order Pop ID batch_1.sumstats_summary.tsv
  Pop_name<-c("Aj","An","Iz","Ma","Out","Pe","Tl","To","Za") 
  
  # Create vector with name of variable (Stacks parameter) used to create data
  param<-rep(i,9)
  
  # Add pop names to dataf
  data<-cbind(param,Pop_name,data)
  
  # Reorder
  data<-data[c(1,2,3,4,6,7,8,9,5),]
  
  # add to data frame
  pop.gen.variant<-rbind(pop.gen.variant,data)
  
  
  ## For all (variant and fixed) positions
  # read part of interest batch_1.sumstats_summary.tsv
  data<-read.delim(file = paste0(WD,outfolder,"PopSamples_", i,"/Popsouts_Rselec/out.noreplicates/batch_1.sumstats_summary.tsv"),
    skip= 12) # to skip the only Variant positions info
    
    # Create vector with pop names in same order Pop ID batch_1.sumstats_summary.tsv
    Pop_name<-c("Aj","An","Iz","Ma","Out","Pe","Tl","To","Za") 
    
    # Create vector with name of variable (Stacks parameter) used to create data
    param<-rep(i,9)
    
    # Add pop names to dataf
    data<-cbind(param,Pop_name,data)
    
    # Reorder
    data<-data[c(1,2,3,4,6,7,8,9,5),]
    
    # add to data frame
    pop.gen.varfix<-rbind(pop.gen.varfix,data)
    
}

# Se resulted data
pop.gen.variant
pop.gen.varfix

# Make nice plots
library(ggplot2)
library(grid)

## Private alleles per pop
# Plot absolute Number of private alleles per population 
plt<-ggplot(data= pop.gen.varfix, mapping=aes(x=Pop_name, y=Private)) + theme_bw()
plt<- plt + geom_point(aes(shape=factor(param))) + xlab("Population") + ylab("Absolute number of private SNPs")
plt1<- plt + scale_shape(name="Parameter") + theme(legend.position="none") 


# Plot Number of private alleles per population in proportion to the variant sites found in that pop
plt<-ggplot(data= pop.gen.varfix, mapping=aes(x=Pop_name, y=Private/Variant.Sites)) + theme_bw()
plt<- plt + geom_point(aes(shape=factor(param))) + xlab("Population") + ylab("Relative number of private SNPs")
plt2<-plt + scale_shape(name="Parameter") 

# plot them in two pannels
grid.draw(cbind(ggplotGrob(plt1), ggplotGrob(plt2), size="last"))


## Plot Observed Ho per population considering only the variant sites
plt<-ggplot(data= pop.gen.variant, mapping=aes(x=param, y=Obs.Hom))
plt + geom_bar(stat="identity", fill="white", colour="black", width=.7) + 
  geom_errorbar(aes(ymin=Obs.Hom-StdErr.2, ymax=Obs.Hom+StdErr.2, width=.2)) +
  facet_wrap(~Pop_name) + theme_bw() + xlab("Parameter") + ylab("Observed homozygosity")

## Plot Observed Ho per population considering Variant and fixed sites
plt<-ggplot(data= pop.gen.varfix, mapping=aes(x=param, y=Obs.Hom))
plt + geom_bar(stat="identity", fill="white", colour="black", width=.7) + 
  geom_errorbar(aes(ymin=Obs.Hom-StdErr.2, ymax=Obs.Hom+StdErr.2, width=.2)) +
  facet_wrap(~Pop_name) + theme_bw() + xlab("Parameter") + ylab("Observed homozygosity")

## Plot % Polymorphic loci
plt<-ggplot(data= pop.gen.varfix, mapping=aes(x=param, y=X..Polymorphic.Loci))
plt + geom_bar(stat="identity", fill="white", colour="black", width=.7) + 
    facet_wrap(~Pop_name) + theme_bw() + xlab("Parameter") + ylab("% Polymorphic loci")

## Plot Fis
plt<-ggplot(data= pop.gen.varfix, mapping=aes(x=param, y=Fis))
plt + geom_point(stat="identity", fill="white", colour="black", width=.7) + 
  geom_errorbar(aes(ymin=Fis-StdErr.7, ymax=Fis+StdErr.7, width=.2)) +
  facet_wrap(~Pop_name) + theme_bw() + xlab("Parameter") + ylab(expression(italic(F)[IS]))




######## Summary of Settings for running this script -----
# Working directory
getwd()

# And R specifications
sessionInfo()


