rm(list = ls())

# Define WD and other directories
WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
setwd(WD) 
list.files()

Popoutsfolder = "/data.out/PopSamples_m3" # where the outputs from running Stacks populations are

### load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

### Source functions to read Stacks populations output
source(paste0(WD,"/bin/read.sumstats_summary.R"))
source(paste0(WD,"/bin/read.sumstats.R"))
source(paste0(WD,"/bin/read.fst_summary.R"))


###### Define colors 

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

############# ----- ANALYSES -----

# In the following analyses the subset of samples correspond to:
# * BerAll: all populations from Berberis alpina (including Za), Berberis moranensis (An population), Berberis trifolia (outgroup).
# * BerwoOut: all populations from Berberis alpina (including Za), Berberis moranensis (An population) but EXCLUDING outgroup (B. trifolia)
# * woZaOut: excluding samples from El Zamorano population (Za) and Berberis trifolia (outgroup)
# * BerSS: Berberis alpina sensu stricto, populations (Aj, Iz, Ma, Pe, Tl, To) ie Berall excluding Za, Out and An.
#
# The subset of loci correspond to:
# ** ExcludingParalogs: excluding the putatively paralogous loci (see scripts 5.*)
# ** IncludingParalogs: all RAD loci


# Function to plots FIS and SNP frequency spectrum 
library(ggplot2)
SNPFis.plots<- function (popsumstats){
  ### Function to plot the Allele frequency spectrum (for all pops and by pop), 
  ### and the FIS values across loci within each population (for all popbs and by pop)
  # popsumstats: output of read.sumstats, ie Stacks populations output file batch_1.sumstats.tsv with PopNames added to it
  
  ##### Allele frequency spectrum distribution for loci in each population that are SNPs in Berberis
  # Plot the frequency of the major allele (P)
  plt <- ggplot(data=popsumstats, aes(x=P)) + theme_bw()
  plt1 <- plt + geom_histogram(aes(y =..density..), colour="black", fill="white") +
    ylab("Percentage of loci") + xlab("Allele Frequency")
  
  
  # Plot the frequency of the major allele (P) by population
  plt <- ggplot(data=popsumstats, aes(x=P)) + theme_bw()
  plt <- plt + geom_histogram(aes(y =..density..), colour="black", fill="white")  +
    ylab("Percentage of loci") + xlab("Allele Frequency")
  plt2 <- plt + facet_wrap(~Pop.Name)
  
  
  #### Frequency distribution of FIS values across loci within each population.
  # Plot FIS for all pops
  plt <- ggplot(data=popsumstats, aes(x=Fis)) + theme_bw()
  plt3 <- plt + geom_histogram(aes(y =..count..), colour="black", fill="white") +
    ylab("Number of loci") + xlab("Fis value")
  
  
  # Plot FIS by population
  plt <- ggplot(data=popsumstats, aes(x=Fis)) + theme_bw()
  plt<- plt + geom_histogram(aes(y =..count..), colour="black", fill="white") +
    ylab("Number of loci") + xlab("Fis value")
  plt4 <- plt + facet_wrap(~Pop.Name)
  
  # plot all plots
  return(list(plt1, plt2, plt3, plt4))
}

# Function to give number of loci and samples and to plot a PCA and a NJ tree
infoPCATree<- function(plinkfile, popNames, outgroup){ 
  ### Function to give number of samples and loci, to plot a NJ tree and PCA
  # plinkfile = path to the plink.raw file
  # popNames = population names that would be used for ploting
  # outgroup = sample name of the outgroup used to root the NJ tree
  
  ####---1) Create genlight object and summarize basic data
  # Inport using the plink file exported from populations Stacks program
  # using whilelist file (i.e. loci in final matrix) and pop map excluding replicates
  require(adegenet)
  liSNPs<- read.PLINK(file=plinkfile)
  
  # Change pop labels to Pop names
  levels(pop(liSNPs))<-popNames
  
  ### Check number of recovered samples and loci
  print("number of recovered samples and total number of SNPs")
  print(nInd(liSNPs)) 
  print(nLoc(liSNPs))
  
  ### Visualize the Matrix as a plot
  par(mfrow = c(1, 1))
  glPlot(liSNPs)
  
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
  
  # root
  tree<-root(tree, outgroup=outgroup)
  
  # plot
  plot(tree, 
    type="phylogram",  #type of tree to draw, try fan, looks cool
    show.tip=TRUE, #take out full sample name
    cex=0.5, underscore=TRUE, # use underscore=TRUE to keep the underscore as it and not a space)                         
    tip.col= pops.col.tree(tip.label) #pops.col is the function defined before to match samples to its population and a specified color
  ) 
  
  ## PCoA display by POPULATION and defined colors
  source(paste0(WD,"/bin/PCoA_pop.r"))
  par(mfrow = c(1, 2))
  PCoA_pop(dat.d, cols.pop.key=cols.pop.key, pop, pile = TRUE) #labels = pops
}


#### Function to estimate the no. of loci and SNPs as well 
# as the mean coverage, the loci, allele and SNP error rates
Summary_perfo <- function(id, final, final.cov, plink){
  # id is a chracter string to identify the dataset 
  # final and final.cov are the SNP.SNPs and COV.COVs matrices procuded by PostCleaning.r
  # plink is the path to the plink.raw file to be used to estimate SNP error rate
  # final, final.cov and plink paths should be given from outfolder onwards
  
  # get files
  final = final
  final.cov= final.cov
  
  #### See number of loci obtained
  n.loci<-nrow(final)
  
  #### General tools
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
  liSNPs<- read.PLINK(file = plink)
  
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
  
    
  #### Return summary of results
  summ.performance<- cbind(id, n.loci, n.SNPs, m.cov, sd.cov,
    m.locus.error, sd.locus.error,
    m.allele.error, sd.allele.error,
    m.SNP.error, sd.SNP.error)
  
  
  return(summ.performance)
  
}


############################ IncludingParalogs
PopoutsfolderINC = paste0(Popoutsfolder, "/IncludingParalogs")

############### Information content and error rates
####### BerAll
# define paths
final=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.SNP.SNPs")
final.cov=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.COV.COVs")
plink=paste0(WD,PopoutsfolderINC,"/AllLoci/BerAll/out.replicates/plink.raw")
# open files
final = read.delim(final, header = T) 
final.cov= read.delim(final.cov, header = T) 
# Run Summary performance function
Summary_perfo(id="BerAll", final, final.cov, plink)

####### BerwoOut
# define paths
final=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.SNP.SNPs")
final.cov=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.COV.COVs")
plink=paste0(WD,PopoutsfolderINC,"/AllLoci/BerwoOut/out.replicates/plink.raw")
# open files
final = read.delim(final, header = T) 
final.cov= read.delim(final.cov, header = T) 
# Subset data
s<-grep("Out", colnames(final), invert= TRUE)
final <- final[,s]
s<-grep("Out", colnames(final.cov), invert= TRUE)
final.cov <- final.cov[,s]
# Run Summary performance function
Summary_perfo(id="BerwoOut", final, final.cov, plink)

####### woZaOut
# define paths
plink=paste0(WD,PopoutsfolderINC,"/AllLoci/woZaOut/out.replicates/plink.raw")
# Subset data (from previously subset)
s<-grep("Za", colnames(final), invert= TRUE)
final <- final[,s]
s<-grep("Za", colnames(final.cov), invert= TRUE)
final.cov <- final.cov[,s]
# Run Summary performance function
Summary_perfo(id="woZaOut", final, final.cov, plink)

####### BerSS
# define paths
plink=paste0(WD,PopoutsfolderINC,"/AllLoci/BerSS/out.replicates/plink.raw")
# Subset data (from previously subset)
s<-grep("An", colnames(final), invert= TRUE)
final <- final[,s]
s<-grep("An", colnames(final.cov), invert= TRUE)
final.cov <- final.cov[,s]
# Run Summary performance function
Summary_perfo(id="BerSS", final, final.cov, plink)


############### Plots of FIS and SNP frequency spectrum 
####### BerAll
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

####### BerwoOut
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

####### woZaOut
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)


####### BerSS
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

########### PCA and NJTree

####### BerAll
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")



####### BerwoOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

####### woZaOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

####### BerSS
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")



############################ ExcludingParalogs
PopoutsfolderEXC = paste0(Popoutsfolder, "/ExcludingParalogs")
paralogs<-read.delim(paste0(WD, "/docs/potentialparalogs"), header=F)[,1]



############### Information content and error rates
####### BerAll
# define paths
final=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.SNP.SNPs")
final.cov=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.COV.COVs")
plink=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerAll/out.replicates/plink.raw")
# open files
final = read.delim(final, header = T) 
final.cov= read.delim(final.cov, header = T) 
# Subset data
l<-!final$CatalogID %in% paralogs # loci not in paralogs
final<-final[l,]
l<-!final.cov$CatalogID %in% paralogs # loci not in paralogs
final.cov<-final.cov[l,]
# Run Summary performance function
Summary_perfo(id="BerAll", final, final.cov, plink)

####### BerwoOut
# define paths
final=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.SNP.SNPs")
final.cov=paste0(WD,Popoutsfolder, "/PopSamples_BeralpBt_m3.COV.COVs")
plink=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerwoOut/out.replicates/plink.raw")
# open files
final = read.delim(final, header = T) 
final.cov= read.delim(final.cov, header = T) 
# Subset data
s<-grep("Out", colnames(final), invert= TRUE)
final <- final[,s]
s<-grep("Out", colnames(final.cov), invert= TRUE)
final.cov <- final.cov[,s]
l<-!final$CatalogID %in% paralogs # loci not in paralogs
final<-final[l,]
l<-!final.cov$CatalogID %in% paralogs # loci not in paralogs
final.cov<-final.cov[l,]
# Run Summary performance function
Summary_perfo(id="BerwoOut", final, final.cov, plink)

####### woZaOut
# define paths
plink=paste0(WD,PopoutsfolderEXC,"/AllLoci/woZaOut/out.replicates/plink.raw")
# Subset data (from previously subset)
s<-grep("Za", colnames(final), invert= TRUE)
final <- final[,s]
s<-grep("Za", colnames(final.cov), invert= TRUE)
final.cov <- final.cov[,s]
l<-!final$CatalogID %in% paralogs # loci not in paralogs
final<-final[l,]
l<-!final.cov$CatalogID %in% paralogs # loci not in paralogs
final.cov<-final.cov[l,]
# Run Summary performance function
Summary_perfo(id="woZaOut", final, final.cov, plink)

####### BerSS
# define paths
plink=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerSS/out.replicates/plink.raw")
# Subset data (from previously subset)
s<-grep("An", colnames(final), invert= TRUE)
final <- final[,s]
s<-grep("An", colnames(final.cov), invert= TRUE)
final.cov <- final.cov[,s]
l<-!final$CatalogID %in% paralogs # loci not in paralogs
final<-final[l,]
l<-!final.cov$CatalogID %in% paralogs # loci not in paralogs
final.cov<-final.cov[l,]
# Run Summary performance function
Summary_perfo(id="BerSS", final, final.cov, plink)


############### Plots of FIS and SNP frequency spectrum 
####### BerAll
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

####### BerwoOut
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

####### woZaOut
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)


####### BerSS
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

########### PCA and NJTree

####### BerAll
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")



####### BerwoOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

####### woZaOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

####### BerSS
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")


