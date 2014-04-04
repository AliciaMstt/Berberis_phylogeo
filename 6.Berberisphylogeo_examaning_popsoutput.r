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
# And inside following subcategories:
# * AllLoci: all RAD loci
# * HitGP: loci that positively blasted against green plants
#   -cpDNA subset of HitGP that blasted to chloroplast
#   -mtDNA subset of HitGP that blasted to mithocondrium 
# * NohitGP: loci that didn't blasted against green plants





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


############################ IncludingParalogs
############### Plots of FIS and SNP frequency spectrum 
PopoutsfolderINC = paste0(Popoutsfolder, "/IncludingParalogs")
####### BerAll
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/NohitGP/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

####### BerwoOut
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/NohitGP/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

####### woZaOut
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/NohitGP/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)


####### BerSS
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/AllLoci/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderINC,"/NohitGP/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)


########### general info PCA and NJTree

####### BerAll
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/NohitGP/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")

####### BerwoOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/NohitGP/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

####### woZaOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/NohitGP/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

####### BerSS
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/AllLoci/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/HitGP/HitGP/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderINC,"/NohitGP/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")



############################ ExcludingParalogs

############# Plots of FIS and SNP frequency spectrum 
PopoutsfolderEXC = paste0(Popoutsfolder, "/ExcludingParalogs")
####### BerAll
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/NohitGP/BerAll/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  
SNPFis.plots(popsumstats=popsumstats)

####### BerwoOut
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/NohitGP/BerwoOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=8, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"))                  
SNPFis.plots(popsumstats=popsumstats)

####### woZaOut
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/NohitGP/woZaOut/out.noreplicates/batch_1.sumstats.tsv"),
  npop=7, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)


####### BerSS
## AllLoci
# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/AllLoci/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## HitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)

## NohitGP
popsumstats <-read.sumstats(file=paste0(WD,PopoutsfolderEXC,"/NohitGP/BerSS/out.noreplicates/batch_1.sumstats.tsv"),
  npop=6, popNames=c("Aj","Iz","Ma","Pe","Tl","To"))                  
SNPFis.plots(popsumstats=popsumstats)


################ general info PCA and NJTree
####### BerAll
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/NohitGP/BerAll/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"),
  outgroup= "OutBtAl214")

####### BerwoOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/NohitGP/BerwoOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za"),
  outgroup= "ZaB06")

####### woZaOut
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/NohitGP/woZaOut/out.noreplicates/plink.raw"),
  popNames=c("Aj","An","Iz","Ma","Pe","Tl","To"),
  outgroup= "AnB01")

####### BerSS
## AllLoci
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/AllLoci/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")

## HitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/HitGP/HitGP/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")

## NohitGP
infoPCATree(plinkfile = paste0(WD,PopoutsfolderEXC,"/NohitGP/BerSS/out.noreplicates/plink.raw"),
  popNames=c("Aj","Iz","Ma","Pe","Tl","To"),
  outgroup= "PeB01")


