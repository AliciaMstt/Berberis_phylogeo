rm(list = ls())

# Define WD and other directories
WD<-"/Volumes/TO_GO_1/BerL_1_2_3/3Berberis_phylogeo" 
setwd(WD) 
list.files()

Popoutsfolder = "/data.out/PopSamples_m3/Popsouts_RselecStructure" # where the outputs from running Structure populations are

### load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

### Source functions to read Stacks populations output
source(paste0(WD,"/bin/read.sumstats_summary.R"))
source(paste0(WD,"/bin/read.sumstats.R"))
source(paste0(WD,"/bin/read.fst_summary.R"))

############# ----- ANALYSES -----

library(ggplot2)

# load Stacks output file with populations summary statistics 
popsumstats <-read.sumstats(file=paste0(WD,Popoutsfolder,"/out.noreplicates/batch_1.sumstats.tsv"),
  npop=9, popNames=c("Aj","An","Iz","Ma","Pe","Tl","To","Za","Out"))                  


##### Allele frequency spectrum distribution for loci in each population that are SNPs in Berberis
# Plot the frequency of the major allele (P)
plt <- ggplot(data=popsumstats, aes(x=P)) + theme_bw()
plt + geom_histogram(aes(y =..density..), colour="black", fill="white") +
  ylab("Percentage of loci") + xlab("SNP Frequency")

# Plot the frequency of the major allele (P) by population
plt <- ggplot(data=popsumstats, aes(x=P)) + theme_bw()
plt <- plt + geom_histogram(aes(y =..density..), colour="black", fill="white")  +
  ylab("Percentage of loci") + xlab("SNP Frequency")
plt + facet_wrap(~Pop.Name)


#### Frequency distribution of FIS values across loci within each population.
# Plot FIS for all pops
plt <- ggplot(data=popsumstats, aes(x=Fis)) + theme_bw()
plt + geom_histogram(aes(y =..count..), colour="black", fill="white") +
  ylab("Number of loci") + xlab("Fis value")

# Plot FIS by population
plt <- ggplot(data=popsumstats, aes(x=Fis)) + theme_bw()
plt<- plt + geom_histogram(aes(y =..count..), colour="black", fill="white") +
  ylab("Number of loci") + xlab("Fis value")
plt + facet_wrap(~Pop.Name)




