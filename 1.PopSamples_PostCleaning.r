rm(list = ls())
WD<-"/gpfs/bio/bxh10nyu/BerL_1_2_3/3Berberis_phylogeo/" #directory to work on the cluster
# WD<-"/Volumes/TO_GO/BerL_1_2_3/2R/3Berberis_phylogeo/"
setwd(WD) 
list.files()

###Load and save data:
workfolder = "savedanalyses" #leave as it is
outfolder = "data.out" #change to folder where out data will be
infolder = "data.in/" # change to where the .snp and .cov matrices are


#load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"info/Ber_06oct13.info", sep = ""), header = T)


# DEBUG params (should be commented in production version)
ncores = 3
# use a low number for bestscorers to recover all the samples 
bestscorers = .3 #parameter to pick-up best specimens. these are defined as having more than bestscorers * mean(NrStacks / specimen)
minTaxDepth = .8
maxAlleles = 5
maxlenSNP = 5
keepmono = 0 #0 = remove monomorphic stacks or 1 leave them in there
#minCov = 3 #added this parameter, desired min coverage of all alles in a locus to make it a valid target locus 


# for log purposes, display input parameters. Leave as it is.
ncores
bestscorers
minTaxDepth
maxAlleles
#minCov


######## GENERAL TOOLS

#### Colors functions that will be used for plots

## to color a plot according to LANE group
#define list of nice colors for lanes
colors <-rainbow(10)
lane.names<-as.character(c(1:10))

# bind pops to colors to make a stable color key
cols.lane.key =cbind(lane.names, colors)
cols.lane.key

# function
lane.col<- function(x){ # x is the dataframe to plot
  plot.samples = names(x) #takes the samples names from the samples in the data we want to plot
  idx = match(plot.samples, matinfo$SampleSEQ.ID) # matches the samples names with the samples names of a matrix info
  info.plot = matinfo[idx,] #to add the information of the matinfo to the plot
  lanes = info.plot$Lane # the lane are now in the colum Lane, 
  lanes = as.factor(lanes)
  levels(lanes) = cols.lane.key[match(levels(lanes), cols.lane.key[, 1]), 2] # reassing the names of the lanes to its corresponding color key defined before
  lanes = as.character(lanes) #this is what is read as a color list for col=
}


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
mcov_cell = function(val){
  nall = NA
  if(is.na(val) == F){
    if(val != ""){
      nall = mean(as.numeric(strsplit(val, "/")[[1]]))
    }
  } else {
    nall = NA
  }
  nall 
}
mcov_row = function(vec) sapply(vec, mcov_cell)

#Same funtion to estimate coverage modified by Alicia (added as.character(val))
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

#Function to estimate min. coverage of the two alles (if present) or of the loci per each loci and sample

min.cov_cell = function(val){
  nall = NA
  if(is.na(val) == F){
    if(val != ""){
      nall = min(as.numeric(strsplit(as.character(val), "/")[[1]]))
    }
  } else {
    nall = NA
  }
  nall 
}
min.cov_row = function(vec) sapply(vec, min.cov_cell)


### Get number of alleles per locus
nall_loc = function(locus){
  nall = nlevels(as.factor(unlist(strsplit(locus, "/"))))
}



######## ANALYSIS
## Create a function to the whole analysis for every different Param data set.  

PostCleaning <- function(tsv, cov, directory,nsamps){ 
  # tsv and cov are the SNP and COV matrices exported with export_sql.pl from stacks
  # directory is the name of the directory inside data.in and data.out where matrices are
    tsv = tsv
    cov = cov
  
    # Ok, now first, we check whether these data have not been already processed
    path_to_workfiles = dir(workfolder, pattern = ".Rdata")
    path_to_workfiles = path_to_workfiles[grep(basename(tsv), path_to_workfiles)] #make sure we pick the right file in there

    if(length(path_to_workfiles) == 0){ #no file saved from a previous session is available, we have to compute everything.
     ## Get data
     # initiate multithreading
    library(doMC)
    library(foreach)
    library(multicore)
    library(R.utils)
    registerDoMC(cores=ncores)
  
    # open files
    # Stacks produces a file with the sample names at the end, to do not read them and load data properly:
    # first count the number of lines
    # and rest according to the number of samples that the matrix has 
    maxl<-countLines(paste(WD,infolder, directory,"/", tsv, sep = "")) - nsamps
    data = read.delim((paste(WD,infolder, directory,"/", tsv, sep = "")), header = T, row.names = 1, nrows=maxl) 
    data.cov = read.delim((paste(WD,infolder, directory,"/", cov, sep = "")), header = T, row.names = 1, nrows=maxl)
      

    cls = colnames(data)
    colnames(data) = cls
  
    # get generic filenames
    tsv = basename(tsv)
    cov = basename(cov)
  
    # get out occurrences
    occs = data[, 12:ncol(data)]
    occs.cov = data.cov[, 12:ncol(data.cov)]
    occs = as.matrix(occs)
    occs.cov = as.matrix(occs.cov)
  
    # Compute number of alleles per sample
    # NALL = t(apply(occs, 1, nall_row)) #usual, single core method
    NALL = foreach(i=1:nrow(occs), .combine=rbind)%dopar%{ #multithreaded
    nall_row(occs[i, ])
      }
    rownames(NALL) = rownames(occs)
  
    # compute average coverage per rad / sample (average over alleles)
    MCOV = foreach(i=1:nrow(occs.cov), .combine=rbind)%dopar%{ #multithreaded
      mcov_row(occs.cov[i, ])
      }
    rownames(MCOV) = rownames(occs.cov)
  
    # compute MIN coverage per rad / sample (min over alleles)
    #  MinCOV = foreach(i=1:nrow(occs.cov), .combine=rbind)%dopar%{ #multithreaded
    #    min.cov_row(occs.cov[i, ])
    #  }
    #  rownames(MinCOV) = rownames(occs.cov)
  
    # compute number of alleles segregating in the sampling at each locus
    NSEG = foreach(i=1:nrow(occs), .combine=rbind)%dopar%{ #multithreaded
      nall_loc(occs[i, ])
    }
    rownames(NSEG) = rownames(occs)
  
    # make sure we save the produced file in a binary format. will save time for future tests
    save(list = c("data", "data.cov", "tsv", "cov", "NALL", "MCOV", "NSEG"), file = paste(workfolder, "/WorkingFiles_", tsv, ".Rdata", sep = ""))
  
    } else { #Bingo, we have a file already in here
      ## Get data
      load(paste(workfolder, path_to_workfiles, sep = "/"))
  
    # get out occurrences
    occs = data[, 12:ncol(data)]
    occs.cov = data.cov[, 12:ncol(data.cov)]
    occs = as.matrix(occs)
    occs.cov = as.matrix(occs.cov)
  
    }

    ##### Collect informations about the dataset
    # Count how many stacks are present in each specimen
    OUT = NALL
    OUT[OUT > 1] = 1
    OUT[is.na(OUT)] = 0
    LocDepth = colSums(OUT)

    # Keep only best scoring specimens (selection made relative to mean)
    passedspecimens = which(LocDepth > bestscorers * mean(LocDepth))
    OUT = OUT[, passedspecimens]
    NALL = NALL[, passedspecimens]

    # Detect loci with more than 2 alleles / specimen (filtering step applied later)
    TaxDepth = rowMeans(OUT)

    # Detect loci with more than 2 alleles / specimen (filtering step applied later)
    nbAlleles = apply(NALL, 1, max, na.rm = T)

    # Get how many SNPs segregate per locus
    lenSNP = data$Num.SNPs

    # Check whether a given locus is polymorphic
    mono = occs == "consensus"
    mono = rowSums(mono, na.rm = T)

    # Detect the min coverage per per RADloci / sample (min over alleles) in the passed specimens
    #LocCov<-MinCOV[, passedspecimens] #Select the min coverage of the passed specimens
    #LocCov<-apply(LocCov, 1, min, na.rm=T) #detect loci min coverage )filtering step applied later)

    ##### Clean dataset
    # identify STACKs of interest
    # Noticie I added a LocCov min value,
    target = which(TaxDepth > minTaxDepth & nbAlleles < maxAlleles & lenSNP < maxlenSNP & mono == keepmono)

    #you can add LocCov >= minCov

    #Visualize Number of stacks per sample
    summary(LocDepth)
    barplot(LocDepth,
    ylab="Number of stacks",
    xlab="sample",
    cex.names=0.3,
    las=3,#for vertical labels
    col= lane.col(LocDepth)) 
    title("LocDepth per sample")

    #See which specimens failed
    failspecimens <- setdiff(names(LocDepth), names(passedspecimens)) #compare passed specimens against allspeciemns
    idx<-match(failspecimens, matinfo$sample)  #match agaist matinfo
    failspecimens<-matinfo[idx,]  #add matinfo to failspecimens
    failspecimens<-as.data.frame(failspecimens)
    failspecimens


    ## filter by STACKS of interest in the dataset
    final = data.frame(CatalogID = rownames(data[target, ]), Consensus=data[target, ]$Consensus.Sequence, SNPs = data[target, ]$SNPs, occs[ target, passedspecimens])
    final.cov = data.frame(CatalogID = rownames(data[target, ]), Consensus=data[target, ]$Consensus.Sequence, SNPs = data[target, ]$SNPs, occs.cov[ target, passedspecimens])
    final[final == ""] = -9
    final.cov[final.cov == ""] = -9

    final.t = data.frame(CatalogID = rownames(data[target, ]), Consensus=data[target, ]$Consensus.Sequence, SNPs = data[target, ]$SNPs, occs[ target, passedspecimens])
    final.cov.t = data.frame(CatalogID = rownames(data[target, ]), Consensus=data[target, ]$Consensus.Sequence, SNPs = data[target, ]$SNPs, occs.cov[ target, passedspecimens])
    final.t[final.t == ""] = -9
    final.cov.t[final.cov.t == ""] = -9

    # save out the result
    write.table(final, file = paste(outfolder,"/", directory, "/", tsv, ".SNPs", sep = ""), quote = F, sep = "\t", row.names = F)
    write.table(final.cov, file = paste(outfolder,"/", directory, "/", cov, ".COVs", sep = ""), quote = F, sep = "\t", row.names = F)

  ### See distribution of SNPs by positon
  final = read.delim(paste(paste(WD, outfolder,"/", directory, "/", tsv, sep = ""), ".SNPs", sep = ""), header = T) 
  library(stringr)
  SNPs_pos<-str_extract_all(final$SNPs, #get the column of the snps positon  
    "[0-9.]+") #this extracts only the numbers
  occur.SNPs<-unlist(SNPs_pos) #unlist to have a list of all total ocurrences
  occur.SNPs<-as.numeric(occur.SNPs) #as numeric to be able to order 
  # create table with number of times that a SNP is observed in each positon, transform to dataframe
  occur.SNPs<-table(occur.SNPs) 
  occur.SNPs<-as.data.frame(occur.SNPs)
  row.names(occur.SNPs)<-occur.SNPs$occur.SNPs #make rows names be SNPs pos 
  ## and plot
  # take from tsv name of parameter
  param<-gsub(tsv, pattern="($SNP){0,1}(BerExplo_){0,1}", replacement="")
  # plot
  par(las=2)  # make label text perpendicular to axis
  barplot(occur.SNPs$Freq,
    cex.names= 0.5, las=3,
    names.arg=occur.SNPs$occur.SNPs,
    xlab= "position in the RAD read",
    ylab= "Number of SNPs",
    main= paste("Number of SNPs per read position for", param)
    )

}


## Run the function with every different dataset
pdf(file=paste(WD,"docs/3Berberis_phylogeo_PostCleaning.pdf", sep = ""))
PostCleaning("3Berberis_phylogeo_BeralpBt_m3.SNP", "3Berberis_phylogeo_BeralpBt_m3.COV", "3Berberis_phylogeo_m3", 92)
dev.off()
