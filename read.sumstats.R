read.sumstats <- function(file, npop, popNames){
  ### Function to read the Stacks populations output file batch_1.sumstats.tsv and add PopNames to it
  # file = path to batch_1.sumstats_summary.tsv file
  # npop = number of populations (PopID field) the datafile contains
  # popNames = vector with population names in the same order than PopID in the file 
  
  ### Get data
  # read part of interest batch_1.sumstats.tsv
    data<-read.delim(file = file, skip= npop) # to skip the first lines of the file that only indicate the pops key
    Pop.Name<-as.factor(data$Pop.ID)
    levels(Pop.Name)<-popNames
  
  # Add pop names to data
    data<-cbind(Pop.Name,data)
  }
  