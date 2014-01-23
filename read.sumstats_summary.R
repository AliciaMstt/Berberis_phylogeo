read.sumstats_summary <- function(file, allP, npop, popNames){
  ### Function to read the Stacks populations output file batch_1.sumstats_summary.tsv
  # file = path to batch_1.sumstats_summary.tsv file
  # allP = wheater to read the information from all variant and fixed positions (TRUE), or only the variant (FALSE)
  # npop = number of populations (PopID field) the datafile contains
  # popNames = vector with population names in the same order than PopID in the file 
  
  if (allP==TRUE){ 
    ## Read data for all (variant and fixed) positions
    # read part of interest batch_1.sumstats_summary.tsv
    data<-read.delim(file = file, skip= npop+3) # to skip the only Variant positions info

    # Add pop names to data
    data<-cbind(popNames,data)
    
  } else {
    ## Read data for only variant postions
    # read part of interest batch_1.sumstats_summary.tsv
    data<-read.delim(file = file, skip= 1, # to skip the first line of the file which only indicates that the following are the variant positions
                     nrows= npop) # to do not read the info corresponding to both fixed and variant positions)  
    
    # Add pop names to data
    data<-cbind(popNames,data)
  }
 
}

