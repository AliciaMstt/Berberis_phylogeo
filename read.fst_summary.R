read.fst_summary <- function(file, popNames){
    ### Function to read the Stacks populations output file batch_1.fst_summary.tsv and add PopNames to it
    # file = path to batch_1.fst_summary.tsv file
    # popNames = vector with population names in the same order than PopID in the file 
    
    ### Get data
    Fstmat<-data.matrix(read.delim(file = file, row.names=1, fill=TRUE)) 
    # add col names
    colnames(Fstmat)<- popNames 
    Fstmat    
  }
  