PopMapEQsz<-function(tsv, writedirectory, matinfo, homogp, dsz, PopOrder){
  # Function to create a list of the samples present in a desired matrix 
  # choosing at random a desired sampling size (dsz) number of individuals per pop for 
  # a desired subgroup (homogp) of samples to be homogenized (therefore leaving rest of populations with original sampling size),
  # and write a file that can be read by Stacks as a Population Map of desired samples to analyse with the populations program
  # Note: if replicates (samples ending with _r or _ir) are present they would be excluded
  #
  ## Variables:
  # tsv: the path to the file of the matrix of selected RADloci from the PostCleaning script
  # writedirectory: The directory where the PopMap file should be saved, whitout starting or ending with /
  # dsz: desired sampling size for each population
  # homogp: numeric vector with Pop ID of populations which sampling sizes should be homogenized to dsz 
  # matinfo: a cvs file with sample names in a column "sample", and population information in a column "Pop" 
  # PopOrder: numeric vector with desired order to sort populations. Use to fit the PopKey desired order for your populations change line 43 (levels(x$Pop)….)
  
  # First load the matrix 
  final = read.delim(paste(tsv, sep = ""), header = T) 
  # Extact name of recovered samples
  samples <- colnames(final[, 4:ncol(final)])
  # Extract these samples from the matinfo
  x <- match(samples, matinfo$sample) #match labels with samples names
  x <- matinfo[x,] #create a new matrix with the samples from dat.d
  
  ### Drop replicates
  s <-grep("_r|_ir", x$sample, value= TRUE, invert= TRUE) # drop the samples ending with _r or _ir
  s <- match(s, x$sample) # select the rest
  x <- x[s,]
  # Keep only samples and population columns
  x <- x[, c(1,4)]
  # Change levels of Pop to integers
  levels(x$Pop) <- PopOrder
  # transform to vectors
  x$sample <- as.vector(x$sample) 
  
  ### Focus on desired group to homogenize sz to
  exc<-x[!x$Pop %in% homogp,] # excluded samples 
  homogp <- x[x$Pop %in% homogp,] # desired group
  homogp$Pop<-as.vector(homogp$Pop)
  
  ### Keep only dsz samples chosen at random for each desired population 
  set.seed(1) 
  df<- lapply(split(homogp, homogp$Pop), #split dataframe by Pop
    function(z) z[sample(1:nrow(z), dsz), ]) #get dsz samples of each group
  df<-do.call(rbind, df) # rbind into a single dataframe again
  
  ### Add data of excluded populations
  x<-rbind(df,exc)
  
   
  
  # Write it to a file with the name of the tvs matrix plus _EQsz (indicating equal samplign size)
  write.table(x, file= paste0(writedirectory, "/", basename(tsv), "_EQsz.tsv"), sep = "\t",
    row.names =FALSE, col.names=FALSE, eol = "\n",
    quote=FALSE)  
}