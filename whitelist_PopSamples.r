rm(list = ls())

# Define WD and other directories
WD<-"/Volumes/TO_GO/BerL_1_2_3/2R/PopSamples" 
setwd(WD) 
list.files()

outfolder = "/data.out/" # tvsdirectory, where the SNP matrices are

### load samples meta information (lane, barcode, pop, etc)
matinfo = read.delim(paste(WD,"/info/Ber_06oct13.info", sep = ""), header = T)

## Fuction to create a list of the CatalogID loci present in a desired matrix 
# and write a file that can be read by Stacks as whitelist of desired loci
# 
whiteRADlist<-function(tsv, writedirectory){
  # Fuction to create a list of the CatalogID loci present in a desired matrix 
  # and write a file that can be read by Stacks as whitelist of desired loci
  # 
  ## Variables:
  # tsv the path to the file of the matrix of selected RADloci from the PostCleaning script
  # writedirectory The directory where the whitelist file should be saved, , whitout starting or ending /
  
  # First load the matrix 
  final = read.delim(paste(tsv, sep = ""), header = T) 
  # Extact name of loci (CatalogID)
  whitelist <- final$CatalogID
  # Write it to a file with the name of the tvs matrix plus _whitelist
  write(whitelist, file= paste0(writedirectory, "/", basename(tsv), "_whitelist.tsv"), ncolumns = 1)
}

# Run for m3, m4, m10 and def

whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_m4/PopSamples_BeralpBt_m4.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_m10/PopSamples_BeralpBt_m10.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whiteRADlist(tsv= paste0(WD,outfolder,"PopSamples_def/PopSamples_BeralpBt_def.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))


### Function to create a list of the samples present in a desired matrix
# and write a file that can be read by Stacks as a Population Map of desired samples to analyse with the populations program
whitePopMap<-function(tsv, writedirectory, drop.rep){
  # Function to create a list of the samples present in a desired matrix
  # and write a file that can be read by Stacks as a Population Map of desired samples to analyse with the populations program
  #
  # The PopKey will be stable Aj:1, An:2, Iz:3, Ma:4, Pe:5, Tl:6, To:7, Za:8, Out:9
  #
  ## Variables:
  # tsv: the path to the file of the matrix of selected RADloci from the PostCleaning script
  # writedirectory: The directory where the PopMap file should be saved, whitout starting or ending /
  # drop.rep: if TRUE replicates (sampled ending with _r or _ir) are discarded, if FALSE they are kept
  
  # First load the matrix 
  final = read.delim(paste(tsv, sep = ""), header = T) 
  # Extact name of recovered samples
  samples <- colnames(final[, 4:ncol(final)])
  # Extract these samples from the matinfo
  x <- match(samples, matinfo$sample) #match labels with samples names
  x <- matinfo[x,] #create a new matrix with the samples from dat.d
  
  if(drop.rep == TRUE){
  # Drop replicates
  s <-grep("_r|_ir", x$sample, value= TRUE, invert= TRUE) # drop the samples ending with _r or _ir
  s <- match(s, x$sample) # select the rest
  x <- x[s,]
  # Keep only samples and population columns
  x <- x[, c(1,4)]
  # Change levels of Pop to integers
  levels(x$Pop) <- c(1,2,3,4,9,5,6,7,8)
  # transform Pop to vectors
  x$Pop <- as.integer(x$Pop)
  x$sample <- as.vector(x$sample)
  str(x)
  # Write it to a file with the name of the tvs matrix plus _PopMap_withrep (indicating reps are included)
  write.table(x, file= paste0(writedirectory, "/", basename(tsv), "_PopMap_norep.tsv"), sep = "\t",
    row.names =FALSE, col.names=FALSE,
    quote=FALSE)
  } else {
    # Don't drop replicates
    # Keep only samples and population columns
    x <- x[, c(1,4)]
    # Change levels of Pop to integers
    levels(x$Pop) <- c(1,2,3,4,9,5,6,7,8)
    # transform Pop to vectors
    x$Pop <- as.integer(x$Pop)
    x$sample <- as.vector(x$sample)
    str(x)
    # Write it to a file with the name of the tvs matrix plus _PopMap_norep (indicating reps are excluded)
    write.table(x, file= paste0(writedirectory, "/", basename(tsv), "_PopMap_withrep.tsv"), sep = "\t",
      row.names =FALSE, col.names=FALSE,
      quote=FALSE)   
  }
  
}

# Run for m3, m4, m10 and def EXCLUDING replicates
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_m4/PopSamples_BeralpBt_m4.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_m10/PopSamples_BeralpBt_m10.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whitePopMap(drop.rep = TRUE, tsv= paste0(WD,outfolder,"PopSamples_def/PopSamples_BeralpBt_def.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))

# Run for m3, m4, m10 and def INCLUDING replicates
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_m3/PopSamples_BeralpBt_m3.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_m4/PopSamples_BeralpBt_m4.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_m10/PopSamples_BeralpBt_m10.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))
whitePopMap(drop.rep = FALSE, tsv= paste0(WD,outfolder,"PopSamples_def/PopSamples_BeralpBt_def.SNP.SNPs"), writedirectory= paste0(WD,"/docs"))


