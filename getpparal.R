getpparal<-function(popsumstats, mp){
  library(plyr)
  ## Function to get a list of potential paralog loci from the batch_1.sumstats.tsv output of Stacks 
  ## paralogs are defined as those loci are P=0.5 in a given number mp populations. Requieres package plyr
  # popsumsstats: dataframe with column info as in batch_1.sumstats.tsv 
  # mp: minimun number of populations in which loci P=0.5 should be to be considered paralogs
  
  ### Function:
  ## Extact which loci are P=0.5 in any pop 
  P0.5<-popsumstats$P==0.5
  P0.5<-popsumstats[P0.5, ]
  
  ## extract which of those are in at least mp populations 
  df<-ddply(transform(P0.5,id =interaction(Locus.ID,BP)),.(id),
    function(x)if(length(unique(x$Pop.ID))>=mp)x)
  
  # output loci ID
  levels(as.factor(df$Locus.ID))
}