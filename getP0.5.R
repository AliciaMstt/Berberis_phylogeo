getP0.5<-function(popsumstats, mp){

  ## Function subset the potential paralog loci from the batch_1.sumstats.tsv output of Stacks 
  ## paralogs are defined as those loci are P=0.5 in at least mp populations. 
  # popsumsstats: dataframe with column info as in batch_1.sumstats.tsv 
  # mp: minimun number of populations in which loci P=0.5 should be to be considered paralogs
  require(plyr)

  # Get all Locus where P=0.5 in at least one pop
  P0.5<-popsumstats$P==0.5
  P0.5<-popsumstats[P0.5, ]
  P0.5<-P0.5$Locus.ID
  
  # Put them in context of the P found for that same locus in other populations
  lc05<-popsumstats$Locus.ID %in% P0.5
  lc05<-popsumstats[lc05,]
  
  # Split data to have managable sizes
  dfs<-split(lc05,lc05$Locus.ID)
  
  # Create column id with locus+BP info 
  dfs<-lapply(dfs, function(x) transform(x,id =interaction(Locus.ID,BP)))
 
  # select locus+bp (id) that are P=0.5 in at least mp populations
  desired<- lapply(dfs, function(x) ddply(x, .(id), function(x)
                                  if(sum(x$P==0.5)>=mp)x) )
  
  
  #put in a single df again
  desired<- do.call("rbind", desired)
  return(desired)
}

