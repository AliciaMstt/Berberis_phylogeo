tsv2Fasta <- function(tsvfile, outputfile, len){
  # This function takes a tsv file produced by the exporting tool of Stacks, or from a modification of it, and exports the
  # consensus sequences as Fasta format
  # tsv = path to the tsv .SNP.SNPs file 
  # outputfile= path to the desired outputfile
  # len = corresponds to the colw term of the write.dna function, ie a numeric specifying the number of nucleotides per column
  
  ##### extract consensus sequence from .SNP.SNP tvs file generated by Stacks
  final = read.delim(paste(tsvfile, sep = ""), header = T) 
  x <- cbind(final$CatalogID, as.vector(final$Consensus))
  rownames(x) <- x[,1] 
  x<-x[,2]
  
  ## Transform to fasta
  require(ape)
  write.dna(x, file=outputfile, format="fasta", colw=len)
}
