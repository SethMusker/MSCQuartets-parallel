library(MSCquartets)
library(future.apply)
plan(multiprocess,gc=TRUE)

source("Parallel_NANUQ.R")

taxa<-read.table("mytaxa.txt") # taxon names in one column w/out header
taxa<-as.vector(taxa$V1)

all_nanuq<-Parallel_NANUQ(genedata="genetrees.tre",
                 taxanames = taxa,
                 outfile = "genetrees_NANUQ.nex",
		 RAM_Gigs=5)


