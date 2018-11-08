#############################################
# Marker enrichment analysis through combinatorial compact
# Author: DGS
# 
# This all assumes that the data objects from "pre_built.RData" are in the environment 
#
# VERSION 3: major changes to the algorithm--importantly the master array is dropped
#   1. COMPACT matrices will be calculated on the fly
#   2. COMPACT time codes will be sorted differently (using Raj's method)
#   3. Enrichment is scored as a p-value calculated from a hypergeometric distribution
#
# General procedure
#   1. Read in the gene list
#   2. Build COMPACT gene lists
#   3. Get the gene list / COMPACT intersections
#   4. calculate p-values
#
#############################################
library(reshape2)
library(ggplot2)
library(gplots)

#read in the gene/marker list to tested
gene_list <- read.table(file = "HSC_test_list.txt", colClasses = "character")
gene_list <- gene_list$V1


#make empty results matrix
enrichMatrix <- data.frame(matrix(0, nrow = length(Tlist), ncol = length(Tlist)), row.names = names(Tlist))
colnames(enrichMatrix) <- names(Tlist)

pairs_for_CB <- combn(length(discreteList),2)

for(p in 1:ncol(pairs_for_CB)){
  #browser()
  ind <- pairs_for_CB[,p]
  print(paste("working on pairings: ", ind[1], ind[2]))
  #make 3D array of COMPACT and individual gene lists
  patternLists <- getPatternGenes(ind[1], ind[2])
  #make an empty 3D array to store temporary COMPACT and gene lists
  #...getting the dimensions
  maxA <- max(sapply(patternLists[["a"]],length, USE.NAMES = F))
  maxB <- max(sapply(patternLists[["b"]],length, USE.NAMES = F))
  maxGenes <- max(c(maxA,maxB))
  COMPACTdims <- c(length(patternLists[["patterns"]]), maxGenes)
  #...making the empty array
  COMPACT <- array(NA, dim =  c(COMPACTdims[1], COMPACTdims[1], COMPACTdims[2]))
  for(i in 1:nrow(COMPACT)){
    for(j in 1:ncol(COMPACT)){
      genes <- intersect(unique(patternLists[["a"]][[i]]), unique(patternLists[["b"]][[j]]))
      if(length(genes) > 0){
        COMPACT[i,j,1:length(genes)] <- genes
      }
    }
  }
  #make an empty matrix to store the Marker occurance in gene lists
  occurance.matrix <- matrix(0, nrow = COMPACTdims[1], ncol = COMPACTdims[1])
  for(i in 1:nrow(COMPACT)){
    for(j in 1:ncol(COMPACT)){
      occurance.matrix[i,j] <- sum(sapply(gene_list, function(x) sum(grepl(x, COMPACT[i,j,]))))
    }
  }
  #flatten the COMPACT array
  COMPACT.flat <- flattenCOMPACT()
  #calculate enrichment p-value
  pvalues <- hypergeometric.pvalue()
  enrichMatrix[ind[1], ind[2]] <- pvalues[["m"]]
  enrichMatrix[ind[2], ind[1]] <- pvalues[["n"]]
}

#make diagonal NAs, I'm concerned that the 0's will influence FDR adjustments
ematrix <- as.matrix(enrichMatrix)
ematrix[ematrix == 0] <- NA
#get rid of negative p-values...also, why are there negative p-values??
if(any(enrichMatrix < 0)){
  ematrix[which(ematrix < 0)] <- NA
}

#Correcting for multiple comparisons
#FDR by ematrix row
enrich.adjusted.pvalue <- vector(mode="numeric")
# for(trmt in 1:nrow(ematrix)){
#   enrich.adjusted.pvalue <- c(enrich.adjusted.pvalue,unname(p.adjust(ematrix[trmt,], method="BH")))
# }

#FDR by whole matrix
enrich.adjusted.pvalue <- p.adjust(ematrix, method="BH")
enrich.adjusted.pvalue[is.na(enrich.adjusted.pvalue)] <- 1

#The p-values seemed similar between by row and by whole matrix, so I'll continue with whole matrix bc it's less lines

#putting the p-values in a matrix
#this method of putting a vector into a matrix only works for square matrices, which I think we will always have
#For FDR by row, the matrix needs to be transformed
enrich.matrix.adjusted <- data.frame(matrix(enrich.adjusted.pvalue, nrow = length(Tlist), ncol = length(Tlist)), row.names = names(Tlist))
colnames(enrich.matrix.adjusted) <- names(Tlist)

#visualizing...
#because there is such a huge range of values and we only care about the low values, the values above 5e-3 are forced to 5e-3
ematrix <- as.matrix(enrich.matrix.adjusted)
ematrix[ematrix > 5e-3] <- 5e-3

heatmap.2(-log(t(ematrix)),
          Colv = F,
          dendrogram = "row",
          #scale = "none",
          trace = "none",
          key = T,
          cexCol = 1,
          cexRow = 1,
          margins = c(5,15))


save(ematrix, file="HSC_enrichment_v3.RData")


