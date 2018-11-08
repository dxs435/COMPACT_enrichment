#Taken from "MD_COMPACT_build_array.R"
#This is the sequence to restructure the gsets and create the discretized matrices

library(limma)
library(dplyr)
library(GEOquery)

#This will load all the functions
#load("pre_built.RData")
source("COMPACT_functions_v3.R")

i <- 1
Tlist = list()
#Genes are sparsely identified in GSE33785
for(f in list.files("Data/")){
  load(paste("Data/",f,sep=""))
  for(trmt in unique(gset$Treatment)){
    #need to check if we need a 0 time or control ("ctrl") point
    if(trmt == "ctrl"){
      next()
    }
    else if(min(gset$Time[gset$Treatment == trmt]) == 0){
      gse.tmp <- gset[,gset$Treatment == trmt]
      Tlist[[i]] <- gse.tmp
      #browser()
      names(Tlist)[i] <- paste(gsub(".RData","",f),paste(trmt),sep="_")
      i = i + 1
    }
    else if(any(gset$Time == 0)){
      gse.tmp <- gset[,(gset$Treatment == trmt | gset$Treatment == "ctrl")]
      gse.tmp$Treatment <- trmt
      Tlist[[i]] <- gse.tmp
      #browser()
      names(Tlist)[i] <- paste(gsub(".RData","",f),paste(trmt),sep="_")
      i = i + 1
    }
    #May need another 'else if' statement for experiments where the baseline is like 1 h...
    else{
      print("One or more GSEs are formatted incorrectly, or a baseline is needed")
      break()
    }
  }
}
#Tlist <- Tlist[-c(10,11)]
print("Finished treatment list")

discreteList <- list()
for(i in 1:length(Tlist)){
  disc.res <- eDcollapse(Tlist[[i]])
  disc.res <- eDcontrast(disc.res)
  disc.res <- eDdiscretize(disc.res)
  discreteList[[i]] <- disc.res
}
names(discreteList) <- names(Tlist)
print("Finished making dicrete time")
