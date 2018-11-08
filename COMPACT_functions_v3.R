#############################################
#   FUNCTIONS START
#   version 3 all functions for on-the-fly COMPACT calculation and removing the master array
#     -all fData commands were replaced by a saved list of fData objects--precomputed
#############################################

library(limma)
library(dplyr)

eDcollapse <- function(g){
  #browser()
  res.avg <- data.frame(matrix(nrow = nrow(fData(g)), ncol = length(unique(g$Time))))
  for(i in 1:length(unique(g$Time))){
    #Does rowmeans work on N=1?
    res.avg[,i] <- rowMeans(exprs(g[,g$Time == unique(g$Time)[i]]))
    colnames(res.avg)[i] <- paste(g$Treatment[1], as.character(unique(g$Time)[i]), sep = "//")
  }
  rownames(res.avg) <- rownames(fData(g))
  return(res.avg)
}

#This function needs a "0". In "GSE63742_sham operation", the first time is "2". I need to generally correct this, but
#for now, I'm just going to set the lowest time point as the 0 point.
eDcontrast <- function(eD){
  #browser()
  ct.matrix <- data.frame(matrix(nrow=nrow(eD), ncol = ncol(eD)-1))
  if(any(grepl("^0$",strsplit2(names(eD),"//")[,2]))){
    zeroi <- grep("^0$",strsplit2(names(eD),"//")[,2])
  }
  else{
    mintime <- min(as.integer(strsplit2(names(eD),"//")[,2]))
    zeroi <- grep(paste("^",mintime,"$",sep = ""),strsplit2(names(eD),"//")[,2])
  }
  zeroMatrix <- eD[,zeroi]
  edMatrix <- eD[,-zeroi]
  if(ncol(eD)>2){
    for(j in 1:ncol(edMatrix)){
      ct.matrix[,j] <- edMatrix[,j] - zeroMatrix 
      colnames(ct.matrix) <- colnames(edMatrix)
    }
  }
  else{
    ct.matrix[,1] <- edMatrix - zeroMatrix 
    names(ct.matrix) <- colnames(eD)[-zeroi]
  }
  
  rownames(ct.matrix) <- rownames(eD)
  return(ct.matrix)
}

eDdiscretize <- function(ctrst, minSet = 5, tn = 2){
  #browser()
  thresh <- log(tn)
  res.unary <- data.frame(matrix(nrow = nrow(ctrst), ncol = ncol(ctrst)))
  colnames(res.unary) <- colnames(ctrst)
  rownames(res.unary) <- rownames(ctrst)
  res.unary <- apply(ctrst,2,function(x) {(abs(x)>thresh)*sign(x)})
  if(ncol(res.unary) > 1){
    res.unary <- res.unary[,order(as.integer(strsplit2(colnames(res.unary),"//")[,2]))]  
  }
  #make a data frame for treatment-specific patterns
  #need to add something so I'll know what the time points are
  return(res.unary)
}

sortPatterns <- function(df){
  #browser()
  #if there is only one fold change...
  if(ncol(df) == 1){
    df.ordering <- order(df, decreasing = T)
    df <- as.data.frame(matrix(df[df.ordering,], dimnames = list(rownames(df)[df.ordering], colnames(df))))
    return(df)
  }
  #if there are multiple fold changes...
  else{
    p1 <- df+1
    p.decimal <- unname(apply(p1,1,function(x,d){sum(x*d^x)},seq(ncol(p1))))
    df <- df[order(p.decimal, decreasing = T),]
    df <- as.data.frame(df)
    return(df)
  }
}

getPatternGenes <- function(m, n){
  #browser()
  a<-discreteList[[m]]
  b<-discreteList[[n]]
  if(!any(strsplit2(colnames(a),"//")[,2] %in% strsplit2(colnames(b),"//")[,2])){
    #just do one contrast with the closest time
    timeOffset = 1000
    for(i in 1:ncol(a)){
      timeOffseti <- min(abs(as.integer(strsplit2(colnames(a),"//")[,2])[i]
                             - as.integer(strsplit2(colnames(b),"//")[,2])))
      if(timeOffseti < timeOffset){
        timeOffset <- timeOffseti
        adjust.a <- i
      }
    }
    adjust.b <- which(abs(as.integer(strsplit2(colnames(a),"//")[,2])[adjust.a] -
                as.integer(strsplit2(colnames(b),"//")[,2])) == timeOffset)
    colnames(b)[adjust.b] <- paste(strsplit2(colnames(b)[adjust.b],"//")[,1], "//",
                                   strsplit2(colnames(a)[adjust.a],"//")[,2], sep = "")
  }
  if(ncol(a) == 1 | ncol(b) == 1){
    #might need to be considerate of a situation where no gene changes twofold
    a.index<-which(strsplit2(colnames(a),"//")[,2] %in% strsplit2(colnames(b),"//")[,2])
    b.index<-which(strsplit2(colnames(b),"//")[,2] %in% strsplit2(colnames(a),"//")[,2])
    a.genes <- list()
    b.genes <- list()
    b.name <- colnames(b)[b.index]
    b <- as.data.frame(b[,b.index])
    colnames(b) <- b.name
    a.name <- colnames(a)[a.index]
    a <- as.data.frame(a[,a.index])
    colnames(a) <- a.name
    a <- sortPatterns(a)
    b <- sortPatterns(b)
    a.patterns <- unique(a)
    b.patterns <- unique(b)
    a.b.patterns <- intersect(unlist(a.patterns),unlist(b.patterns))
    for(i in 1:nrow(a.patterns)){
      probeIDs <- rownames(a)[a==a.patterns[i,]]
      a.genes[[i]] <- fData(Tlist[[m]])$Gene.symbol[match(probeIDs,fData(Tlist[[m]])$ID)]
      a.genes[[i]] <- a.genes[[i]][a.genes[[i]]!=""]
      names(a.genes)[i] <- as.character(a.patterns[i,])
    }
    a.genes <- a.genes[sort(names(a.genes))]
    for(i in 1:nrow(b.patterns)){
      probeIDs <- rownames(b)[b==b.patterns[i,]]
      #This needs to be changed for running on Gondolin
      b.genes[[i]] <- fData(Tlist[[n]])$Gene.symbol[match(probeIDs,fData(Tlist[[n]])$ID)] 
      b.genes[[i]] <- b.genes[[i]][b.genes[[i]]!=""]
      names(b.genes)[i] <- as.character(b.patterns[i,])
    }
    b.genes <- b.genes[sort(names(b.genes))]
  }
  if(ncol(a) > 1){
    a.index<-which(strsplit2(colnames(a),"//")[,2] %in% strsplit2(colnames(b),"//")[,2])
    b.index<-which(strsplit2(colnames(b),"//")[,2] %in% strsplit2(colnames(a),"//")[,2])
    a.genes <- list()
    b.genes <- list()
    a <- matrix(a[,a.index], nrow = nrow(a), ncol = length(a.index), dimnames = list(rownames(a), colnames(a)[a.index]))
    b <- matrix(b[,b.index], nrow = nrow(b), ncol = length(b.index),dimnames = list(rownames(b), colnames(b)[b.index]))
    a <- sortPatterns(a)
    b <- sortPatterns(b)
    #the next two lines run kind of slow
    a$pattern <- unlist(apply(a,1,paste,collapse=","), use.names = F)
    b$pattern <- unlist(apply(b,1,paste,collapse=","), use.names = F)
    #drop rare patterns (< 2)
    #I think there is an 'aggregate' function to this...
    a.patterns <- unique(a$pattern)
    b.patterns <- unique(b$pattern)
    a.b.patterns <- intersect(a.patterns,b.patterns)
    a <- a[(a$pattern %in% a.b.patterns),]
    b <- b[(b$pattern %in% a.b.patterns),]
    pattern.index <- names(table(c(a$pattern,b$pattern)))[table(c(a$pattern,b$pattern))>10]
    a <- a[(a$pattern %in% pattern.index),]
    b <- b[(b$pattern %in% pattern.index),]
    a.b.patterns <- unique(a$pattern)
    for(i in 1:length(a.b.patterns)){
      a.probeIDs <- rownames(a)[a$pattern == a.b.patterns[i]]
      a.genes[[i]] <- fData(Tlist[[m]])$Gene.symbol[match(a.probeIDs,fData(Tlist[[m]])$ID)]
      a.genes[[i]] <- a.genes[[i]][a.genes[[i]]!=""]
      names(a.genes)[i] <- i
      b.probeIDs <- rownames(b)[b$pattern == a.b.patterns[i]]
      b.genes[[i]] <- fData(Tlist[[n]])$Gene.symbol[match(b.probeIDs,fData(Tlist[[n]])$ID)]
      b.genes[[i]] <- b.genes[[i]][b.genes[[i]]!=""]
      names(b.genes)[i] <- i
    }
  }
  patternGenes <- list("a" = a.genes, "b" = b.genes, "patterns" = a.b.patterns)
  return(patternGenes)
}

#function to convert 3D COMPACT to 2D matrix of gene list lengths
flattenCOMPACT <- function(mat = COMPACT){
  #browser()
  mat[which(!is.na(mat))] <- 1
  #COMPACT <- (COMPACT == "1")*1
  COMPACT.flat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  for(i in 1:nrow(COMPACT.flat)){
    for(j in 1:ncol(COMPACT.flat)){
      #just trying to make it go faster...
      if(!is.na(mat[i,j,1])){
        COMPACT.flat[i,j] = sum(!is.na(mat[i,j,]))
      }
    }
  }
  return(COMPACT.flat)
}

binomialCoef <- function(n,k){
  #these numbers are too large to use factorial
  value <- lgamma(n+1)-((lgamma(k+1)+lgamma(n-k+1)))
  value <- exp(value)
  return(value)
}


hypergeometric.pvalue <- function(cf = COMPACT.flat, om = occurance.matrix, n = length(gene_list)){
  #browser()
  m.filter <- upper.tri(om, diag = F)*1
  n.filter <- lower.tri(om, diag = F)*1
  m.e <- sum(om * m.filter)
  print(paste("m.e=",m.e))
  n.e <- sum(om * n.filter)
  print(paste("n.e=",n.e))
  m.f <- sum(cf * m.filter)
  n.f <- sum(cf * n.filter)
  c.total <- sum(cf)
  total.prob = 0
  #get p-value for m
  #Should we sum over m.e or m.e - 1 ??
  if(m.e != 0){
    for(i in 0:(m.e-1)){
      hyper.prob = (binomialCoef(m.f, i)*binomialCoef(c.total-m.f, n-i))/binomialCoef(c.total, n)
      total.prob = total.prob + hyper.prob
    }
  }
  m.pvalue = 1 - total.prob
  
  total.prob = 0
  #get p-value for n
  if(n.e != 0){
    for(i in 0:(n.e-1)){
      hyper.prob = (binomialCoef(n.f, i)*binomialCoef(c.total-n.f, n-i))/binomialCoef(c.total, n)
      total.prob = total.prob + hyper.prob
    }
  }
  n.pvalue = 1 - total.prob
  #Ultimately these p-values will have to be corrected for multiple comparisons
  p.values = list("m" = m.pvalue, "n" = n.pvalue)
  return(p.values)
}

#i and j are the indices in the rotated matrix
getGeneNames <- function(m,n,i,j,mat){
  target_list <- masterArray[m,n,ncol(mat)+i-j,i,]
  #something like this may be more useful,
  genes <- sapply(gene_list, function(x) which(grepl(x,target_list)))
  print(names(genes)[which(genes > 0)])
  #print(gene_list[which(gene_list %in% target_list)])
}

#############################################
#   FUNCTIONS END // Psych?
#############################################