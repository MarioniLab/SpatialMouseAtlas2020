get_k.mapped = function(joint_pcs, queryCells, referenceCells, ...) {
  # joint_pcs is a matrix of joint PC coordinates
  # referenceCells is character vector (corresponding to rownames of joint_pcs)
  # queryCells is character vector(corresponding to rownames of joint_pcs)
  # ... passes to BiocNeighbors::queryKNN
  require(BiocNeighbors)
  
  knns = BiocNeighbors::queryKNN(
    joint_pcs[referenceCells,],
    joint_pcs[queryCells,],
    get.index = TRUE,
    get.distance = FALSE,
    k = 10,
    ...)
  
  k.mapped = vectorSubset(referenceCells, knns$index)
  rownames(k.mapped) <- queryCells
  
  return(k.mapped)
}


weightedMeanMatrix = function(E, D, sub = NULL) {
  # calculate a matrix of weighted means
  # where E is a genes x cells expression matrix
  # D is a cells x cells distance matrix (larger means more far apart)
  # and sub is a index or character vector to subset the cells by
  
  if (!is.null(sub)) {
    return(t(t(E[,sub, drop = FALSE] %*% (1/D[sub,sub, drop = FALSE])) / colSums(1/D[sub,sub, drop = FALSE])))
  }
  return(t(t(E %*% (1/D)) / colSums(1/D)))
}

genesetGOtest = function(set, universe, termsList) {
  # set is the character vector of genes in geneset
  # universe is character vector of genes to include in universe
  # termsList is a list of the pathways
  
  termsListfiltered = lapply(termsList, function(x) return(intersect(universe, x)))
  keepForTesting = unlist(lapply(termsListfiltered, length)) >= 8 & unlist(lapply(termsListfiltered, length)) <= 500
  
  pval = sapply(1:length(termsList), function(i){
    
    if (!keepForTesting[i]) return(1)
    
    termGenes = termsListfiltered[[i]]
    return(fisher.test(table(universe %in% set, universe %in% termGenes), alt = "g")$p.value)
  })
  names(pval) <- names(termsListfiltered)
  return(pval)
}


diagonalise = function(mat) {
  # reorder rows and columns so that it's closest to diagonal as possible
  # mat = table(commA[commonGenes],commB[commonGenes])
  newmat = mat
  mis = NULL
  # i = 1
  while (sum(newmat)>0) {
    nonzeromin = min(newmat[newmat>0])
    bestind = which(newmat==nonzeromin,arr=TRUE)
    bestrow = bestind[1,1]
    bestcol = bestind[1,2]
    newmat = newmat[c(bestrow,setdiff(1:nrow(newmat),bestrow)),
                    c(bestcol,setdiff(1:ncol(newmat),bestcol))]
    mis = c(mis,sum(newmat[2:nrow(newmat),1],
                    newmat[1,2:ncol(newmat)]))
    newmat[1,1] <- 0
    
    # print(sum(newmat))
  }
  
  matOrdered = mat[rownames(newmat),colnames(newmat)]
  return(matOrdered)
}
