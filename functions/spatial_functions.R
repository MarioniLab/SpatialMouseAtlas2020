split.data.frame.by = function(df, by) {
  split.data.frame(df, df[, by])
}

updateMesenchymeLabel = function(mdt) {
  # mdt is a data.frame with factors containing Mesenchyme
  oldname = "Mesenchyme"
  newname = "Mixed mesenchymal mesoderm"
  
  # get names of columns in need of name change
  colsToChange = unique(which(mdt == oldname, arr = TRUE)[,2])
  for (i in colsToChange) {
    fac = mdt[,i]
    lvl = levels(fac)
    lvl[lvl == oldname] <- newname
    fac2 = factor(fac, levels = c(lvl, oldname))
    fac2[fac2 == oldname & !is.na(fac2)] <- newname
    fac2 = droplevels(fac2)
    mdt[,i] <- fac2
  }
  return(mdt)
}

rev_df = function(df){
  df[rev(seq_len(nrow(df))),]
}

sample_df = function(df, var, na.last = FALSE) {
  # var is a character of the column
  ind = which(is.na(df[,var]))
  ind_not = which(!is.na(df[,var]))
  if (na.last) {
    newind = c(sample(ind_not),ind)
  } else {
    newind = c(ind, sample(ind_not))
  }
  return(df[newind,])
}

add_scalebar = function(dist_um = 250, x = 2.75, y = -3.20, ...) {
  # this is a quantity to add to an existing ggplot
  
  # dist_um is the distance for the scalebar in um, default 250um
  # x and y are the coordinates to place the scalebar
  # usage: to add to a ggplot object like "g + add_scalebar()"
  
  # useful optional argument is box.col = "white"
  
  # need to make sure that ggplot of interest doesn't have group aesthetic
  # defined in the ggplot, but in the geom_polygon itself
  
  require(ggsn)
  add = ggsn::scalebar(location = "bottomright",
                       dist = dist_um/227.74, dist_unit = "units",
                       transform = FALSE,
                       x.min = x, x.max = x,
                       y.min = y, y.max = y,
                       height = 0.2,
                       box.fill = "black",
                       st.size = 0,
                       inherit.aes = FALSE,
                       ...)
  return(add)
}

reclusterGraph = function(graph, labels) {
  # graph is a graph object that has already been clustered
  # and the cluster labels are given as labels
  # this must be a named character vector
  # graph = snn
  # labels = clusters_sub
  
  # snn = buildSNNGraph(sce_filt_sub, use.dimred = "mbPCA_cosine_MNN")
  # clusters_graph = cluster_louvain(snn)
  # clusters_sub = paste0(cluster_val_name, ".", as.numeric(membership(clusters_graph)))
  # names(clusters_sub) <- sce_filt_sub$uniqueID
  
  labels_sub_all = c()
  
  for (label in sort(unique(labels))) {
    
    g_sub = induced.subgraph(graph, which(labels %in% label))
    clusters_graph_sub = cluster_louvain(g_sub)
    labels_sub = paste0(label, ".", as.numeric(membership(clusters_graph_sub)))
    names(labels_sub) <- names(labels)[which(labels %in% label)]
    
    labels_sub_all = c(labels_sub_all, labels_sub)
  }
  
  return(labels_sub_all[names(labels)])
}


reNormalise = function(sce_filt) {
  assay(sce_filt, "logcounts") <- logcounts(scater::logNormCounts(
    sce_filt, size_factors = sce_filt$total - counts(sce_filt)["Xist",]))
  return(sce_filt)
}

calculateEuclideanDistanceCells = function(sce,
                                           referenceCells,
                                           split = interaction(sce$embryo, sce$z),
                                           x = "x_global_affine",
                                           y = "y_global_affine") {
  # e.g.
  # referenceCells = colnames(sce)[sce$celltype_mapped == "Cardiomyocytes"]
  
  # select a set of cells that are at most r euclidean distance away from a
  # specified set of cells
  
  # input
  # sce object
  # grouping factor
  # x and y coords
  # subset of cells, character vector
  # radius
  
  # output
  # character vector of cells within the area
  require(class)
  
  sceList = splitSCE(sce, split)
  names(sceList) <- NULL
  
  distsList = lapply(sceList, function(sce_sub) {
    
    # use KNN to calculate the closest referenceCell, and then calculate
    # the euclidean distance of these to the nearest cell
    
    referenceCells_sub = intersect(colnames(sce_sub), referenceCells)
    if (length(referenceCells_sub) == 0) return(NA)
    
    out = class::knn(train = colData(sce_sub)[referenceCells_sub, c(x,y)],
                     test = colData(sce_sub)[, c(x,y)],
                     cl = factor(referenceCells_sub),
                     k = 1,
                     use.all = FALSE
    )
    pairs1 = colnames(sce_sub)
    pairs2 = as.character(out)
    
    coords1 = as.matrix(colData(sce_sub)[pairs1, c(x,y)])
    coords2 = as.matrix(colData(sce_sub)[pairs2, c(x,y)])
    
    d = sqrt(rowSums((coords1 - coords2)^2))
    names(d) <- colnames(sce_sub)
    
    return(d)
  })
  
  dists = unlist(distsList)
  dists <- dists[colnames(sce)]
  
  return(dists)
}

splitSCE = function(sce, split) {
  if (length(split) != ncol(sce)) stop("need split to be same length as ncol(sce)")
  
  out = sapply(unique(split), function(s){
    sce[,which(split == s)]
  }, simplify = FALSE)
  names(out) <- unique(split)
  return(out)
}

get_umap_graph <- function(umap_output) {
  # with thanks to Karsten Bach
  # umap_output is output from umap::umap() function
  
  require(Matrix)
  require(igraph)
  # extracts the graph information and outputs a weighted igraph
  dists.knn <- umap_output$knn[["distances"]]
  indx.knn <- umap_output$knn[["indexes"]]
  m.adj <- Matrix(0, nrow=nrow(indx.knn), ncol=nrow(indx.knn), sparse=TRUE)
  rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn)
  for (i in seq_len(nrow(m.adj))) {
    m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- dists.knn[i,]
  }
  igr <- graph_from_adjacency_matrix(m.adj, weighted=TRUE)
  return(igr)
}

get_MVS_vals = function(predmat) {
  MVS_vals_svm = t(apply(predmat,1,function(x) {
    # first column is the mapped celltype
    # second column is the runner-up celltype
    # mode_prop is the proportion of nearest neighbours with either mode1 or mode2
    # d is the nodepath distance for those celltypes
    
    if (all(is.na(x))) return(c(NA,NA,NA,NA))
    
    mode1 = getmode(x, 1:length(x))
    xx = x[x != mode1]
    if (length(xx) == 0) {
      mode2 = NA
      
    } else {
      mode2 = getmode(xx, 1:length(xx))
      if (is.null(mode2)) {
        mode2 <- NA
      }
      
    }
    return(c(mode1, mode2, mean(x %in% mode1), mean(x %in% c(mode1,mode2))))
  }))
  colnames(MVS_vals_svm) <- c("celltype_mapped",
                              "celltype_alternative",
                              "mapping_score",
                              "mapping_alternative_score"
  )
  return(MVS_vals_svm)
  
}



getMarkerArray = function(sce,
                          group,
                          assayName = "logcounts",
                          block = NULL,
                          subset = NULL,
                          verbose = TRUE,
                          pseudobulk = FALSE) {
  require(scran)
  
  # sce is a SingleCellExperiment object
  # group is a character for the group to split by
  # block is a character vector of the block to use
  # subset is a logical of the same length as ncol(sce) for which cells to
  # use for testing
  # verbose if TRUE
  # pseudobulk whether you perform pseudobulk by averaging values first
  
  if (!is.null(subset)) {
    sce_sub = sce[,subset]
  } else {
    sce_sub <- sce
  }
  
  genes = rownames(sce_sub)
  groups = sort(unique(colData(sce_sub)[, group]))
  
  if (!is.null(block)) {
    blockFactor = getInteractionFactor(sce_sub, levels = block)
  } else {
    blockFactor = NULL
  }
  
  marker_array = array(NA, dim = c(nrow(sce_sub), length(groups), 5),
                       dimnames = list(
                         rownames(sce_sub),
                         groups,
                         NULL
                       ))
  
  # per group level, output a matrix of genes and test outputs
  
  exprs = assay(sce_sub, assayName)
  
  for (grouplevel in groups) {
    
    if (verbose) print(grouplevel)
    
    grps = colData(sce_sub)[, group] == grouplevel
    
    if (pseudobulk) {
      
      colData(sce_sub)[, "grouplevel"] <- grps
      # then replace SCE with the pseudobulk one
      sce_pseudo <- sumCountsAcrossCells(assay(sce_sub, assayName),
                                         ids = colData(sce_sub)[,c("grouplevel", "grouplevel", block)],
                                         average = TRUE)
      blockFactor <- NULL
      
      exprs <- assay(sce_pseudo, "sum")
      grps = colData(sce_pseudo)[,"grouplevel"]
    }
    
    out = findMarkers(exprs,
                      groups = grps,
                      block = blockFactor,
                      test.type = "t",
                      direction = "any",
                      log.p = FALSE,
                      sorted = FALSE)
    out_TRUE <- out[["TRUE"]]
    
    meanIn = rowMeans(exprs[, grps, drop = FALSE])
    meanOut = rowMeans(exprs[, !grps, drop = FALSE])
    
    out_mat = cbind("p.value" = out_TRUE[, "p.value"],
                    "FDR" = out_TRUE[, "FDR"],
                    "LFC" = out_TRUE[, "logFC.FALSE"],
                    "meanIn" = meanIn,
                    "meanOut" = meanOut)
    
    marker_array[,grouplevel,] <- out_mat
    
  }
  
  dimnames(marker_array)[[3]] <- colnames(out_mat)
  
  return(marker_array)
  
  # e.g. marker_array[,"Gut",] to see all genes for this grouplevel
  # e.g. marker_array["Ttn",,] to see all p-values for this gene
}

markerArrayPlot = function(marker_array,
                           grouplevel,
                           otherGroupLabel = NULL,
                           LFC = 0.5,
                           p.value = 0.05,
                           FDR = NULL,
                           onlyLabelTopRanked = NULL) {
  # marker_array is the output from getMarkerArray
  # grouplevel is the 2nd dimension name to graph
  # onlyLabelTopRanked can be NULL, FALSE, or a number
  
  require(ggplot2)
  require(ggrepel)
  
  mk = as.data.frame(marker_array[,grouplevel,])
  if (is.null(FDR)) {
    mk$sig <- mk$p.value < p.value & abs(mk$LFC) > LFC
  } else {
    mk$sig <- mk$FDR < FDR & abs(mk$LFC) > LFC
  }
  mk$p.valueLFC <- mk$p.value
  mk$p.valueLFC[abs(mk$LFC) < LFC] <- 1
  mk$rank <- rank(mk$p.valueLFC)
  mk$lab <- rownames(mk)
  mk$gene <- rownames(mk)
  if (is.numeric(onlyLabelTopRanked)) {
    mk$lab[mk$rank > onlyLabelTopRanked] <- NA
  } else {
    mk$lab[!mk$sig] <- NA
  }
  
  if (!is.null(otherGroupLabel)) {
    gtitle = paste0("High LFC indicates higher expression in ", grouplevel,
                    ", compared to ", otherGroupLabel)
  } else {
    gtitle = paste0("High LFC indicates higher expression in ", grouplevel)
  }
  
  if (is.null(FDR)) {
    ltitle = paste0("Significantly DE (P-value < ", p.value,")")
  } else {
    ltitle = paste0("Significantly DE (FDR-adjusted P-value < ", FDR,")")
  }
  
  g = ggplot(mk, aes(x = LFC, y = -log10(p.value), label = gene)) +
    geom_point(aes(colour = sig)) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ylab("-log10(P-value)") +
    xlab("Log Fold Change") +
    ggtitle(gtitle) +
    guides(colour = guide_legend(title = ltitle,
                                 override.aes = list(size = 5))) +
    NULL
  
  if (!isFALSE(onlyLabelTopRanked)) {
    g = g + geom_text_repel(aes(label = lab))
  }
  
  if (is.null(FDR)) {
    g = g + geom_hline(yintercept = -log10(p.value), linetype = "dotted")
  }
  
  return(g)
}

markerArrayPlot2 = function(marker_array,
                            grouplevel,
                            otherGroupLabel = NULL,
                            LFC = 0.5,
                            p.value = 0.05,
                            FDR = NULL,
                            onlyLabelTopRanked = NULL,
                            sigColours = c("red","red")) {
  # marker_array is the output from getMarkerArray
  # grouplevel is the 2nd dimension name to graph
  # onlyLabelTopRanked can be NULL, FALSE, or a number
  
  require(ggplot2)
  require(ggrepel)
  
  mk = as.data.frame(marker_array[,grouplevel,])
  if (is.null(FDR)) {
    mk$sig <- mk$p.value < p.value & abs(mk$LFC) > LFC
  } else {
    mk$sig <- mk$FDR < FDR & abs(mk$LFC) > LFC
  }
  mk$p.valueLFC <- mk$p.value
  mk$p.valueLFC[abs(mk$LFC) < LFC] <- 1
  mk$rank <- rank(mk$p.valueLFC)
  mk$lab <- rownames(mk)
  mk$gene <- rownames(mk)
  if (is.numeric(onlyLabelTopRanked)) {
    mk$lab[mk$rank > onlyLabelTopRanked] <- NA
  } else {
    mk$lab[!mk$sig] <- NA
  }
  
  if (!is.null(otherGroupLabel)) {
    gtitle = paste0("High LFC indicates higher expression in ", grouplevel,
                    ", compared to ", otherGroupLabel)
  } else {
    gtitle = paste0("High LFC indicates higher expression in ", grouplevel)
  }
  
  if (is.null(FDR)) {
    ltitle = paste0("Significantly DE (P-value < ", p.value,")")
  } else {
    ltitle = paste0("Significantly DE (FDR-adjusted P-value < ", FDR,")")
  }
  
  g = ggplot(mk, aes(x = LFC, y = -log10(p.value), label = gene)) +
    geom_point(aes(fill = interaction(sig, LFC > 0),
                   colour = sig),
               shape = 21) +
    theme_classic() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("TRUE.TRUE" = sigColours[1],
                                 "TRUE.FALSE" = sigColours[2],
                                 "FALSE.TRUE" = "lightgrey",
                                 "FALSE.FALSE" = "lightgrey")) +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ylab("-log10(P-value)") +
    xlab("Log Fold Change") +
    ggtitle(gtitle) +
    guides(colour = guide_legend(title = ltitle,
                                 override.aes = list(size = 5))) +
    NULL
  
  if (!isFALSE(onlyLabelTopRanked)) {
    g = g + geom_text_repel(aes(label = lab))
  }
  
  if (is.null(FDR)) {
    g = g + geom_hline(yintercept = -log10(p.value), linetype = "dotted")
  }
  
  return(g)
}



getInteractionFactor = function(sce, levels) {
  # sce is a SingleCellExperiment object
  # levels is a character vector of the colData name to extract a factor for
  
  do.call(interaction,
          sapply(levels, function(x)
            colData(sce)[,x], simplify = FALSE))
}

getJointPCA = function(atlas, spatial,
                       assayNameAtlas = "cosineNorm",
                       assayNameSpatial = "cosineNorm",
                       inc_z = TRUE,
                       atlas_batch = "sample",
                       spatial_batch = c("embryo", "pos", "z"),
                       irlba = FALSE,
                       multibatch_2 = FALSE) {
  
  all_assay = cbind(assay(spatial, assayNameSpatial),
                    assay(atlas, assayNameAtlas)[rownames(spatial),])
  
  if (!inc_z) {
    spatial_batch <- setdiff(spatial_batch, "z")
  }
  
  spatial_batchFactor = do.call(interaction,
                                sapply(spatial_batch, function(x)
                                  colData(spatial)[,x], simplify = FALSE))
  
  atlas_batchFactor = do.call(interaction,
                              sapply(atlas_batch, function(x)
                                colData(atlas)[,x], simplify = FALSE))
  
  batchFactor = factor(c(as.character(spatial_batchFactor),
                         as.character(atlas_batchFactor)))
  print(table(batchFactor))
  
  merge.ordering = c(match(unique(as.character(atlas_batchFactor)), levels(batchFactor)),
                     match(unique(as.character(spatial_batchFactor)), levels(batchFactor)))
  print(merge.ordering)
  
  if (!irlba) {
    out = fastMNN(all_assay,
                  batch = batchFactor,
                  subset.row = !(rownames(all_assay) %in% "Xist"),
                  d = 50,
                  merge.order = merge.ordering)
    
    print(head(out@metadata$merge.info))
    
    joint_mnn_pca <- reducedDim(out, "corrected")
  } else {
    
    if (multibatch_2) {
      referenceBatchFactor = factor(ifelse(colnames(all_assay) %in% colnames(atlas), "reference", "query"))
      
      print(table(referenceBatchFactor))
      
      pca_all = multiBatchPCA(all_assay,
                              batch = referenceBatchFactor,
                              weights = c("query" = 0, "reference" = 1),
                              subset.row = !(rownames(all_assay) %in% "Xist"),
                              d = 50,
                              preserve.single = TRUE)[[1]]
      
    } else {
      pca_all <- irlba::prcomp_irlba(t(all_assay), n = 50)$x
      rownames(pca_all) <- colnames(all_assay)
    }
    
    out = reducedMNN(pca_all,
                     batch = batchFactor,
                     # subset.row = !(rownames(all_assay) %in% "Xist"),
                     # d = 50,
                     merge.order = merge.ordering)
    
    joint_mnn_pca <- out$corrected
    
  }
  
  
  return(joint_mnn_pca)
}

# function to predict within the joint PCA embedding usng SVM
getSVMPrediction = function(atlas, spatial,
                            joint_mnn_pca,
                            outcome = "celltype_parsed",
                            rebalance = FALSE) {
  
  require(e1071)
  
  X_spatial = joint_mnn_pca[colnames(spatial), , drop = FALSE]
  
  if (nrow(X_spatial) == 0) return(NULL)
  
  X = joint_mnn_pca[colnames(atlas),]
  y = factor(colData(atlas)[,outcome])
  
  if (length(unique(y)) == 1) {
    print("only one class in reference, assigning all query to that class")
    
    return(rep(unique(y), nrow(X_spatial)))
    
  }
  
  if (rebalance) {
    # repeatedly subset (with min provided) each class and fit the model
    # get majority rules from each model fit
    
    nval = max(c(min(table(y)), 100))
    n_rebalance = 50
    
    pred.celltypesList = sapply(1:n_rebalance, function(i){
      
      print(paste0("rebalance... ", i, " of ", n_rebalance))
      
      y_ind = unlist(tapply(seq_along(y), y, function(x) {
        if (length(x) < nval) return(x)
        sample(x, nval)
      },simplify = FALSE))
      
      X_sub = X[y_ind, ]
      y_sub = y[y_ind]
      
      fit = svm(x = X_sub, y = y_sub)
      
      pred.celltypes = predict(fit, X_spatial)
      
      return(pred.celltypes)
    }, simplify = FALSE)
    pred.celltypesMat = do.call(cbind, lapply(pred.celltypesList, as.character))
    
    pred.celltypes = apply(pred.celltypesMat, 1, getmode, 1:ncol(pred.celltypesMat))
    
  } else {
    
    # normal SVM for all cells
    
    fit = svm(x = X, y = y)
    
    print("fitted model")
    
    pred.celltypes = predict(fit, X_spatial)
    
  }
  
  print("predicted")
  
  return(pred.celltypes)
  
}

vectorSubset = function(vec, mat){
  # vec is a named vector
  # mat is a matrix containing the names or indices for which you want
  # to get the entries of vec
  
  vmat = c(mat)
  vvec = vec[vmat]
  
  vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
  colnames(vecmat) <- colnames(mat)
  rownames(vecmat) <- rownames(mat)
  
  return(vecmat)
}

vectorMatch = function(vec, mat, vecnames){
  # vec is an unnamed vector
  # vecnames is the names of vec
  # mat is a matrix containing the names or indices for which you want
  # to get the entries of vec, matching vecnames
  
  vmat = c(mat)
  
  vecind = match(vmat,vecnames)
  
  vvec = vec[vecind]
  
  vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
  colnames(vecmat) <- colnames(mat)
  rownames(vecmat) <- rownames(mat)
  
  return(vecmat)
}

expressionDotPlot = function(SCE, genes, group, assayName = "logcounts") {
  # SCE single cell experiment object
  # genes character vector of genes corresponding to rownames of SCE
  # group is the grouping taken from colData(SCE)
  # get expression dotplot for these genes and celltypes
  
  propExpressed = apply(assay(SCE, assayName)[genes,],1,function(x){
    tapply(x, colData(SCE)[,group], function(x) mean(x>0))
  })
  
  genes_ordered = genes[hclust(dist(t(propExpressed)))$order]
  
  meanExpressed = apply(assay(SCE, assayName)[genes,],1,function(x){
    tapply(x, colData(SCE)[,group], function(x) mean(x))
  })
  
  meanExpressedNorm <- t(t(meanExpressed)/apply(meanExpressed,2,max))
  
  expressedDF = cbind(melt(propExpressed),
                      melt(meanExpressed)[,3],
                      melt(meanExpressedNorm)[,3])
  colnames(expressedDF) <- c("celltype", "gene", "proportion", "mean", "mean_norm")
  expressedDF$gene <- factor(expressedDF$gene, levels = genes_ordered)
  expressedDF$celltype <- factor(expressedDF$celltype)
  
  g = ggplot(expressedDF, aes(x = celltype, y = gene)) +
    geom_point(aes(size = proportion, fill = celltype,
                   alpha = mean_norm, stroke = mean_norm), colour = "black", shape = 21) +
    theme_classic() +
    theme(legend.position = "right") +
    #  scale_colour_manual(values = celltype_colours, aesthetics = c("fill", "colour")) +
    scale_alpha_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
    scale_size_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") +
    ylab("") +
    guides(fill = FALSE) +
    guides(size = guide_legend(title = "Proportion of cells\nexpressing")) +
    guides(alpha = guide_legend(title = "Mean standardised\nexpression")) +
    #ggtitle("Pijuan Sala et al data") +
    NULL
  print(g)
  return(g)
}

# function to rotate the embryos
rotateDF = function(DF,
                    xname = "segmentation_vertices_x_global_affine",
                    yname = "segmentation_vertices_y_global_affine",
                    ang = 0) {
  # ang is a numeric vector named three values corresponding to embryo1, embryo2, and embryo3
  
  ang_long = ang[as.character(DF$embryo)]
  ang_rad = ang_long/180
  
  x = DF[,xname]
  y = DF[,yname]
  
  x_turn = x*cos(ang_rad) - y*sin(ang_rad)
  y_turn = x*sin(ang_rad) + y*cos(ang_rad)
  
  # reset the columns and then return the DF
  
  DF[,xname] <- x_turn
  DF[,yname] <- y_turn
  
  return(DF)
  
}

getAreaFactorsDirectly = function(sce, transform = NULL) {
  # sce is a SingleCellExperiment object
  # transform is a function
  meta = colData(sce)
  counts = assay(sce, "counts")
  if (!is.null(transform)) {
    sizeFactors.area <- transform(meta$Area)/mean(transform(meta$Area))
  } else {
    sizeFactors.area <- meta$Area/mean(meta$Area)
  }
  logcounts.area <- log2(t(t(counts)/sizeFactors.area) + 1)
  return(list(sizeFactors = sizeFactors.area,
              logcounts.area = logcounts.area))
}

filterAreas = function(areas) {
  # areas is a numeric vector, output is a logical vector giving TRUE (keep)
  # or FALSE (remove)
  
  sqrt_centre = median(sqrt(areas))
  sqrt_scale = mad(sqrt(areas))
  
  pnorm.fit = pnorm(sqrt(areas), mean = sqrt_centre,
                    sd = sqrt_scale, lower.tail = F)
  area_thresh.sqrt = min(sqrt(areas)[which(p.adjust(pnorm.fit, method = "BH") < 0.01)])
  area_thresh = area_thresh.sqrt^2
  
  return(areas < area_thresh)
}

getmode <- function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

matvec2numeric = function(vec) {
  # vec is a single character object, that's in array format from matlab
  # e.g. "[1462 1462 1525 1525]"
  as.numeric(unlist(strsplit(gsub("\\[|\\]","",vec), " ")))
}

spatialPlotPointsGene = function(gene,
                                 sce,
                                 assayname = "logcounts",
                                 save = FALSE,
                                 save.dir = NULL,
                                 all = TRUE,
                                 binary = FALSE) {
  d = data.frame(cbind(colData(sce), "gene" = assay(sce, assayname)[gene,]))
  if (!all) {
    d <- subset(d, zstack == "Z1")
  }
  d <- sort_df(d, "gene")
  g = ggplot(d,
             aes(x = x_global, y = -y_global)) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "bottom") +
    ggtitle(gene) +
    NULL
  
  gs_gradient = scale_colour_gradient2(low = "grey", mid = "yellow", high = "red")
  gs_manual = scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey"))
  
  gp_gradient = geom_point(aes(colour = gene), size = 0.25)
  gp_manual = geom_point(aes(colour = gene > 0), alpha = 1, size = 1)
  
  if (binary) {
    g = g + gp_manual + gs_manual
  } else {
    g = g + gp_gradient + gs_gradient
  }
  
  if (save) ggsave(g, file = paste0(save.dir,gene,"_allcoords.png"), height = 4, width = 12)
  return(g)
}


getSegmentationVerticesDF = function(DF,
                                     xname = "segmentation_vertices_x_global",
                                     yname = "segmentation_vertices_y_global",
                                     othercols = c("uniqueID","z")) {
  # DF is a DataFrame object
  # othercols is the others to keep
  
  long_x = unlist(DF[,xname])
  long_y = unlist(DF[,yname])
  
  if (length(long_x) != length(long_y)) stop("x and y need to be same length")
  
  long_xy = data.frame(
    long_x,
    long_y
  )
  colnames(long_xy) <- c(xname, yname)
  
  long_DF = cbind(
    rep(DF[,othercols], times = unlist(lapply(DF[,xname], length))),
    long_xy
  )
  
  return(as.data.frame(long_DF))
}


processSeqFISH = function(params) {
  
  list2env(params, globalenv())
  
  #   params = list(
  #   fov_info_file = fov_info_file,
  #   mRNA_files_dir = mRNA_files_dir,
  #   cell_files_dir = cell_files_dir,
  #   segmentation_files_dir = segmentation_files_dir,
  #   processed_dir = processed_dir,
  #   processed_suffix = processed_suffix,
  #   fov_info = fov_info,
  #   mRNA_files = mRNA_files,
  #   cell_files = cell_files,
  #   embryo_label = embryo_label,
  #   Estage_label = Estage_label,
  #   allgenes = allgenes
  # )
  
  mRNA_sceList = list()
  mRNA_dfList = list()
  
  for (mRNA_file in mixedsort(mRNA_files)) {
    
    print(mRNA_file)
    # write("started!", "mRNA_files_completed.txt", append = TRUE)
    
    pos = grep("Pos", unlist(strsplit(mRNA_file, "-")), value = TRUE)[1]
    fov = paste0("fov.", as.numeric(gsub("Pos","",pos))+1)
    embryo = embryo_label[fov]
    Estage = Estage_label[fov]
    
    cell_file = grep(paste0(pos, "-"), cell_files, value = TRUE)[1]
    
    mRNA_df_raw = read.csv(mRNA_file, header = TRUE)
    head(mRNA_df_raw)
    
    if ("inensity" %in% colnames(mRNA_df_raw)) {
      cnames = colnames(mRNA_df_raw)
      cnames[cnames == "inensity"] <- "intensity"
      colnames(mRNA_df_raw) <- cnames
    }
    
    dim(mRNA_df_raw)
    length(unique(mRNA_df_raw$cellID))
    length(unique(mRNA_df_raw$geneID))
    
    mRNA_df_raw$geneID = factor(mRNA_df_raw$geneID, levels = allgenes)
    mRNA_df_raw$pos = pos
    mRNA_df_raw$fov = fov
    mRNA_df_raw$embryo = embryo
    mRNA_df_raw$Estage = Estage
    
    length(unique(mRNA_df_raw$z))
    hist(mRNA_df_raw$intensity, 100, xlab = pos)
    hist(mRNA_df_raw$x, 100, xlab = pos)
    hist(mRNA_df_raw$y, 100, xlab = pos)
    
    sort(table(mRNA_df_raw$geneID))
    
    # cell info
    cell_info = read.csv(cell_file, header = TRUE)
    head(cell_info)
    
    # give a unique name
    uniqueID = paste0(
      embryo, "_",
      pos, "_",
      "cell", cell_info$cellID, "_",
      "z", cell_info$z
    )
    
    cell_info$uniqueID <- uniqueID
    rownames(cell_info) <- uniqueID
    
    uniqueID_mRNA = paste0(
      mRNA_df_raw$embryo, "_",
      mRNA_df_raw$pos, "_",
      "cell", mRNA_df_raw$cellID, "_",
      "z", mRNA_df_raw$z
    )
    
    all(uniqueID_mRNA %in% uniqueID) # this must be TRUE
    mRNA_df_raw$uniqueID <- uniqueID_mRNA
    
    fov_info_sub = fov_info[fov,,drop = FALSE]
    
    bounds_x = as.numeric(fov_info_sub[,paste0("bound_x_",1:4)])
    bounds_y = as.numeric(fov_info_sub[,paste0("bound_y_",1:4)])
    
    bound_x = range(bounds_x)
    bound_y = range(bounds_y)
    
    subcoords = mRNA_df_raw[,c("x","y")]
    
    range_x = c(0,2048)
    range_y = c(0,2048)
    
    newcoords_x = (subcoords[,"x"] - range_x[1])*(bound_x[2] - bound_x[1])/range_x[2] + bound_x[1]
    newcoords_y = (subcoords[,"y"] - range_y[1])*(bound_y[2] - bound_y[1])/range_y[2] + bound_y[1]
    
    cell_info$x_global = (cell_info$Centroid_1 - range_x[1])*(bound_x[2] - bound_x[1])/range_x[2] + bound_x[1]
    cell_info$y_global = (cell_info$Centroid_2 - range_y[1])*(bound_y[2] - bound_y[1])/range_y[2] + bound_y[1]
    
    mRNA_df_raw$x_global = newcoords_x
    mRNA_df_raw$y_global = newcoords_y
    
    mRNA_counts_arr = tapply(mRNA_df_raw$intensity,
                             list(mRNA_df_raw$uniqueID,
                                  mRNA_df_raw$geneID),
                             function(x) length(x))
    mRNA_counts_arr[is.na(mRNA_counts_arr)] <- 0
    
    mRNA_counts = apply(mRNA_counts_arr, length(dim(mRNA_counts_arr)), c)
    
    # only keep rows with corresponding cell info, although this should all be TRUE
    mRNA_counts = mRNA_counts[rownames(mRNA_counts) %in% rownames(cell_info),]
    
    coords = cell_info[rownames(mRNA_counts),]
    coords$pos = pos
    coords$fov = fov
    coords$embryo = embryo
    coords$Estage = Estage
    
    # grab boundary segmentation
    segmentation_files = list.files(segmentation_files_dir,
                                    pattern = paste0("boundaryseg_", pos,"-"),
                                    full.names = TRUE)
    zvals = unlist(lapply(unlist(lapply(lapply(strsplit(segmentation_files, pos), "[",2), strsplit, "-|_"), recursive = FALSE), "[",2))
    names(segmentation_files) <- zvals
    
    # only keep zvals that appear in our count data
    zvals <- intersect(zvals, sprintf("%04.f",unique(mRNA_df_raw$z)))
    
    boundary_polygons_zvals = list()
    
    for (i in zvals) {
      
      print(i)
      
      boundary_polygons = extractPolygonsFromSegmentationH5(
        segmentation_files[i],
        modes = c("boundaryPolygons"),
        cellFilterFactor = 5,
        minPixels = 50,
        cellPrefix = paste0(embryo, "_", pos, "_cell"),
        cellSuffix = paste0("_z", as.numeric(i)),
        plot = FALSE,
        verbose = FALSE
      )
      colnames(boundary_polygons) <- c("x","y","cell","uniqueID")
      
      boundary_polygons <- pruneStraightLines(boundary_polygons, group = "uniqueID")
      
      boundary_polygons$z <- as.numeric(i)
      
      boundary_polygons$x_global = (boundary_polygons$x - range_x[1])*(bound_x[2] - bound_x[1])/(range_x[2]/4) + bound_x[1]
      boundary_polygons$y_global = (boundary_polygons$y - range_y[1])*(bound_y[2] - bound_y[1])/(range_y[2]/4) + bound_y[1]
      
      boundary_polygons_zvals[[i]] <- boundary_polygons
    }
    
    boundary_polygons <- do.call(rbind, boundary_polygons_zvals)
    
    g = ggplot(boundary_polygons, aes(x = x, y = -y)) +
      geom_polygon(aes(fill = uniqueID), show.legend = FALSE, colour = "black") +
      theme_classic() +
      facet_wrap(~z)
    print(g)
    
    mRNA_sce = SingleCellExperiment(assays = list(counts = t(mRNA_counts)),
                                    colData = DataFrame(coords)
    )
    
    colData(mRNA_sce)$segmentation_vertices_x = IRanges::NumericList(split(as.numeric(boundary_polygons$x),
                                                                           factor(boundary_polygons$uniqueID, levels = rownames(colData(mRNA_sce)))
    ))
    
    colData(mRNA_sce)$segmentation_vertices_y = IRanges::NumericList(split(as.numeric(boundary_polygons$y),
                                                                           factor(boundary_polygons$uniqueID, levels = rownames(colData(mRNA_sce)))
    ))
    
    colData(mRNA_sce)$segmentation_vertices_x_global = IRanges::NumericList(split(as.numeric(boundary_polygons$x_global),
                                                                                  factor(boundary_polygons$uniqueID, levels = rownames(colData(mRNA_sce)))
    ))
    
    colData(mRNA_sce)$segmentation_vertices_y_global = IRanges::NumericList(split(as.numeric(boundary_polygons$y_global),
                                                                                  factor(boundary_polygons$uniqueID, levels = rownames(colData(mRNA_sce)))
    ))
    
    g = ggplot(getSegmentationVerticesDF(colData(mRNA_sce)),
               aes(x = segmentation_vertices_x_global,
                   y = -segmentation_vertices_y_global)) +
      geom_polygon(aes(group = uniqueID), fill = "grey", colour = "black") +
      facet_wrap(~z) +
      theme_classic() +
      NULL
    print(g)
    
    mRNA_sceList[[mRNA_file]] <- mRNA_sce
    mRNA_dfList[[mRNA_file]] <- mRNA_df_raw
    
    # write(mRNA_file, "mRNA_files_completed.txt", append = TRUE)
    
  }
  
  out = list(mRNA_sceList = mRNA_sceList,
             mRNA_dfList = mRNA_dfList)
  
  return(out)
}
