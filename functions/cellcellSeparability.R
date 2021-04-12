

cellCellSeparability = function(sce, split, group, coordNames, plot = FALSE) {
  
  ############## boundary sharpness measure
  # cell cell separability measure
  # i.e. how well can SVM distinguish between these two groups in space?
  # averaged per embryo and z-stack
  # need to finish writing this function, and then save in neighbourhood.R
  
  # for a pair of cell_sets, calculate their boundary sharpness
  # require(e1071)
  
  # plot = FALSE
  # group = "cluster_joint_sub"
  # split = interaction(sce$embryo, sce$z)
  # coordNames = c("x_global_affine", "y_global_affine")
  
  sceSplit = splitSCE(sce, split)
  
  subgroups = sort(unique(colData(sce)[, group]))
  dist_array = array(NA, dim = c(length(subgroups), length(subgroups), length(sceSplit)),
                     dimnames = list(subgroups, subgroups, names(sceSplit)))
  
  for (sceSplitName in names(sceSplit)) {
    
    sce_sub = sceSplit[[sceSplitName]]
    
    for (ctype1 in subgroups) {
      
      sce_set1 = sce_sub[, colData(sce_sub)[,group] == ctype1]
      
      for (ctype2 in subgroups) {
        
        print(sceSplitName)
        print(ctype1)
        print(ctype2)
        
        if (ctype1 == ctype2) {
          dist_array[ctype1,ctype2,sceSplitName] <- 0
          next
        }
        
        sce_set2 = sce_sub[, colData(sce_sub)[,group] == ctype2]
        
        if (ncol(sce_set1) <= 1 | ncol(sce_set2) <= 1) {
          dist_array[ctype1,ctype2,sceSplitName] <- NA
          next
        }
        
        if (!is.na(dist_array[ctype2,ctype1,sceSplitName])) {
          dist_array[ctype1,ctype2,sceSplitName] <- dist_array[ctype2,ctype1,sceSplitName]
          next
        }
        
        coords_sub = colData(sce_sub)[,coordNames]
        coords1 = colData(sce_set1)[,coordNames]
        coords2 = colData(sce_set2)[,coordNames]
        
        coords = as.data.frame(rbind(cbind(coords1, group = ctype1),
                                     cbind(coords2, group = ctype2)))
        coords$group <- factor(coords$group)
        
        fit = svm(x = coords[,1:2], y = coords[,3])
        resub = predict(fit, coords[,1:2])
        table(resub, coords[,3])
        # balanced accuracy, higher value means they are more distinct
        # thus can directly be treated as a distance
        bal = mean(tapply(coords[,3] == resub, coords[,3], mean))
        bal
        
        if (plot) {
          plot(coords_sub[,1], -coords_sub[,2], col = factor(predict(fit, coords_sub)), pch = 16, asp = 1)
          points(coords1[,1], -coords1[,2], col = "blue", pch = 16, asp = 1)
          points(coords2[,1], -coords2[,2], col = "yellow", pch = 16, asp = 1)
        }
        
        dist_array[ctype1,ctype2,sceSplitName] <- bal
        
      }
    }
  }
  
  dist_mat = apply(dist_array, c(1,2), mean, na.rm = TRUE)
  dist_mat[is.nan(dist_mat)] <- NA
  
  dist_mat
  range(dist_mat[lower.tri(dist_mat)], na.rm = TRUE)
  
  dist_mat_mask = dist_mat
  dist_mat_mask[is.na(dist_mat_mask)] <- 1
  diag(dist_mat_mask) <- 0.5
  
  return(dist_mat_mask)
}



cellCellSeparabilityMap = function(dist_mat, ...) {
  dist_mat_out = list(obs = dist_mat, pmat = dist_mat)
  out = cellCellContactMap(dist_mat_out, ...) +
    scale_fill_gradient2(low = "#ed6495", high = "white", na.value = "white",
                         midpoint = median(c(dist_mat)),
                         limits = c(0.5, 1),
                         breaks = c(0.5, 1),
                         labels = c("Not Separable", "Separable")) +
    NULL
  return(out)
}
