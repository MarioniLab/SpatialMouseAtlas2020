pruneStraightLines = function(df, minRadians = 0, groups = NULL){
  # df is a data.frame that must have columns x and y
  # and should be grouped by cell
  # the difference is angles will range from -pi to pi
  
  if (!is.null(groups)) {
    groups_vals = df[,groups]
    groups_diff = groups_vals[1:(length(groups_vals)-1)] != groups_vals[2:length(groups_vals)]
  }
  
  return(df[c(1, which(
    c(groups_diff[2:length(groups_diff)]) | 
      (abs(diff(atan(diff(df$y)/diff(df$x)))) > minRadians)
  )+1),])
  
}



# function that takes in the segmentation filename and outputs a list
# of dataframes, one for chull and another for boundary points
# with a flag for plotting

extractPolygonsFromSegmentationH5 = function(filename,
                                             modes = c("chull",
                                                       "boundaryPoints",
                                                       "boundaryPolygons"),
                                             cellFilterFactor = 10,
                                             minPixels = 20,
                                             cellPrefix = NULL,
                                             cellSuffix = NULL,
                                             plot = TRUE,
                                             verbose = TRUE) {
  
  # cellPrefix is the character to put before for naming convention
  # cellSuffix si the character to put after for naming convention
  # cellFilterFactor is the multiplier of the median pixel size to remove 
  # large segmented regions
  # modes is a character vector of the modes of boundary points or polygons 
  # to extract from this segmentation
  # minPixels is the minimum number of pixels to be included
  
  # cellPrefix = "Pos0_"
  # cellSuffix = "_4"
  
  require(reshape)
  require(rhdf5)
  if (plot) {
    require(ggplot2)
    require(patchwork)
  }
  require(grDevices) # chull package
  
  out = list()
  
  
  if (!all(modes %in% c("chull","boundaryPoints","boundaryPolygons"))) {
    stop("modes can only be the chull and one of boundaryPoints or boundaryPolygons")
  }
  if (!any(modes %in% c("chull","boundaryPoints","boundaryPolygons"))) {
    stop("Need at least one mode for extraction")
  }
  
  if (all(c("boundaryPoints", "boundaryPolygons") %in% modes)) {
    stop("Choose only one of boundaryPoints or boundaryPolygons to include in modes")
  }
  
  data_name = h5ls(filename)$name
  
  data = h5read(file = filename, data_name)
  
  # convert into a polygon
  allcells_raw = sort(unique(c(data)))
  
  # filter some 'cells' out based on too high pixel count
  # e.g. 10x the number of pixels compared to median
  cell_pixels = unclass(table(c(data)))
  allcells = as.numeric(names(cell_pixels)[cell_pixels <= cellFilterFactor*median(cell_pixels, na.rm = TRUE)])
  
  data_cens = data
  data_cens[! data_cens %in% allcells] <- NA
  
  data_long = reshape::melt(data_cens)
  data_long$valueFac = factor(as.character(data_long$value),
                              levels = sample(unique(as.character(data_long$value))))
  
  data_long_split = split.data.frame(data_long, data_long$valueFac)
  
  if ("chull" %in% modes) {
    chull_points_split = lapply(data_long_split, function(df){
      df_chull = df[chull(df$X1, df$X2),]
      df_chull = df_chull[chull(df_chull$X1, df_chull$X2),]
      return(df_chull)
    })
    
    chull_points = do.call(rbind, chull_points_split)
    
    if (plot) {
      g = ggplot(chull_points, aes(x = X1, y = -X2)) + 
        geom_point(aes(colour = valueFac), show.legend = FALSE) + 
        geom_polygon(aes(fill = valueFac), alpha = 0.2, show.legend = FALSE) +
        theme_classic() +
        ggtitle("Convex hull") +
        NULL
      g
    }
    
    chull_points$cell = paste0(cellPrefix, as.character(chull_points$valueFac), cellSuffix)
    
    out[["chull"]] <- chull_points
    
    if (plot) {
      g = ggplot(chull_points, aes(x = X1, y = -X2)) + 
        geom_point(aes(colour = valueFac), show.legend = FALSE) + 
        geom_polygon(aes(fill = valueFac), alpha = 0.7, show.legend = FALSE) +
        theme_classic() +
        ggtitle("Convex hull") +
        NULL
      
      g + geom_text(aes(x = X1, y = -X2, label = cell), 
                    data = data.frame(
                      X1 = tapply(chull_points$X1,chull_points$cell, mean),
                      X2 = tapply(chull_points$X2,chull_points$cell, mean),
                      cell = names(tapply(chull_points$X1,chull_points$cell, mean))),
                    size = 3)
    }
  }
  
  if (any(c("boundaryPoints", "boundaryPolygons") %in% modes)) {
    
    mode = grep("boundaryPoints|boundaryPolygons", modes, value = TRUE)[1]
    
    xmax = max(data_long$X1)
    ymax = max(data_long$X2)
    
    boundary_points_all_split = sapply(allcells, function(cell) {
      
      if (verbose) {
        print(paste0("cell ", cell))
      }
      
      
      # create a sparse array with only that cell
      data_cell = 1*(data[,,1] == cell)
      
      if (sum(data_cell) < minPixels) {
        if (verbose) print("Not enough pixels!")
        return(NULL)
        
        
        # for all pixels with 1s, calculate if at least one neighbour is a zero
        nonzero_cells = which(data_cell == 1, arr = TRUE)
        colnames(nonzero_cells) <- c("X1","X2")
      } 
      
      
      
      if (nrow(nonzero_cells) == 0) return(NULL)
      
      diffneighbour = apply(as.matrix(nonzero_cells[,c("X1","X2")]), 1, function(coord) {
        neighbour_coords = t(coord + t(matrix(c(-1,0,1,-1,1,-1,0,1,
                                                1,1,1,0,0,-1,-1,-1), nrow = 8, ncol = 2)))
        neighbour_coords = neighbour_coords[
          neighbour_coords[,1] > 0 & neighbour_coords[,1] <= xmax & neighbour_coords[,2] > 0 & neighbour_coords[,2] <= ymax,                             ]
        
        
        
        
        # are any of the neighbours zero
        if (sum(diag(data_cell[neighbour_coords[,1], neighbour_coords[,2]]) == 0) >=2 ) {
          out = 1
        } else {
          out = 0
        }
        
        
        # if the pixel is on the edge of the range then it should automatically be included
        if (any(coord == 1) | coord[1] == xmax | coord[2] == ymax) {
          out = 1
        }
        
        return(out)
      })
      
      
      nonzero_cells_df = data.frame(
        X1 = nonzero_cells[,1],
        X2 = nonzero_cells[,2],
        neighbour = diffneighbour
      )
      
      
      if (plot) {
        ggplot(nonzero_cells_df, aes(x = X1, y = X2)) + 
          geom_point(aes(colour = factor(neighbour)), size = 10) + 
          theme_classic()
      }
      
      boundary_points_df = subset(nonzero_cells_df, neighbour == 1)
      
      if (!"boundaryPolygons" %in% modes) {
        boundary_points_df$valueFac <- factor(cell, levels = levels(data_long$valueFac))
        return(boundary_points_df)
      }
      
      
      
      # order cells based on closest
      dmat = as.matrix(dist(boundary_points_df[,c("X1","X2")], method = "manhattan"))
      boundary_points_df$closestOrdering = c(1, rep(NA, nrow(boundary_points_df)-1))
      nties = 0
      for (i in 2:nrow(boundary_points_df)) {
        dvec = dmat[which.max(boundary_points_df$closestOrdering),]
        dvec[!is.na(boundary_points_df$closestOrdering)] <- NA
        
        
        if (sum(dvec == min(dvec, na.rm = TRUE), na.rm = TRUE) > 1) {
          # need to break the tie using angles
          # maybe there's three way ties etc
          j = order(dvec)[1]
          if (verbose) {
            print("Tie detected!")
          }
          nties = nties + (sum(dvec == min(dvec, na.rm = TRUE), na.rm = TRUE)) - 1
          
        } else {
          j = order(dvec)[1]
        }
        boundary_points_df[j,"closestOrdering"] <- i
      }
      
      # pruning step
      # if the distance between last and first point is more
      # than last and second last, then remove the last point
      # and repeat
      
      if (nties > 0) {
        for (i in 1:max(1,nties+1)) {
          last = which.max(boundary_points_df$closestOrdering)
          first = which.min(boundary_points_df$closestOrdering)
          nthlast = order(boundary_points_df$closestOrdering)[nrow(boundary_points_df) - nties]
          
          if (dmat[first, nthlast] < dmat[first,last]) {
            if (verbose) {
              print("pruning a point!")
            }
            dmat <- dmat[!boundary_points_df$closestOrdering == last,
                         !boundary_points_df$closestOrdering == last]
            boundary_points_df <- boundary_points_df[!boundary_points_df$closestOrdering == max(boundary_points_df$closestOrdering),]
          }
        }
      }
      
      
      boundary_points_df_sorted = reshape::sort_df(boundary_points_df, "closestOrdering")
      
      
      boundary_points_df_sorted$valueFac <- factor(cell, levels = levels(data_long$valueFac))
      
      
      return(boundary_points_df_sorted[,!colnames(boundary_points_df_sorted) %in% c("neighbour","closestOrdering")])
    }, simplify = FALSE)
    boundary_points_all <- do.call(rbind, boundary_points_all_split)
    
    # chull_points$cell = paste0("Pos0_", as.character(chull_points$valueFac), "_4")
    boundary_points_all$cell = paste0(cellPrefix, as.character(boundary_points_all$valueFac), cellSuffix)
    
    out[[mode]] <- boundary_points_all
    
    if (plot) {
      g2 = ggplot(boundary_points_all, aes(x = X1, y = -X2)) + 
        geom_point(aes(colour = valueFac), show.legend = FALSE) + 
        geom_polygon(aes(fill = valueFac), alpha = 0.7, show.legend = FALSE) +
        theme_classic() +
        ggtitle("Boundary points") +
        NULL
      g2
      
      g2 + geom_text(aes(x = X1, y = -X2, label = cell), 
                     data = data.frame(
                       X1 = tapply(boundary_points_all$X1,boundary_points_all$cell, mean),
                       X2 = tapply(boundary_points_all$X2,boundary_points_all$cell, mean),
                       cell = names(tapply(boundary_points_all$X1,boundary_points_all$cell, mean))),
                     size = 3)
    }
    
    if (plot & "chull" %in% modes) {
      g = ggplot(chull_points, aes(x = X1, y = -X2)) + 
        geom_point(aes(colour = valueFac), show.legend = FALSE) + 
        geom_polygon(aes(fill = valueFac), alpha = 0.7, show.legend = FALSE) +
        theme_classic() +
        ggtitle("Convex hull") +
        NULL
      g
    }
    
  }
  
  if (length(out) == 1) {
    return(out[[1]])
  }
  return(out)
}
