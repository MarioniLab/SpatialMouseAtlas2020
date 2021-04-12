
processSeqFISH2 = function(params) {
  
  # port parameters into global environment
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
    
    pos = grep(posName, unlist(strsplit(mRNA_file, nameSep)), value = TRUE)[1]
    fov = paste0("fov.", as.numeric(gsub(posName,"",pos))+1)
    embryo = embryo_label[fov]
    Estage = Estage_label[fov]
    
    cell_files_sub = grep(paste0(pos, nameSep), cell_files, value = TRUE)
    cell_file = grep(paste0(pos, nameSep), cell_files, value = TRUE)[1]
    
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
    
    
    if (!"z" %in% colnames(cell_info)) {
      # if z value is provided in the cell_info file
      # if it's not needs to be extracted from filename
      
      # first make sure the correct file is pulled out
      z_val = substring(grep("z[0-9]", unlist(strsplit(mRNA_file, nameSep)), value = TRUE)[1],2,2)
      
      # reload the correct file
      if (!grepl(paste0("z",z_val), cell_file)) {
        # then load the correct one
        cell_file <- grep(paste0("z", z_val), cell_files_sub, value = TRUE)[1]
        cell_info = read.csv(cell_file, header = TRUE)
        head(cell_info)
      }
      
    }
    
    cell_info$z <- z_val
    
    # change mRNA_file z-values
    mRNA_df_raw$z <- z_val
    
    # give a unique name
    uniqueID = paste0(
      embryo, "_",
      pos, "_",
      "cell", cell_info[,cellIDName], "_",
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
    # if it's not, likely corresponds to mRNAs not identified to be
    # within a given segmented cell, i.e. cell0
    print(setdiff(uniqueID_mRNA, uniqueID))
    
    mRNA_df_raw$uniqueID <- uniqueID_mRNA
    
    fov_info_sub = fov_info[fov,,drop = FALSE]
    
    bounds_x = as.numeric(fov_info_sub[,paste0("bound_x_",1:4)])    
    bounds_y = as.numeric(fov_info_sub[,paste0("bound_y_",1:4)])
    
    bound_x = range(bounds_x)
    bound_y = range(bounds_y)
    
    subcoords = mRNA_df_raw[,c("x","y")]
    
    # setting the range as what's observed assumes that we identify cells 
    # at the edge
    # but can be the best estimate if we have no other info
    
    range_x = c(0,2048)
    range_y = c(0,2048)
    
    newcoords_x = (subcoords[,"x"] - range_x[1])*(bound_x[2] - bound_x[1])/range_x[2] + bound_x[1]
    newcoords_y = (subcoords[,"y"] - range_y[1])*(bound_y[2] - bound_y[1])/range_y[2] + bound_y[1]
    
    
    cell_info$x_global = (cell_info[, cell_infoX] - range_x[1])*(bound_x[2] - bound_x[1])/range_x[2] + bound_x[1]
    cell_info$y_global = (cell_info[, cell_infoY] - range_y[1])*(bound_y[2] - bound_y[1])/range_y[2] + bound_y[1]
    
    mRNA_df_raw$x_global = newcoords_x
    mRNA_df_raw$y_global = newcoords_y
    
    mRNA_counts_arr = tapply(mRNA_df_raw$intensity, 
                             list(mRNA_df_raw$uniqueID,
                                  mRNA_df_raw$geneID),
                             function(x) length(x))
    mRNA_counts_arr[is.na(mRNA_counts_arr)] <- 0
    
    mRNA_counts <- mRNA_counts_arr
    
    # only keep rows with corresponding cell info, although this should all be TRUE
    # losing mRNAs that arents associated with a specific cell
    mRNA_counts = mRNA_counts[rownames(mRNA_counts) %in% rownames(cell_info),]
    
    coords = cell_info[rownames(mRNA_counts),]
    coords$pos = pos
    coords$fov = fov
    coords$embryo = embryo
    coords$Estage = Estage
    
    # grab boundary segmentation
    # copy segmentation files into main folder and grab their pattern
    segmentation_files = list.files(segmentation_files_dir,
                                    pattern = paste0(pos, "_dapi_", as.numeric(z_val)-1, "z_labels.tif"),
                                    full.names = TRUE)
    
    # only keep zvals that appear in our count data
    
    boundary_polygons_zvals = list()
    
    zvals = z_val
    
    i = 1
    print(i)
    
    boundary_polygons = extractPolygonsFromSegmentationTIFF(
      segmentation_files[i],
      modes = c("boundaryPolygons"),
      cellFilterFactor = 5,
      minPixels = 50,
      cellPrefix = paste0(embryo, "_", pos, "_cell"),
      cellSuffix = paste0("_z", as.numeric(z_val)),
      plot = FALSE,
      verbose = FALSE
    )
    colnames(boundary_polygons) <- c("x","y","cell","uniqueID")
    
    
    boundary_polygons <- pruneStraightLines(boundary_polygons, group = "uniqueID")
    
    boundary_polygons$z <- as.numeric(z_val)
    
    boundary_polygons$x_global = (boundary_polygons$x - range_x[1])*(bound_x[2] - bound_x[1])/(range_x[2]) + bound_x[1]
    boundary_polygons$y_global = (boundary_polygons$y - range_y[1])*(bound_y[2] - bound_y[1])/(range_y[2]) + bound_y[1]
    
    boundary_polygons_zvals[[i]] <- boundary_polygons
    
    boundary_polygons <- do.call(rbind, boundary_polygons_zvals)
    
    g = ggplot(boundary_polygons, aes(x = x, y = -y)) + 
      geom_polygon(aes(fill = uniqueID), show.legend = FALSE, colour = "black") + 
      theme_classic() +
      coord_fixed() +
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
      coord_fixed() +
      NULL
    print(g)
    
    mRNA_sceList[[mRNA_file]] <- mRNA_sce
    mRNA_dfList[[mRNA_file]] <- mRNA_df_raw
    
    write(mRNA_file, "mRNA_files_completed.txt", append = TRUE)
    
  }
  
  out = list(mRNA_sceList = mRNA_sceList,
             mRNA_dfList = mRNA_dfList)
  
  return(out)
}




# function that takes in the segmentation filename and outputs a list
# of dataframes, one for convex hull (chull) and another for boundary points
# with a flag for plotting

extractPolygonsFromSegmentationTIFF = function(filename,
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
  require(tiff)
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
  
  # from tiff package
  img = readTIFF(filename,
                 all = TRUE,
                 info = TRUE,
                 indexed = TRUE,
                 as.is = TRUE)
  data = t(as.array(img[[1]]))
  
  print(length(unique(c(data))))
  
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
  # take a given cell
  # cell = 50
  
  
  if (any(c("boundaryPoints", "boundaryPolygons") %in% modes)) {
    
    mode = grep("boundaryPoints|boundaryPolygons", modes, value = TRUE)[1]
    
    xmax = max(data_long$X1)
    ymax = max(data_long$X2)
    
    boundary_points_all_split = sapply(allcells, function(cell) {
      
      if (verbose) {
        print(paste0("cell ", cell))
      }
      
      # create a sparse array with only that cell
      data_cell = 1*(data == cell)
      
      if (sum(data_cell) < minPixels) {
        if (verbose) print("Not enough pixels!")
        return(NULL)
      }
      
      # for all pixels with 1s, calculate if at least one neighbour is a zero
      nonzero_cells = which(data_cell == 1, arr = TRUE)
      colnames(nonzero_cells) <- c("X1","X2")
      
      
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
        # print(i)
        dvec = dmat[which.max(boundary_points_df$closestOrdering),]
        # dvec[i] <- NA
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

