# # inpolygon function
# 
# sce = readRDS("/Users/ghazan01/Dropbox/Backup/SpatialEmbryos/analysis_output/E8.5/E8.5_sce_filt.Rds")
# 
# sce
# 
# source("/Users/ghazan01/Dropbox/Backup/SpatialEmbryos/scripts/spatial_functions.R")
# source("/Users/ghazan01/Dropbox/Backup/SpatialEmbryos/scripts/segmentation_functions.R")
# 
# polygon_df = getSegmentationVerticesDF(
#     colData(sce),
#     xname = "segmentation_vertices_x_global_affine",
#     yname = "segmentation_vertices_y_global_affine",
#     othercols = c("uniqueID","z","embryo","pos")
# )
# 
# polygon_df_single = subset(polygon_df, uniqueID == "embryo3_Pos39_cell99_z4")
# # polygon_df_single = subset(polygon_df, uniqueID == "embryo1_Pos0_cell10_z2")
# 
# # outside
# xval = 2.89
# yval = -0.5
# 
# # inside
# xval = 2.87
# yval = -0.48
# 
# g = ggplot(polygon_df_single, aes(x = segmentation_vertices_x_global_affine,
#                                   y = segmentation_vertices_y_global_affine)) + 
#     geom_polygon(alpha = 0.5) + 
#     geom_point(x = xval, y = yval)
# g
# 
# plot(polygon_df_single$segmentation_vertices_x_global_affine,
#      polygon_df_single$segmentation_vertices_y_global_affine, type = "n")
# text(polygon_df_single$segmentation_vertices_x_global_affine,
#      polygon_df_single$segmentation_vertices_y_global_affine,
#      labels = 1:length(polygon_df_single$segmentation_vertices_y_global_affine))
# 
# library(pracma)
# inpolygon(xval, yval, 
#           polygon_df_single$segmentation_vertices_x_global_affine,
#           polygon_df_single$segmentation_vertices_y_global_affine,
#           boundary = TRUE)
# 
# poly_length(polygon_df_single$segmentation_vertices_x_global_affine,
#             polygon_df_single$segmentation_vertices_y_global_affine)
# 
# abs(polyarea(polygon_df_single$segmentation_vertices_x_global_affine,
#              polygon_df_single$segmentation_vertices_y_global_affine))
# 
# # need to close the polygon for this to work
# center = poly_center(c(polygon_df_single$segmentation_vertices_x_global_affine,
#                        polygon_df_single$segmentation_vertices_x_global_affine[1]),
#                      c(polygon_df_single$segmentation_vertices_y_global_affine,
#                        polygon_df_single$segmentation_vertices_y_global_affine[1])) 
# center
# 
# # given a search polygon, find overlaps
# # search_x = c(2.87,2.87, 2.86, 2.86)
# # search_y = c(-0.55, -0.48, -0.48, -0.55)
# 
# # search_x = c(2.87,2.87, 2.86, 2.86)
# # search_y = c(-0.55, -0.525, -0.525, -0.55)
# 
# search_x = c(2.87,2.87, 2.86, 2.86)
# search_y = c(-0.475, -0.478, -0.478, -0.475)
# 
# search_df = data.frame(
#     search_x = search_x,
#     search_y = search_y
# )
# 
# g + geom_polygon(aes(x = search_x, y = search_y), data = search_df, alpha = 0.5, fill = "blue")
# 
# # check if inside polygon, if yes then overlaps
# inpolygon(search_x, search_y, 
#           polygon_df_single$segmentation_vertices_x_global_affine,
#           polygon_df_single$segmentation_vertices_y_global_affine,
#           boundary = TRUE)
# 
# # if any crossings then yes overlaps
# poly_crossings(rbind(polygon_df_single$segmentation_vertices_x_global_affine,
#                      polygon_df_single$segmentation_vertices_y_global_affine),
#                rbind(search_df$search_x,
#                      search_df$search_y))

# given the center of a polygon, expand or contract the polygon vertices 
# by a factor
expandPoint = function(x0,y0,x,y,c) {
    # x0, y0 is the centre x,y-coords
    # x, y is the point x,y-coords
    # c is the factor to multiply vector by
    
    dy = y - y0
    dx = x - x0
    tan_sign = tan(atan2(dy,dx))
    
    dc = ifelse(c >= 1, 1-c, c-1)
    csign = ifelse(c >= 1, 1, -1)
    
    if (dx == 0) {
        
        y_e = y0 + c*dy
        
        return(c(x, y_e))    
    }
    
    if (dy == 0) {
        
        x_e = x0 + c*dx
        
        return(c(x_e, y))    
    }
    
    Dx = 1 + (dy/dx)^2
    
    z2 = dx^2 + dy^2
    
    x_e = x + csign*sign(dx)*sqrt((dc^2)*z2/Dx)
    
    y_e = y + (x_e-x)*tan_sign
    
    return(c(x_e,y_e))
}

expandPointV = Vectorize(expandPoint)

closeVertices = function(x) {c(x,x[1])}

# pt_all = t(expandPointV(center[1],center[2],
#                         polygon_df_single$segmentation_vertices_x_global_affine,
#                         polygon_df_single$segmentation_vertices_y_global_affine,
#                         1))

# pt_all = NULL
# for (i in 1:length(polygon_df_single$segmentation_vertices_x_global_affine)) {
#     pt = expandPoint(center[1],center[2],
#                      polygon_df_single$segmentation_vertices_x_global_affine[i],
#                      polygon_df_single$segmentation_vertices_y_global_affine[i],
#                      2)
#     pt_all <- rbind(pt_all, pt)
# }

# g + geom_point(colour = "red") + 
#     geom_point(x = pt_all[,1], y = pt_all[,2],colour = "blue") +
#     geom_polygon(x = pt_all[,1], y = pt_all[,2],colour = "blue", fill = NA) +
#     xlim(c(2.8, 2.95)) + 
#     ylim(c(-0.55, -0.40))



# given a dataframe - with a column to split by - output two extra columns for
# expanded vertices
addExpandedVertices = function(verticesDF,
                               xname = "segmentation_vertices_x_global_affine",
                               yname = "segmentation_vertices_y_global_affine",
                               group = "uniqueID",
                               expansionFactor = 1.1,
                               new_xname = paste0(xname, "_expanded"),
                               new_yname = paste0(yname, "_expanded")) {
    
    require(pracma)
    
    verticesDFList = split.data.frame(verticesDF,
                                      verticesDF[,group])
    
    # get expanded vertices
    verticesDFList <- lapply(verticesDFList, function(verticesDF) {
        
        closedVertices_x = closeVertices(verticesDF[,xname])
        closedVertices_y = closeVertices(verticesDF[,yname])
        
        centre = poly_center(closedVertices_x,
                             closedVertices_y)
        
        newvertices = t(expandPointV(centre[1], centre[2],
                                     verticesDF[,xname],
                                     verticesDF[,yname],
                                     c = expansionFactor))
        
        verticesDF[,new_xname] <- newvertices[,1]
        verticesDF[,new_yname] <- newvertices[,2]
        
        return(verticesDF)
    })
    
    # plot
    verticesDF_long = do.call(rbind, verticesDFList)
    
    return(verticesDF_long)
}









# given a (named) list of vertex dataframes, output a graph of overlapping
# polygons according to an expansion factor

# polygon_df_sub = subset(polygon_df, embryo == "embryo1" & pos %in% c("Pos0", "Pos1") & z == 3)
# 
# verticesDFList = split.data.frame(polygon_df_sub,
#                                   polygon_df_sub$uniqueID)

neighbourVertices = function(verticesDFList,
                             xname = "segmentation_vertices_x_global_affine",
                             yname = "segmentation_vertices_y_global_affine",
                             expansionFactor = 1.1,
                             plot = FALSE,
                             plot2 = TRUE,
                             full = FALSE,
                             verbose = FALSE) {
    
    # full means check overlaps of polygons not just vertices,
    # this might be (much?) slower
    # note this should ideally be done within a single z-stack,
    # might be useful to do over multiple but will need to think how to do it
    
    require(pracma)
    require(igraph)
    require(gtools)
    
    # get expanded vertices
    if (expansionFactor != 1) {
    verticesDFList <- lapply(verticesDFList, function(verticesDF) {
        
        if (verbose) {
            print(verticesDF$uniqueID[1])
            }
        
        if (plot) {
            g = ggplot(verticesDF, aes(x = get(xname),
                                       y = get(yname))) + 
                geom_polygon(fill = NA, colour = "black")
            print(g)
        }
        
        closedVertices_x = closeVertices(verticesDF[,xname])
        closedVertices_y = closeVertices(verticesDF[,yname])
        
        centre = poly_center(closedVertices_x,
                             closedVertices_y)
        
        newvertices = t(expandPointV(centre[1], centre[2],
                                     verticesDF[,xname],
                                     verticesDF[,yname],
                                     c = expansionFactor))
        
        verticesDF[,"expanded_vertices_x"] <- newvertices[,1]
        verticesDF[,"expanded_vertices_y"] <- newvertices[,2]
        
        if (plot) {
            
            g + geom_polygon(aes(x = expanded_vertices_x, y = expanded_vertices_y),
                             data = verticesDF,
                             inherit.aes = FALSE,
                             fill = NA, colour = "blue")
            
        }
        
        
        return(verticesDF)
    })
    } else {
        verticesDFList <- lapply(verticesDFList, function(verticesDF) {
            verticesDF[,"expanded_vertices_x"] <- verticesDF[,xname]
            verticesDF[,"expanded_vertices_y"] <- verticesDF[,yname]
            return(verticesDF)
        })
    }
    
    # plot
    verticesDF_long = do.call(rbind, verticesDFList)
    
    if (verbose) {
        print("Expanded polygons")
        write("Expanded polygons", file = "status.txt", append = TRUE)
    }
    
    if (plot2) {
        g = ggplot(verticesDF_long, aes(x = get(xname),
                                    y = get(yname))) + 
            geom_polygon(aes(group = uniqueID), fill = NA, colour = "black") + 
            geom_polygon(aes(x = expanded_vertices_x, y = expanded_vertices_y,
                             group = uniqueID),
                         inherit.aes = FALSE, fill = NA, colour = "blue") + 
            theme_classic() + 
            NULL
        print(g)
    }
    
    # build overlap graph
    
    # do this in a loop
    verticesDFList_nms <- mixedsort(names(verticesDFList))
    neighbourSegmentsList = as.list(verticesDFList_nms)
    names(neighbourSegmentsList) <- verticesDFList_nms
    for (segment in verticesDFList_nms) {
        
        if (verbose) {
            print(segment)
        }
        
        # ask if any cell segmentation vertices are within expanded region
        # this is fast-ish
        
        inPoly = inpolygon(verticesDF_long[,xname], 
                           verticesDF_long[,yname], 
                           verticesDFList[[segment]][,"expanded_vertices_x"],
                           verticesDFList[[segment]][,"expanded_vertices_y"],
                           boundary = TRUE)
        
        neighbourSegments = unique(verticesDF_long[inPoly,"uniqueID"])
        
        if (full) {
            # check using polygon crossings
            polycrossingList = lapply(verticesDFList, function(verticesDF) {
                polycrossing = poly_crossings(t(verticesDFList[[segment]][,c("expanded_vertices_x","expanded_vertices_y")]),
                                              t(verticesDF[,c(xname,yname)]))
                if (is.null(polycrossing)) return(FALSE)
                return(TRUE)
            })
            neighbourSegmentsFull = names(verticesDFList)[unlist(polycrossingList)]
            
            neighbourSegments <- union(neighbourSegments, neighbourSegmentsFull)
        }
        
        neighbourSegmentsList[[segment]] <- sort(setdiff(neighbourSegments,segment))
        
    }
    
    neighbourSegmentsPairs = cbind(rep(names(neighbourSegmentsList),
                                       times = unlist(lapply(neighbourSegmentsList,length))),
                                   unlist(neighbourSegmentsList))
    
    if (verbose) {
        print("Got neighbour segment pairs")
        write("Got neighbour segment pairs", file = "status.txt", append = TRUE)
    }
    
    if (plot2) {
        
        centres = do.call(rbind,lapply(verticesDFList, function(verticesDF) {
            poly_center(closeVertices(verticesDF[,xname]),
                        closeVertices(verticesDF[,yname]))
        }))
        
        centresDF = data.frame(
            segment = rownames(centres),
            segment_x = centres[,1],
            segment_y = centres[,2]
        )
        
        neighbourSegmentsDF = data.frame(
            segment1 = neighbourSegmentsPairs[,1],
            segment2 = neighbourSegmentsPairs[,2],
            segment1_x = centres[neighbourSegmentsPairs[,1],1],
            segment1_y = centres[neighbourSegmentsPairs[,1],2],
            segment2_x = centres[neighbourSegmentsPairs[,2],1],
            segment2_y = centres[neighbourSegmentsPairs[,2],2]
        )
        
        g = ggplot(verticesDF_long, aes(x = get(xname),
                                        y = get(yname))) +
            geom_polygon(aes(group = uniqueID), fill = NA, colour = "grey", size = 0.5) + 
            geom_point(aes(x = segment_x, y = segment_y), inherit.aes = FALSE,
                       data = centresDF) +
            geom_segment(aes(x = segment1_x, xend = segment2_x,
                             y = segment1_y, yend = segment2_y),
                         inherit.aes = FALSE,
                         data = neighbourSegmentsDF,
                         size = 0.5, colour = "blue") +
            theme_classic() + 
            NULL
        
        print(g)
        
    }
    
    neighbourSegmentsGraph = simplify(graph.edgelist(neighbourSegmentsPairs, 
                                                     directed = FALSE))
    
    if (verbose) {
        print("Got neighbour segment graph")
        write("Got neighbour segment graph", file = "status.txt", append = TRUE)
    }
    
    return(neighbourSegmentsGraph)
}


# out = neighbourVertices(verticesDFList,
#                         xname = "segmentation_vertices_x_global_affine",
#                         yname = "segmentation_vertices_y_global_affine",
#                         expansionFactor = 1.1,
#                         plot = FALSE,
#                         plot2 = TRUE,
#                         full = FALSE)

# given a set of 2D coordinates, calculate the network distance
# trim extra long edges using 1st quartile

get_delaunay = function(coords,
                        euc_filter = "quantile",
                        euc_filter_value = 0.25,
                        plot = TRUE) {
    
    # coords must be 2d, columns x- and y-
    # euc_filter either "quantile" or "absolute" or "none" 
    # output is a list containing the network and the network distances
    
    require(spatstat)
    require(igraph)
    
    if (euc_filter != "none") {
        euc_dist = as.matrix(dist(coords))
    }
    
    if (euc_filter == "quantile") {
        euc_dist_max = quantile(euc_dist[lower.tri(euc_dist, diag = FALSE)], euc_filter_value)
    }
    
    if (euc_filter == "absolute") {
        euc_dist_max = euc_filter_value
    }
    
    # euc_dist_max
    
    win_val = owin(xrange = range(coords[,1]), yrange = range(coords[,2]))
    
    # win_val = do.call(owin, unlist(apply(coords,2,function(x) list(range(x))), recursive = FALSE))
    
    # ppp only in 2D
    Y = ppp(coords[,1], coords[,2], window = win_val)
    
    out_net = delaunayNetwork(Y)
    
    full_edgelist = cbind(rownames(coords)[out_net[["from"]]],
                          rownames(coords)[out_net[["to"]]])
    
    out_edgelist = rbind(full_edgelist,
                         cbind(rownames(coords),
                               rownames(coords)))
    
    # prune the edgelist based on euclidean distance
    dist_edgelist = diag(euc_dist[out_edgelist[,1],out_edgelist[,2]])
    
    # table(dist_edgelist < euc_dist_max)
    
    if (euc_filter != "none") {
        out_edgelist <- out_edgelist[dist_edgelist < euc_dist_max,]
    }
    
    out_graph = igraph::simplify(igraph::graph.edgelist(out_edgelist, directed = FALSE))
    
    out_dist = igraph::distances(out_graph)[rownames(coords), rownames(coords)]
    
    if (plot) {
        par(mfrow=c(1,2))
        out_graph_plot = out_graph
        ly = coords[V(out_graph_plot)$name,]
        V(out_graph_plot)$size = 0.05
        V(out_graph_plot)$label <- ""
        
        # it gets weird along the boundaries
        # dev.off()
        plot(delaunayNetwork(Y), main = "delaunay triangulation")
        # plot(coords)
        cexval = 1
        i = 1
        points(coords[i,1], coords[i,2], pch = 16, cex = cexval)
        points(coords[out_dist[i,] == 1,1],coords[out_dist[i,] == 1,2], pch = 16, col = "red", cex = cexval)
        points(coords[out_dist[i,] == 2,1],coords[out_dist[i,] == 2,2], pch = 16, col = "orange", cex = cexval)
        points(coords[out_dist[i,] == 3,1],coords[out_dist[i,] == 3,2], pch = 16, col = "yellow", cex = cexval)
        points(coords[out_dist[i,] == 4,1],coords[out_dist[i,] == 4,2], pch = 16, col = "green", cex = cexval)
        points(coords[out_dist[i,] == 5,1],coords[out_dist[i,] == 5,2], pch = 16, col = "blue", cex = cexval)
        
        plot(out_graph_plot, layout = ly, main = "pruned graph")
        
    }
    
    # out_weight = 1 - out_dist/max(c(out_dist))
    
    return(list(dist = out_dist, graph = out_graph, delaunay = full_edgelist))
    
}






makeTriangular = function(mat) {
    # take a matrix assumed to be square and return triangular
    # based on the ordering given
    mat[is.na(mat)] <- 0
    
    matLower = mat
    matLower[upper.tri(matLower, diag = TRUE)] <- 0
    
    matUpper = mat
    matUpper[lower.tri(matUpper, diag = FALSE)] <- 0
    
    matSym = matLower + t(matUpper)
    matSym[upper.tri(matSym)] <- NA
    
    return(matSym)
}


splitChunk = function(N, N_chunk) {
    # N is a number
    # N_chunk is the number to split into
    # output is a list containing start and end indices
    
    # add an error if either is not positive integer
    
    if (N <= N_chunk) {
        return(list(c(1,N)))
    }
    
    # e.g.
    # N = 300
    # N_chunk = 100
    
    out = list()
    out[[1]] <- c(1)
    
    while (rev(unlist(out))[1] < N) {
        
        n = length(out)
        out[[n]][2] <- out[[n]][1] + (N_chunk - 1)
        if (out[[n]][2] >= N) {
            out[[n]][2] <- N
        } else {
            out[[n+1]] <- out[[n]][2] + 1
        }
        # print(out)
        
    }
    
    return(out)
}

getRandomConnectivity = function(sce, 
                                 neighbours,
                                 group,
                                 option = "observed", 
                                 x_name = "x_global_affine",
                                 y_name = "y_global_affine",
                                 splitgroups = c("embryo", "z"),
                                 chunksize = 10000) {
    ######## new function with more functionality
    # option either observed or random
    # sce singlecellexperiment
    # neighbours a two column matrix with cell names
    # here group is a matrix of prior probabilities of cell types
    
    if (length(group) == 1 & class(group)[1] == "character") {
        # convert to matrix with 1s and 0s
        
        fac = colData(sce)[,group]
        levels = as.character(sort(unique(fac)))
        
        group = do.call(cbind,sapply(levels, function(ctype){
            1*(fac == ctype)
        },simplify = FALSE))
        rownames(group) <- colnames(sce)
    }
    
    group <- group[colnames(sce),]
    
    if (option != "observed") {
        # randomly sample but within the splitgroups
        group_split = split.data.frame(group,
                                       interaction(as.data.frame(colData(sce)[,splitgroups])))
        group_split_randomised = lapply(group_split, function(subgroup) {
            subgroup_randomised <- subgroup[sample(rownames(subgroup)),]
            rownames(subgroup_randomised) <- rownames(subgroup)
            return(subgroup_randomised)
        })
        group_randomised = do.call(rbind, group_split_randomised)[rownames(group),]
        group <- group_randomised
    }
    
    
    nb_arr = 0
    chunks = splitChunk(nrow(neighbours), chunksize)
    
    for (nb in 1:length(chunks)) {
        # print(nb)
        nb_arr_i = t(group[neighbours[chunks[[nb]][1]:chunks[[nb]][2],1],]) %*% group[neighbours[chunks[[nb]][1]:chunks[[nb]][2],2],]
        nb_arr <- nb_arr + nb_arr_i
    }
    
    neighbourSimilaritySym = makeTriangular(nb_arr)
    
    return(neighbourSimilaritySym)
}


getRandomConnectivity_dep = function(sce, 
                                 neighbours,
                                 group,
                                 option = "observed", 
                                 x_name = "x_global_affine",
                                 y_name = "y_global_affine",
                                 splitgroups = c("embryo", "z")) {
    ##### deprecated function!!!
    # option either observed or random
    # sce singlecellexperiment
    # neighbours a two column matrix with cell names
    
    sce$factor_random = unsplit(tapply(colData(sce)[,group],
                                       interaction(as.data.frame(colData(sce)[,splitgroups])),
                                       sample),
                                interaction(as.data.frame(colData(sce)[,splitgroups])))
    
    neighbours_df = data.frame(
        x1 = colData(sce)[neighbours[,1], x_name],
        x2 = colData(sce)[neighbours[,2], x_name],
        y1 = colData(sce)[neighbours[,1], y_name],
        y2 = colData(sce)[neighbours[,2], y_name],
        ct1 = colData(sce)[neighbours[,1], 
                           ifelse(option == "observed", group, "factor_random")],
        ct2 = colData(sce)[neighbours[,2], 
                           ifelse(option == "observed", group, "factor_random")]
    )
    for (i in 1:length(splitgroups)) {
        neighbours_df[, splitgroups[i]] <- colData(sce)[neighbours[,1],splitgroups[i]]
    }
    
    neighbours_df$ctSame = ifelse(
        as.character(neighbours_df$ct1) == as.character(neighbours_df$ct2),
        "Same","Different")
    
    neighbourSimilarity = unclass(
        table(factor(neighbours_df$ct1, levels = sort(unique(colData(sce)[,group]))),
              factor(neighbours_df$ct2, levels = sort(unique(colData(sce)[,group])))))
    
    neighbourSimilaritySym = makeTriangular(neighbourSimilarity)
    
    return(neighbourSimilaritySym)
}


cellCellContact = function(sce,
                           group,
                           graph,
                           nperm = 500,
                           plot = TRUE,
                           cellID = "uniqueID",
                           ...) {
    require(igraph)
    require(abind)
    
    graph_sub = induced_subgraph(graph,
                                 which(V(graph)$name %in% as.character(colData(sce)[,cellID])))
    
    graph_pairs = get.edgelist(graph_sub)
    
    mat_obs = getRandomConnectivity(sce,
                                    neighbours = graph_pairs,
                                    group = group,
                                    option = "observed", 
                                    ...)
    if (FALSE) { # plot the observed number neighbours
        Heatmap(log10(1+mat_obs[,rev(colnames(mat_obs))]), cluster_rows = FALSE, 
                cluster_columns = FALSE,
                col = c("white","gray75","cornflowerblue","black","black"),
                name = "Number\nneighbours\n(log10)",
                na_col = "white",
        )
    }
    
    mat_randomList = replicate(nperm,
                               getRandomConnectivity(sce,
                                                     neighbours = graph_pairs,
                                                     group = group,
                                                     option = "random",
                                                     ...))
    
    mat_random_both = mat_randomList
    # mat_random_both[1,1,n_perm+1] <- mat_obs
    mat_random_both = abind::abind(mat_random_both, mat_obs, along = 3)
    
    mat_p = apply(mat_random_both,c(1,2),
                  function(x) mean(x[nperm+1] >= x[1:nperm]))
    mat_p_sym = mat_p
    mat_p_sym[upper.tri(mat_p_sym)] <- t(mat_p_sym)[upper.tri(mat_p_sym)]
    
    mat_obs_sym = mat_obs
    mat_obs_sym[upper.tri(mat_obs_sym)] <- t(mat_obs_sym)[upper.tri(mat_obs_sym)]
    
    out = list(obs = mat_obs_sym, pmat = mat_p_sym)
    
    if (plot) {
        g = cellCellContactMap(out)
        print(g)
    }
    
    return(out)
}

cellCellContactMap = function(out, order = NULL, exclude = NULL) {
    # this function takes output from cellCellContact() 
    # and gives a ggplot object of the graph
    # out = mat_p_sym
    # if order not given then perform hclust to get ordering
    
    require(reshape)
    require(ggplot2)
    
    obs = out[["obs"]]
    pmat = out[["pmat"]]
    
    mat_p_df = melt(pmat)
    colnames(mat_p_df) <- c("subcluster_1","subcluster_2", "pvalue")
    mat_p_df$Sig <- ifelse(abs(mat_p_df$pvalue - 0.5) > (0.5 - 0.01), "*", "")
    mat_p_df$Obs <- c(obs)
    mat_p_df <- na.omit(mat_p_df)
    
    if (is.null(order)) {
    hc = hclust(dist(pmat))
    mat_p_df$subcluster_1 <- factor(mat_p_df$subcluster_1, levels = 
                                        hc$labels[hc$order])
    mat_p_df$subcluster_2 <- factor(mat_p_df$subcluster_2, levels = 
                                        hc$labels[hc$order])
    } else {
        mat_p_df$subcluster_1 <- factor(mat_p_df$subcluster_1, levels = 
                                            order)
        mat_p_df$subcluster_2 <- factor(mat_p_df$subcluster_2, levels = 
                                            order)
    }
    
    mat_p_df$keep = as.numeric(mat_p_df$subcluster_1) >= as.numeric(mat_p_df$subcluster_2) 
    
    g = ggplot(subset(mat_p_df, keep & 
                          (!subcluster_1 %in% exclude) & 
                          (!subcluster_2 %in% exclude)),
               aes(x = subcluster_1, y = subcluster_2, fill = pvalue)) + 
        geom_tile() +
        # geom_text(aes(label = Sig, y = subcluster_2 - 0.025), size = 20) +
        # geom_point(pch = "*", size = 10, data = subset(mat_p_df, Sig == "*" & keep)) +
        theme_classic() +
        theme(axis.line = element_blank()) +
        theme(axis.ticks = element_blank()) +
        # scale_y_continuous(position = "right") +
        scale_y_discrete(position = "right") +
        theme(axis.text = element_text(size = 10)) +
        theme(axis.title = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        coord_fixed() +
        # guides(fill = guide_legend(title = "")) +
        scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "#ed6495", na.value = "white",
                             midpoint = 0.5, limits = c(0,1),
                             breaks = c(0,1), labels = c("Segregated", "Integrated")) +
        theme(legend.position = "top") + 
        theme(legend.text = element_text(size = 15)) +
        theme(legend.key.width = unit(1, "in")) +
        guides(fill = guide_colourbar(title.position = "top",
                                      title = "",
                                      title.hjust = 0.5,
                                      ticks = FALSE,
                                      reverse = TRUE)) +
        NULL
    # print(g)
    
    return(g)
    
}


cellCellContactHeatmap = function(out,
                                  col,
                                  split_n = NULL,
                                  order = NULL,
                                  exclude = NULL,
                                  factor = 1) {
    # out is the output from cellcellcontactmap
    # col is a named character vector of the names of the columns
    
    require(ComplexHeatmap)
    
    pmat = out[["pmat"]]
    
    if (!is.null(exclude)) {
        pmat <- pmat[setdiff(rownames(pmat), exclude), setdiff(colnames(pmat), exclude)]
    }
    
    if (!is.null(order)) {
        pmat <- pmat[order,order]
        cluster = FALSE
    } else {
        cluster = TRUE
    }
    
    h = Heatmap(pmat,
                
                col = c("cornflowerblue","white","#ed6495"),
                
                top_annotation = columnAnnotation(foo = colnames(pmat),
                                                  col = list(foo = col),
                                                  show_annotation_name = FALSE,
                                                  show_legend = FALSE),
                
                right_annotation = rowAnnotation(foo = colnames(pmat),
                                                 col = list(foo = col),
                                                 show_annotation_name = FALSE,
                                                 show_legend = FALSE),
                
                show_heatmap_legend = FALSE,
                
                width = factor*unit(ncol(pmat), "cm"),
                height = factor*unit(ncol(pmat), "cm"),
                
                rect_gp = gpar(col = "grey", lwd = 0.5),
                
                row_split = split_n,
                column_split = split_n,
                
                column_title = NULL,
                row_title = NULL,
                column_names_rot = 45,
                row_names_rot = 0,
                row_names_side = "left",
                row_dend_side = "right",
                
                cluster_rows = cluster,
                cluster_columns = cluster,
                
                column_dend_gp = gpar(lwd = 2),
                row_dend_gp = gpar(lwd = 2)
                
    )
    
    return(h)
}


cellCellContactHeatmapTriangle <- function(out,
                                           col_ann,
                                           split_n = NULL,
                                           order = NULL,
                                           exclude = NULL,
                                           factor = 1, ...) {
    
    require(ComplexHeatmap)
    
    pmat = out[["pmat"]]
    
    if (!is.null(exclude)) {
        pmat <- pmat[setdiff(rownames(pmat), exclude), setdiff(colnames(pmat), exclude)]
    }
    
    if (!is.null(order)) {
        pmat <- pmat[order,order]
        cluster = FALSE
    } else {
        cluster = TRUE
    }
    
    mat2 = pmat
    
    if (cluster) {
        hc = hclust(dist(mat2), ...)
        od = hc$order
        # if (is.null(split_n)) {
        #     split_n = dynamicTreeCut::cutreeHybrid(hc)
        # }
        # mat2 = mat2[od, od]
    } else {
        hc = FALSE
        od = 1:nrow(mat2)
        split_n <- NULL
    }
    
    # col_ann = brain_cols
    # col_ann = celltype_colours
    
    # Now we can self-define the heatmap. In following code, col controls how to map values to colors. Since clustering is already applied, we set cluster_rows and cluster_columns to FALSE. When type is set to none in gpar(), the heatmap is initialized but nothing is added, then we can use the self-defined function cell_fun to add graphic into it.
    # 
    # cell_fun will be applied to each small grid in the heatmap. In cell_fun, there are seven arguments:
    #     
    #     j: column index in mat2
    # i: row index in mat2
    # x: position on the plot
    # y: position on the plot
    # w: width of the current grid
    # h: height of the current grid
    # col: filled color for the current grid
    # Use Heatmap() to generate the final heatmap:
    # library(ComplexHeatmap)
    
    # factor = 0.5
    # split_n = 11
    
    mat2_mask = mat2
    for (i in 1:nrow(mat2_mask)) {
        for (j in 1:nrow(mat2_mask)) {
            if (which(od == i) < which(od == j)) {
                mat2_mask[i,j] <- NA 
            }
        }
    }
    
    # hc_rev = hc
    # hc_rev$labels <- rev(hc_rev$labels)
    
    h = Heatmap(mat2_mask,
                col = c("cornflowerblue", "white", "#ed6495"),
                
                # cluster_rows = FALSE,
                # cluster_columns = FALSE,
                
                cluster_rows = hc,
                cluster_columns = hc,
                
                column_dend_side = "bottom",
                show_row_names = FALSE,
                show_column_names = FALSE,
                
                border = FALSE,
                
                show_heatmap_legend = FALSE,
                
                width = factor*unit(ncol(mat2), "cm"),
                height = factor*unit(ncol(mat2), "cm"),
                
                column_title = NULL,
                row_title = NULL,
                
                na_col = NA,
                
                bottom_annotation = columnAnnotation(foo = colnames(mat2),
                                                     col = list(foo = col_ann),
                                                     show_annotation_name = FALSE,
                                                     show_legend = FALSE),
                
                left_annotation = rowAnnotation(foo = colnames(mat2),
                                                col = list(foo = col_ann),
                                                show_annotation_name = FALSE,
                                                show_legend = FALSE),
                
                column_split = split_n,
                row_split = split_n,
                
                column_dend_gp = gpar(lwd = 2),
                row_dend_gp = gpar(lwd = 2),
                
                
                cell_fun = function(j, i, x, y, w, h, col) {
                    if (j == i) {
                        grid.rect(x, y, w, h, gp = gpar(fill = col, col = "grey"))
                        grid.text(rownames(mat2)[i], x + factor*unit(0.5, "cm"), y + factor*unit(0.5, "cm"), rot = 45, hjust = 0)
                    } else if (which(od == i) < which(od == j)) {
                        # grid.rect(x, y, w, h, gp = gpar(fill = NULL, col = NA))
                        grid.rect(x, y, w, h, gp = gpar(fill = col, col = NA))
                    } else {
                        grid.rect(x, y, w, h, gp = gpar(fill = col, col = "grey"))
                    }
                })
    
    # draw(h, padding = unit(c(15, 2, 2, 50), "mm"))
    return(h)

# example usage
# cellCellContactHeatmapTriangle(out, col_ann = celltype_colours,
#                                exclude = "Low quality",
#                                factor = 0.5,
#                                split_n = 10,
#                                method = "complete")
}
