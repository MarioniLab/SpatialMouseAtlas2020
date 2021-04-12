# For R versions lower thatn R4, use the following commit version for scHOT:
# devtools::install_github("shazanfar/scHOT", ref = "ec16c1f")

plotHigherOrderSequence = function(scHOT, 
                                   gene, 
                                   positionType = NULL, 
                                   branches = NULL,
                                   positionColData = NULL,
                                   voronoi_max.radius = 1) {
    if (!("higherOrderSequence" %in% colnames(scHOT@scHOT_output))) {
        stop("higherOrderSequence is not found in scHOT_output")
    }
    namescols = grep("gene_", colnames(scHOT@scHOT_output), value = TRUE)
    wcor <- as.matrix(scHOT@scHOT_output$higherOrderSequence)
    rownames(wcor) <- apply(scHOT@scHOT_output[, namescols, drop = FALSE], 
                            1, paste0, collapse = "_")
    colnames(wcor) <- NULL
    if (ncol(wcor) != ncol(scHOT)) {
        warning("Not all the cell position has higherOrderSequence statistics,\n            set nrow.out = NULL in scHOT_setWeightMatrix to calculate\n            higherOrderSequence for all positions!")
    }
    if (is.null(positionType)) {
        if (is.null(scHOT@positionType)) {
            stop("Both positionType and scHOT@positionType are NULL.")
        }
        else {
            positionType <- scHOT@positionType
        }
    }
    positionType <- match.arg(positionType, c("trajectory", "spatial"), 
                              several.ok = FALSE)
    if (positionType == "trajectory") {
        if (is.null(branches)) {
            message("branches information is not provided")
            wcorsList = list(Branch = wcor)
        }
        else {
            if (!(branches %in% colnames(colData(scHOT)))) {
                stop("branches provided is not found in colData(scHOT)")
            }
            branch_info <- SummarizedExperiment::colData(scHOT)[, 
                                                                branches]
            wcorsList <- lapply(split(seq_along(branch_info), 
                                      branch_info), function(idx) wcor[, idx])[order(unique(branch_info))]
        }
        if (is.null(names(wcorsList))) {
            names(wcorsList) <- paste0("Branch_", seq_len(length(wcorsList)))
        }
        if (is.logical(gene[1])) {
            if (length(unique(unlist(lapply(wcorsList, nrow)))) > 
                1) {
                stop("cannot use logical subset when weighted higher\n             order statistic matrices have differing rows")
            }
            if (length(gene) != nrow(wcorsList[[1]])) {
                stop("cannot use logical subset when length of gene\n             doesn't match nrow of wcorsList matrices")
            }
            wcors_longList = lapply(wcorsList, function(branch) {
                reshape::melt(t(branch[gene, ]))
            })
            gene = ""
        }
        else {
            gene = paste0(sort(gene), collapse = "|")
            wcors_longList = lapply(wcorsList, function(branch) {
                reshape::melt(t(branch[grepl(gene, rownames(branch)), 
                                       , drop = FALSE]))
            })
        }
        branch_long = do.call(rbind, wcors_longList)
        branch_long = cbind(rep(names(wcors_longList), unlist(lapply(wcors_longList, 
                                                                     nrow))), branch_long)
        colnames(branch_long) = c("branch", "SampleOrder", "GenePair", 
                                  "WeightedCorrelation")
        if (max(abs(branch_long$WeightedCorrelation)) < 1) {
            ylimit <- ylim(c(-1, 1))
        }
        else {
            ylimit <- NULL
        }
        g <- ggplot(branch_long, aes(x = branch_long$SampleOrder, 
                                     y = branch_long$WeightedCorrelation, group = branch_long$GenePair, 
                                     col = branch_long$GenePair)) + 
            geom_line(size = 2, 
                      alpha = 0.6) + facet_grid(~branch, scales = "free_x") + 
            theme_classic() + 
            ylimit + 
            geom_hline(yintercept = 0, 
                       size = 1, colour = "grey") + 
            ggtitle(gene) + 
            xlab("Sample Order") + 
            ylab("Weighted Higher Order Statistic") + 
            labs(col = "Test") + 
            NULL
    }
    if (positionType == "spatial") {
        if (is.logical(gene[1])) {
            if (length(unique(unlist(lapply(wcorsList, nrow)))) > 
                1) {
                stop("cannot use logical subset when weighted higher\n             order statistic matrices have differing rows")
            }
            if (length(gene) != nrow(wcor)) {
                stop("cannot use logical subset when length of gene\n             doesn't match nrow of wcorsList matrices")
            }
            gene = ""
        }
        else {
            if (!all(gene %in% rownames(wcor))) {
                if (length(gene) == 2) {
                    if (!paste0(sort(gene), collapse = "_") %in% 
                        rownames(wcor)) {
                        stop("gene pairs has no higherOrderSequence ")
                    }
                    else {
                        gene <- paste0(sort(gene), collapse = "_")
                    }
                }
                else {
                    stop("gene pairs has no higherOrderSequence ")
                }
            }
        }
        if (is.null(positionColData)) {
            if (is.null(scHOT@positionColData)) {
                stop("Both positionColData and scHOT@positionColData are NULL.")
            }
            else {
                positionColData <- scHOT@positionColData
            }
        }
        coords_info <- data.frame(colData(scHOT)[, positionColData])
        colnames(coords_info) <- positionColData
        wcor_all <- matrix(NA, nrow = length(gene), ncol = ncol(scHOT))
        rownames(wcor_all) <- gene
        colnames(wcor_all) <- seq_len(ncol(scHOT))
        wcor_all[gene, rownames(scHOT@weightMatrix)] <- wcor[gene, 
                                                             , drop = FALSE]
        branch_long <- reshape::melt(cbind(coords_info, t(wcor_all)), 
                                     id.vars = positionColData)
        colnames(branch_long) <- c("x", "y", "genepair", "value")
        g <- ggplot(branch_long, aes(x = branch_long$x, y = branch_long$y, 
                                     fill = branch_long$value)) + 
            ggforce::geom_voronoi_tile(max.radius = voronoi_max.radius) + 
            # geom_point(size = 0.5, colour = "black") + 
            theme_classic() + 
            facet_wrap(~genepair) + 
            scale_alpha_continuous(range = c(0,0.5)) + 
            scale_fill_gradient2(low = "blue", mid = "white",
                                 high = "red", limits = c(-1, 1)) + 
            theme(panel.grid = element_blank()) + 
            theme(axis.ticks = element_blank()) + 
            theme(axis.text = element_blank()) + 
            xlab("") + 
            ylab("") + 
            coord_fixed() + 
            labs(fill = "Weighted Higher Order Statistic") + 
            NULL
    }
    return(g)
}

plotOrderedExpression = function(scHOT, genes, positionType = NULL, branches = NULL, 
                                 ranked_by = NULL, xvals = NULL, subsetBranch = NULL, facet = FALSE, 
                                 positionColData = NULL, assayName = NULL) {
    if (is.null(positionType)) {
        if (is.null(scHOT@positionType)) {
            stop("Both positionType and scHOT@positionType are NULL.")
        }
        else {
            positionType <- scHOT@positionType
        }
    }
    if (is.null(assayName)) {
        assayName <- "expression"
    }
    branchData <- SummarizedExperiment::assay(scHOT, assayName)
    genes = genes[genes %in% rownames(branchData)]
    if (length(genes) == 0) {
        stop("No genes found in rownames of dataset!")
    }
    if (scHOT@positionType == "trajectory") {
        if (is.null(ranked_by)) {
            message("ranked_by information is not provided,\n              the expression data is ranked by the branches")
        }
        else {
            if (!(ranked_by %in% colnames(SummarizedExperiment::colData(scHOT)))) {
                stop("ranked_by provided is not found in colData(scHOT)")
            }
            branchData <- branchData[, order(colData(scHOT)[, 
                                                            ranked_by])]
        }
        if (is.null(branches)) {
            message("branches information is not provided")
            branchData = list(Branch = branchData)
        }
        else {
            if (!(branches %in% colnames(colData(scHOT)))) {
                stop("branches provided is not found in colData(scHOT)")
            }
            branch_info <- SummarizedExperiment::colData(scHOT)[, 
                                                                branches]
            branchData <- lapply(split(seq_along(branch_info), 
                                       branch_info), function(idx) branchData[, idx])[order(unique(branch_info))]
        }
        gdf_list = sapply(genes, function(g) {
            gdf = do.call(rbind, lapply(branchData, function(branch) {
                gdf_list_1 = data.frame(Sample = colnames(branch), 
                                        order = seq_len(ncol(branch)), ExpressionGene = branch[g, 
                                        ], gene = g)
                return(gdf_list_1)
            }))
            gdf$branch = rep(names(branchData), times = unlist(lapply(branchData, 
                                                                      ncol)))
            return(gdf)
        }, simplify = FALSE)
        gdf = do.call(rbind, gdf_list)
        if (!is.null(subsetBranch)) {
            gdf_sub = subset(gdf, gdf$branch %in% subsetBranch)
            if (nrow(gdf_sub) == 0) 
                stop("no branches with names in subsetBranch,\n                                 please re-run with correct names\n                                   (should match names of branchData)")
        }
        else {
            gdf_sub = gdf
        }
        if (!is.null(xvals)) {
            if (!(xval %in% colnames(SummarizedExperiment::colData(scHOT)))) {
                stop("xval provided is not found in colData(scHOT)")
            }
            xval <- lapply(split(seq_along(branch_info), branch_info), 
                           function(idx) xval[, idx])[order(unique(branch_info))]
            xval <- apply(as.matrix(gdf_sub), 1, function(x) xvals[[x["branch"]]][as.numeric(x["order"])])
            gdf_sub$order <- xval
        }
        g = ggplot(gdf_sub, aes(x = order, y = gdf_sub$ExpressionGene, 
                                colour = gdf_sub$gene, fill = gdf_sub$gene, linetype = gdf_sub$branch, 
                                shape = gdf_sub$branch)) + geom_point() + labs(fill = "Gene", 
                                                                               col = "Gene", linetype = "Branch", shape = "Branch") + 
            theme_classic() + geom_smooth() + ylab("Expression") + 
            ggtitle(paste0(genes, collapse = ", ")) + NULL
        if (facet == "branch") {
            g = g + facet_grid(~branch)
        }
        if (facet == "gene") {
            g = g + facet_grid(gene ~ .)
        }
        if (facet == "both") {
            g = g + facet_grid(gene ~ branch)
        }
    }
    if (scHOT@positionType == "spatial") {
        if (is.null(positionColData)) {
            if (is.null(scHOT@positionColData)) {
                stop("Both positionColData and scHOT@positionColData are NULL.")
            }
            else {
                positionColData <- scHOT@positionColData
            }
        }
        coords_info <- data.frame(SummarizedExperiment::colData(scHOT)[, 
                                                                       positionColData])
        colnames(coords_info) <- positionColData
        branch_long <- reshape::melt(cbind(coords_info, t(branchData[genes, 
                                                                     , drop = FALSE])), id.vars = positionColData)
        colnames(branch_long) <- c("x", "y", "genes", "value")
        g <- ggplot(branch_long, aes(x = branch_long$x, y = branch_long$y, 
                                     color = branch_long$value)) + 
            geom_point(size = 0.5) + 
            theme_classic() + facet_wrap(~genes) + scale_alpha_continuous(range = c(0, 
                                                                                    0.5)) + scale_color_viridis_c(breaks = c(0, max(branch_long$value)), 
                                                                                                                  limits = c(0, max(branch_long$value)), labels = c("Low", 
                                                                                                                                                                    "High")) + theme(panel.grid = element_blank()) + 
            theme(axis.ticks = element_blank()) + theme(axis.text = element_blank()) + 
            xlab("") + ylab("") + coord_fixed() + labs(col = "Expression value") + 
            NULL
    }
    return(g)
}

weightedVarMeanSubtract = function(x, w) {
    weightedVariance(x - weightedMean(x,w), w = w)
}
