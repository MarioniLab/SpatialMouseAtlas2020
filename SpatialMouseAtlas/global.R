library(ggplot2)
library(patchwork)
library(ggpubr)
library(DelayedArray)
library(HDF5Array)
library(shinycssloaders)
library(ggiraph)
library(Matrix)

source("celltype_colours.R")

colours = c("gray75","cornflowerblue", "black")

# function to check if colours are legit colours
areColours <- function(x) {
    sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), 
                 error = function(e) FALSE)
    })
}
subsetColours = function(x) {
    # given a vector of colours remove the ones that aren't colours
    # give defaults if none or only one are given
    if (x[1] == "jonny") return(c("gray75", "cornflowerblue", "black"))
    if (x[1] == "default") return(c("yellow", "red"))
    xx = x[areColours(x)]
    if (length(xx) >= 2) return(xx)
    if (length(xx) == 1) return(c(xx, "white"))
    if (length(xx) == 0) return(c("white", "white"))
}

meta = readRDS("data/meta.Rds")
meta$selected = factor("Unselected", levels = c("Group A", "Group B", "Unselected"))

celltypes = readRDS("data/celltypes.Rds")
embryos = c("embryo1", "embryo2", "embryo3")
zvals = c(2,5)

genes = readRDS("data/genes.Rds")

embryo_coords_range_x = c(-3.57,3.09)
embryo_coords_range_y = c(-4.16,4.22)

file_combined = "data/combined_compressed.h5"
cnames = readRDS("data/combined_cnames.Rds")
rnames = readRDS("data/combined_rnames.Rds")
genes_imp = rnames

subsetCellsLogical = function(textToParse, input) {
    
    require(reshape)
    # meta, exprs, rnames, and imp are taken from the global environment
    # example inputs, remember quotes need to be escaped
    # textToParse = "Tbx18 < 0.5 & Shh > 0.5"
    # textToParse = "rank(x) < 10"
    # textToParse = "cluster == 4"
    # textToParse = "uniqueID %in% c(\"embryo3_Pos3_cell348_z2\")"
    
    # check if any are viable gene names
    textSplit = strsplit(textToParse, "[ \t\r\n]|=|<|>|!|\\(|\\)|&")[[1]]
    genesParse = intersect(textSplit, union(rnames, rownames(exprs)))
    print(paste(c("subsetting according to these genes:", genesParse), collapse = " "))
    
    # append these genes to the colData, prioritise measured over imputed
    colsToAddSeqFISH = as.matrix(exprs[intersect(genesParse, rownames(exprs)),,drop = FALSE])
    colsToAddImputed = as.matrix(imp[intersect(setdiff(genesParse, rownames(exprs)), rnames),,drop = FALSE])
    dfParse = cbind(meta, t(colsToAddSeqFISH), t(colsToAddImputed))
    
    dfParse_sub = subset(dfParse, embryo %in% input$embryo_subset & z %in% input$z_subset)
    
    selectedCells <- tryCatch(rownames(subset(dfParse_sub, eval(parse(text = textToParse)))), error = function(e) e)
    
    if (any(class(selectedCells) == "error")) {
        showNotification("Error in logical statement, e.g. gene not found, please try again...")
        return(NULL)
    }
    return(selectedCells)
}


embryolabeller = function(string) {
    newstring = string
    newstring[newstring == "embryo1"] <- "Embryo 1"
    newstring[newstring == "embryo2"] <- "Embryo 2"
    newstring[newstring == "embryo3"] <- "Embryo 3"
    return(newstring)
}

add_boundary_polygons = function() {
    if (!"boundary_polygons" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading cell segmentation...")
        boundary_polygons <<- readRDS("data/boundary_polygons.Rds")
        showNotification("Loading cell segmentation... done!")
    }
}

add_exprs = function() {
    if (!"exprs" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading expression...")
        exprs <<- readRDS("data/exprs.Rds")
        showNotification("Loading expression... done!")
    }
}

add_imp = function() {
    if (!"imp" %in% ls(envir = .GlobalEnv)) {
        showNotification("loading imputed data...")
        imp <<- HDF5Array(filepath = file_combined, name = "logcounts")
        rownames(imp) <<- rnames
        colnames(imp) <<- cnames
        showNotification("loading imputed data... done!")
    }
}

add_exprs_norm = function() {
    if (!"exprs_norm" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading batch-corrected expression...")
        exprs_norm <<- readRDS("data/exprs_norm.Rds")
        showNotification("Loading batch-corrected expression... done!")
    }
}

add_mRNA_df = function() {
    if (!"mRNA_df" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading mRNA data...")
        mRNA_df <<- readRDS("data/mRNA.Rds")
        showNotification("Loading mRNA data... done!")
    }
}

g_leg_list = readRDS("data/g_leg_list.Rds")