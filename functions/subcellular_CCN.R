# write into a function, that takes in the data frame, and outputs the
# gene-gene graph and the pmat output of the cell-cell contact map
calculateSubcellularContactMap = function(
  mRNA_table,
  coordNames = c("x_global_affine","y_global_affine"),
  genesName = "geneID",
  distK = 5,
  nPerm = 1000,
  verbose = FALSE
) {
  
  require(bluster)
  require(igraph)
  require(SingleCellExperiment)
  
  # for a given cell, generate a NN network
  cell_NN = makeKNNGraph(mRNA_table[, coordNames], k = distK)
  V(cell_NN)$name <- rownames(mRNA_table)
  
  if (verbose) {
    print(mRNA_table[1,])
  }
  
  # generate a new SingleCellExperiment object
  mRNA_sce = SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 2, ncol = nrow(mRNA_table),
                                  dimnames = list(c("A", "B"), rownames(mRNA_table)))))
  mRNA_sce$geneID = droplevels(mRNA_table[,genesName])
  mRNA_sce$uniqueID = colnames(mRNA_sce)
  mRNA_sce$type = "mRNA"
  mRNA_sce$embryo = "mRNA"
  mRNA_sce$z = "mRNA"
  
  out = cellCellContact(sce = mRNA_sce, group = "geneID", cellID = "uniqueID",
                        graph = cell_NN, nperm = nPerm, plot = FALSE,
                        splitgroups = "type")
  
  return(list(graph = cell_NN, pmat = out["pmat"]))
}

split.data.frame.by = function(df, by) {
  # split a data.frame by providing the name of the column
  # of the data.frame to split by, rather than a factor
  split.data.frame(df, df[,by])
}