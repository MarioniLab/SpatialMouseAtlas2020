# function
# sort cells by their euclidean distance to the centroid for the type,
# and return the sorting quantile value
getSortedQuantile = function(coords,class) {
  df_split = split.data.frame(coords, class)
  quants_split = lapply(df_split, function(df_sub) {
    centroid = colMeans(df_sub)
    dists = colSums((t(df_sub) - centroid)^2)
    quants = rank(dists)/length(dists)
    names(quants) <- rownames(df_sub)
    return(quants)
  })
  names(quants_split) <- NULL
  quants_all = unlist(quants_split)[rownames(coords)]
}


# write function to output xspline data.frame
getClusterBorder = function(df, xname = "x", yname = "y", class = "class", subsetQuantile = 1, ...) {
  # ... passed to xspline
  
  # subset based on quantile value from centroid
  quants = getSortedQuantile(df[,c(xname, yname)], df[,class])
  df <- df[quants <= subsetQuantile, ]
  
  # take the convex hull
  chull.obj = lapply(split.data.frame(df, df[,class]), function(x)x[chull(x[,xname], x[,yname]),])
  
  # for some reason, call plot.new
  plot.new()
  
  # take xspline of the convex hull
  xspline.obj = lapply(chull.obj, xspline, draw = FALSE, open = FALSE, ...)
  
  # bind the x and y into two columns
  xspline.obj2 = lapply(xspline.obj, function(x) do.call(cbind, x))
  
  # convert to df with the class info
  xspline.df = data.frame(do.call(rbind, xspline.obj2),
                          class = rep(names(xspline.obj2), times = unlist(lapply(xspline.obj2, nrow))))
  colnames(xspline.df) <- c(xname, yname, class)
  
  return(xspline.df)
}
