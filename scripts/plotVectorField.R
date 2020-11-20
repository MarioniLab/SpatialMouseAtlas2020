plotVectorField = function(x, y, val, nbins = 15, direction = 1, mask = FALSE, plot = TRUE) {
    
    # x is a numeric vector, e.g. UMAP1
    # y is a numeric vector, e.g. UMAP2
    # val is the outcome of interest to plot vector fields, e.g. log fold change
    # x,y,val need to have the same length, val can contain NAs
    # nbins is the number of bins for the overlaid vector field
    # direction is 1 or -1, if 1 it points from negative to positive
    # mask if TRUE will only plot vector fields for regions with points (problematic for sparse)
    # plot if TRUE will also output the ggplot object
    
    # output is a list containing
    # arrows_df a data.frame with x, y, magnitude and angle (in degrees)
    # plot is a ggplot object of the scatterplot with arrows overlaid
    
    
    # params
    # nbins = 80 # nice for spatial
    # nbins = 15
    # direction = 1 # one of -1 or 1, 1 means towards larger values
    # mask = FALSE # whether to plot arrows for all positions or only for sampled
    
    require(raster) # for raster manipulation
    require(fields) # for interpolation
    require(pracma) # for modulo
    
    if (plot) {
        require(ggplot2)
        require(metR) # for arrows geom
    }
    
    # inside function
    
    # get visual coords and outcome
    df = data.frame(x = x,
                    y = y,
                    val = val)
    
    xcut = seq(from = min(df$x), to = max(df$x), length.out = nbins + 1)
    ycut = seq(from = min(df$y), to = max(df$y), length.out = nbins + 1)
    
    m = tapply(df$val, list(cut(df$x, breaks = xcut, include.lowest = TRUE),
                            cut(df$y, breaks = ycut, include.lowest = TRUE)),
               mean, na.rm = TRUE)
    m <- m/max(abs(m), na.rm = TRUE)
    
    # generate raster object
    r=raster(t(m[,ncol(m):1]), xmn=0.5,xmx=nrow(m)+.5, ymn=0.5,ymx=ncol(m)+0.5)
    projection(r)=CRS("+init=epsg:27700")
    
    # interpolate regions with NAs
    {
        xy <- data.frame(xyFromCell(r, 1:ncell(r)))
        v <- getValues(r)
        # remove NAs
        i <- !is.na(v)
        xy <- xy[i,]
        v <- v[i]
    }
    
    #### Thin plate spline model
    tps <- Tps(xy, v)
    p <- raster(r)
    
    # use model to predict values at all locations
    p <- interpolate(p, tps)
    if (mask) {
        p <- mask(p, r)
    }
    
    # get vector field
    out = raster::terrain(p, opt = c("slope","aspect"), unit = "degrees")
    out_slope = t(as.matrix(out$slope)[(nrow(out$slope)-1):2, 2:(ncol(out$slope)-1)])
    out_aspect = t(as.matrix(out$aspect)[(nrow(out$aspect)-1):2, 2:(ncol(out$aspect)-1)])
    
    arrows_df = data.frame(x = rep(xcut[2:(length(xcut)-2)] + 0.5*(xcut[2]-xcut[1]),
                                   times = ncol(out_slope)),
                           y = rep(ycut[2:(length(ycut)-2)] + 0.5*(ycut[2]-ycut[1]),
                                   each = nrow(out_slope)),
                           magnitude = c(out_slope),
                           direction = c(out_aspect))
    if (direction == 1) {
        arrows_df$direction <- mod(-(arrows_df$direction + 90), 360)
    } else {
        arrows_df$direction <- mod(arrows_df$direction + 90, 360)
    }
    
    
    if (plot) {
        # plot
        g = ggplot(df, aes(x = x, y = y, fill = val)) + 
            geom_point(colour = "lightgrey", shape = 21) +
            theme_classic() + 
            scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                 na.value = "grey")
        # g
        
        # add the vectors
        g2 = g + 
            geom_arrow(aes(x = x, y = y, mag = magnitude, angle = direction), 
                       data = arrows_df, inherit.aes = FALSE,
                       colour = "black")
        g2
    } else {
        g2 = NULL
    }
    
    return(list(arrows_df = arrows_df, plot = g2))
}
