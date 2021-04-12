shinyServer(function(input, output, session) {
    
    values <- reactiveValues(meta = meta)
    
    gene_name_impGenerator = function() {
        if (input$gene_name_imp == "") {
            values$gene_name_imp_name <<- NA
            out = "Please type in a gene name"
            return(out)
        }
        
        ind = match(tolower(input$gene_name_imp), tolower(genes_imp))[1]
        if (is.na(ind)) {
            
            # genes_imp
            ind = grepl(tolower(input$gene_name_imp), tolower(genes_imp))
            
            # then plot for the first one
            values$gene_name_imp_name <<- genes_imp[ind][1]
            
            # display genes matching
            out = paste(c("Selecting first of multiple matching genes:",
                          sort(genes_imp[ind])[1:pmin(5,sum(ind))],
                          ifelse(sum(ind) > 5, "and more.", "")), collapse = " ")
        } else {
            values$gene_name_imp_name <<- genes_imp[ind][1]
            out = paste0("Selecting gene: ", values$gene_name_imp_name)
            return(out)
        }
        
        if (is.na(values$gene_name_imp_name)) {
            out = "no genes match input"
        }
        
        return(out)
    }
    
    
    output$gene_name_imp_parse_status <- renderText({
        out = gene_name_impGenerator()
        return(out)
    })
    
    
    spatialPlotGenerator = function() {
        
        add_boundary_polygons()
        
        meta_sub = subset(
            meta,
            embryo %in% input$embryo_subset & 
                z %in% input$z_subset
        )
        selectedNames = meta_sub$uniqueID
        
        # subset the boundary_polygons dataframe
        boundary_polygons_bkg = subset(
            boundary_polygons,
            uniqueID %in% selectedNames
        )
        
        if (!input$show_segmentation) {
            boundary_polygons_bkg <- subset(
                boundary_polygons_bkg,
                uniquePoint
            )
        }
        
        if (input$tabs == "Spatial plot (Imputed)") {
            
            add_imp()
            
            gene_name_impGenerator()
            
            validate(
                need(isolate(values$gene_name_imp_name) %in% genes_imp, "No valid gene selected.")
            )
            
            expr_imp = imp[isolate(values$gene_name_imp_name),]
            names(expr_imp) <- colnames(imp)
            expr <- expr_imp[rownames(meta)]
            
            pc_vals = expr[selectedNames]
            
            
            validate(
                need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
            )
            
            pc_cols = colorRampPalette(colours)(100)[as.numeric(cut(pc_vals,100))]
            names(pc_cols) <- names(pc_vals)
            
            yl = " (Imputed expression)"
            
        } else {
            
            
            if (input$colour_by == "Batch-corrected gene expression") {
                
                add_exprs_norm()
                
                pc_vals = exprs_norm[input$gene_name, selectedNames]
                
                validate(
                    need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
                )
                
                pc_cols = colorRampPalette(colours)(100)[as.numeric(cut(pc_vals,100))]
                names(pc_cols) <- names(pc_vals)
                
                yl = ""
            }
            if (input$colour_by == "Mapped cell type") {
                
                pc_vals = as.character(meta_sub$celltype_mapped_refined)
                
                validate(
                    need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
                )
                
                pc_cols = celltype_colours[pc_vals]
                names(pc_cols) <- rownames(meta_sub)
                
            }
            
        }
        
        
        # subset further (need to update this)
        if (input$tabs == "Spatial plot (Imputed)") {
            pc_cols[rownames(meta)[!meta$celltype_mapped_refined %in% input$celltype_subset_imp]] <- "grey95"
        } else {
            pc_cols[rownames(meta)[!meta$celltype_mapped_refined %in% input$celltype_subset]] <- "grey95"
        }
        
        pc_cols_border = ifelse(!is.na(pc_cols), "black", "grey90")
        
        g_base = ggplot(boundary_polygons_bkg, aes(x = segmentation_vertices_x_global_affine,
                                                   y = segmentation_vertices_y_global_affine_neg,
                                                   group = uniqueID,
                                                   fill = uniqueID,
                                                   tooltip = celltype_mapped_refined,
                                                   colour = uniqueID))
        
        
        if (input$show_segmentation) {
            g = g_base +
                geom_polygon_interactive(
                    size = 0.05
                )
        } else {
            g = g_base +
                geom_point_interactive(
                    colour = "grey20",
                    size = 0.75,
                    shape = 21
                )
        }
        
        g <- g + 
            scale_fill_manual(values = pc_cols) +
            scale_colour_manual(values = pc_cols_border) +
            theme_classic() +
            coord_fixed() +
            xlab("") +
            # legend needs to be "none" because of way colour is defined 
            theme(legend.position = "none") + 
            theme(axis.line = element_blank()) + 
            theme(axis.text = element_blank()) + 
            theme(axis.title.y = element_blank()) +
            theme(axis.title.x = element_text(face = "italic")) +
            theme(axis.ticks = element_blank()) +
            theme(strip.text.y = element_blank()) +
            theme(strip.background.x = element_rect(colour = "white")) +
            NULL
        
        if (length(input$embryo_subset) > 1) {
            
            if (input$show_segmentation) {
                g <- g +
                    facet_grid(z~embryo, labeller = labeller(embryo = embryolabeller)) +
                    NULL
            } else {
                g <- g +
                    facet_grid(~embryo, labeller = labeller(embryo = embryolabeller)) +
                    NULL
            }
        } else {
            if (length(input$z_subset) > 1) {
                if (input$show_segmentation) {
                    g <- g +
                        facet_grid(z~embryo, labeller = labeller(embryo = embryolabeller)) +
                        NULL
                } else {
                    g <- g +
                        facet_grid(~embryo, labeller = labeller(embryo = embryolabeller)) +
                        NULL
                }
            }
        }
        
        g_leg = NULL
        
        if (input$colour_by %in% c("Area-standardised gene expression", "Batch-corrected gene expression") | input$tabs == "Spatial plot (Imputed)") {
            
            heights = c(10,1)
            
            if (input$tabs == "Spatial plot (Imputed)") {
                gene_name_impGenerator()
                validate(
                    need(isolate(values$gene_name_imp_name) %in% genes_imp, "No valid gene selected.")
                )
                g <- g + 
                    xlab(paste0(isolate(values$gene_name_imp_name), yl)) + 
                    NULL
                ttl = "Imputed gene expression"
            } else {
                g <- g + 
                    xlab(paste0(input$gene_name, yl)) + 
                    NULL
                ttl = input$colour_by
            }
            
            g_leg = g_leg_list[[ttl]]
        }
        
        if (input$colour_by == "Mapped cell type" & input$tabs != "Spatial plot (Imputed)") {
            
            heights = c(10,2)
            
            # get the legend
            g_leg <- g_leg_list[["Mapped cell type"]]
        }
        
        
        if (input$colour_by == "Mapped gut tube subtype" & input$tabs != "Spatial plot (Imputed)") {
            
            heights = c(10,1)
            
            # get the legend
            pc_df = data.frame(Expression = factor(unique(pc_vals), levels = names(AP_pseudo_colours)))
            
            g_leg_raw = ggplot(pc_df, aes(x = Expression, y = Expression)) +
                geom_point(aes(colour = Expression), size = 10) +
                theme_transparent() +
                scale_colour_manual(values = AP_pseudo_colours) +
                theme(legend.position = "bottom") +
                theme(legend.text = element_text(size = 15),
                      legend.title = element_text(size = 20)) +
                guides(color = guide_legend(title.position = "top",
                                            title.hjust = 0.5,
                                            ticks = FALSE,
                                            title = "Mapped gut tube cell type")) +
                NULL
            
            g_leg = as_ggplot(get_legend(g_leg_raw))
        }
        
        
        list(g = g, g_leg = g_leg)
    }
    
    output$spatialPlot <- renderGirafe({
        
        
        g = spatialPlotGenerator()[["g"]]
        
        gi = girafe(code = print(g),
                    options = list(
                        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                    ))
        gi
        
        return(gi)
    })
    
    
    output$spatialPlot_imp <- renderGirafe({
        
        g = spatialPlotGenerator()[["g"]]
        
        gi = girafe(code = print(g),
                    options = list(
                        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                    ))
        gi
        
        return(gi)
    })
    
    output$spatialPlot_leg <- renderPlot({
        return(spatialPlotGenerator()[["g_leg"]])
    })
    
    output$spatialPlot_leg_imp <- renderPlot({
        return(spatialPlotGenerator()[["g_leg"]])
    })
    
    
    
    
    violinPlotGenerator = function() {
        
        meta_sub = subset(
            meta,
            embryo %in% input$embryo_subset & 
                z %in% input$z_subset
        )
        selectedNames = rownames(meta_sub)
        
        if (input$tabs == "Spatial plot (Imputed)") {
            
            add_imp()
            
            # example gene
            # gene = "Pcdh19"
            # expr_imp = imp[input$gene_name_imp,]
            gene_name_impGenerator()
            validate(
                need(isolate(values$gene_name_imp_name) %in% genes_imp, "No valid gene selected.")
            )
            expr_imp = imp[isolate(values$gene_name_imp_name),]
            names(expr_imp) <- colnames(imp)
            expr <- expr_imp[rownames(meta)]
            
            pc_vals = expr[rownames(meta) %in% selectedNames]
            
            
            validate(
                need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
            )
            
            yl = paste0(isolate(values$gene_name_imp_name), " (Imputed expression)")
            
        } else {
            
            add_exprs_norm()
            
            pc_vals = exprs_norm[input$gene_name, ][colnames(exprs_norm) %in% selectedNames]
            validate(
                need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
            )
            
            yl = input$gene_name
            
            
        }
        
        ct_df = data.frame(pc_val = pc_vals, ct = meta[names(pc_vals), "celltype_mapped_refined"])
        ct_count = unclass(table(ct_df$ct))
        ct_df$ct_count = ct_count[as.character(ct_df$ct)]
        
        
        g_violin = ggplot(ct_df, aes(x = ct,
                                     y = pc_val,
                                     group = ct,
                                     fill = ct
        )) + 
            geom_text(aes(label = ct_count), angle = 60, data = data.frame(ct_count = ct_count, ct = names(ct_count), pc_val = 0.5 + max(na.omit(ct_df$pc_val)))) +
            geom_violin(draw_quantiles = 0.5, colour = "black", scale = "width") + 
            theme_classic() + 
            scale_fill_manual(values = celltype_colours) + 
            xlab("") +
            ylab(yl) + 
            theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
            theme(axis.title.y = element_text(face = "italic")) + 
            theme(legend.position = "none") +
            NULL
        
        g_violin
    }
    
    output$violinPlot <- renderGirafe({
        
        
        g = violinPlotGenerator()
        
        gi = girafe(code = print(g),
                    options = list(
                        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                    ))
        gi
        
        return(gi)
    })
    
    output$violinPlot_imp <- renderGirafe({
        
        
        g = violinPlotGenerator()
        
        gi = girafe(code = print(g),
                    options = list(
                        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                    ))
        gi
        
        return(gi)
    })
    
    
    mRNAPlotGenerator = function() {
        
        validate(need(nchar(input$colours_mRNA) > 0,
                      "Please type in colours (separated by spaces) for digital in situ."
        ))
        validate(need(length(input$gene_names_mRNA)>0,
                      "Please select some genes for digital in situ."        
        ))
        
        add_boundary_polygons()
        
        gene_names_mRNA <- input$gene_names_mRNA[1:pmin(5,length(input$gene_names_mRNA))]
        gene_colours_mRNA <- subsetColours(unlist(strsplit(input$colours_mRNA, " ")))
        if (length(gene_colours_mRNA) < length(gene_names_mRNA)) {
            gene_colours_mRNA <- c(gene_colours_mRNA, rep("white", length(gene_names_mRNA) - length(gene_colours_mRNA)))
        }
        
        names(gene_colours_mRNA) <- gene_names_mRNA
        
        embryo2_shift_x = -(input$embryo2_centre_x - input$embryo1_centre_x)
        embryo3_shift_x = -(input$embryo3_centre_x - input$embryo1_centre_x)
        
        embryo2_shift_y = -(input$embryo2_centre_y - input$embryo1_centre_y)
        embryo3_shift_y = -(input$embryo3_centre_y - input$embryo1_centre_y)
        
        heart_shift_x = c(0, embryo2_shift_x, embryo3_shift_x)
        names(heart_shift_x) <- c("embryo1","embryo2","embryo3")
        
        heart_shift_y = c(0, embryo2_shift_y, embryo3_shift_y)
        names(heart_shift_y) <- c("embryo1","embryo2","embryo3")
        
        range_x = c(input$embryo1_centre_x - input$radius, input$embryo1_centre_x + input$radius)
        range_y = c(input$embryo1_centre_y - input$radius, input$embryo1_centre_y + input$radius)
        add_mRNA_df()
        
        mRNA_genes = mRNA_df[mRNA_df$geneID %in% gene_names_mRNA &
                                 mRNA_df$embryo %in% input$embryo_subset & 
                                 mRNA_df$z %in% input$z_subset,]
        
        boundary_polygons_sub = subset(boundary_polygons,
                                       boundary_polygons$embryo %in% input$embryo_subset & 
                                           boundary_polygons$z %in% input$z_subset)
        
        if (!input$mRNA_full) {
            
            mRNA_genes$x_global_affine <- mRNA_genes$x_global_affine + heart_shift_x[mRNA_genes$embryo]
            mRNA_genes$y_global_affine <- mRNA_genes$y_global_affine - heart_shift_y[mRNA_genes$embryo]
            
            boundary_polygons_sub$segmentation_vertices_x_global_affine <- boundary_polygons_sub$segmentation_vertices_x_global_affine + heart_shift_x[boundary_polygons_sub$embryo]
            boundary_polygons_sub$segmentation_vertices_y_global_affine_neg <- boundary_polygons_sub$segmentation_vertices_y_global_affine_neg + heart_shift_y[boundary_polygons_sub$embryo]
            
        }
        
        validate(
            need(nrow(boundary_polygons_sub) > 0,
                 "No cells are selected for subsetting, please tick at least one option for each of the embryo and z-slice categories.")
        )
        
        
        g = ggplot(boundary_polygons_sub,
                   aes(x = segmentation_vertices_x_global_affine,
                       y = segmentation_vertices_y_global_affine_neg,
                       group = uniqueID),
                   fill = "lightgrey",
                   size = 0.1) + 
            geom_polygon() + 
            
            geom_polygon(data = subset(boundary_polygons_sub, 
                                       celltype_mapped_refined %in% input$celltype_outline_mRNA),
                         colour = "white",
                         size = 0.1) +
            
            theme_classic() +
            
            facet_grid(.~embryo, labeller = labeller(embryo = embryolabeller)) +
            guides(fill = guide_legend(title = "normalised logcounts")) +
            coord_fixed() +
            xlab("") +
            ylab("") +
            theme(plot.background = element_rect(fill = "black"),
                  panel.background = element_rect(fill = "black"),
                  strip.background = element_rect(fill = "black"),
                  strip.text = element_text(colour = "white"),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  plot.title = element_text(colour = "white", hjust = 0.5)
            ) +
            theme(legend.background = element_rect(fill = "black"),
                  legend.text = element_text(colour = "white", face = "italic", size = 15)) +
            NULL
        
        if (!input$mRNA_full) {
            g <- g + 
                xlim(range_x) + 
                ylim(range_y) +
                NULL
        }
        
        g_dots = g + 
            geom_point(aes(x = x_global_affine, y = -y_global_affine,
                           colour = geneID), 
                       data = mRNA_genes,
                       size = 0.01,
                       inherit.aes = FALSE) + 
            scale_colour_manual(values = gene_colours_mRNA) +
            guides(colour = guide_legend(title = "Gene",
                                         override.aes = list(size = 10))) +
            NULL
        
        return(g_dots)
    }
    
    g_mRNA <- eventReactive(input$go, {
        mRNAPlotGenerator()
    })
    
    output$mRNAPlot <- renderPlot({
        g_mRNA_val = g_mRNA()
        return(g_mRNA_val)
    })
    
    mRNARegionPlotGenerator = function() {
        
        polydf = data.frame(
            embryo = rep(c("embryo1", "embryo2","embryo3"), each = 4),
            x_global_affine = c(input$embryo1_centre_x - input$radius,
                                input$embryo1_centre_x + input$radius,
                                input$embryo1_centre_x + input$radius,
                                input$embryo1_centre_x - input$radius,
                                
                                input$embryo2_centre_x - input$radius,
                                input$embryo2_centre_x + input$radius,
                                input$embryo2_centre_x + input$radius,
                                input$embryo2_centre_x - input$radius,
                                
                                input$embryo3_centre_x - input$radius,
                                input$embryo3_centre_x + input$radius,
                                input$embryo3_centre_x + input$radius,
                                input$embryo3_centre_x - input$radius
            ),
            y_global_affine = -c(input$embryo1_centre_y - input$radius,
                                 input$embryo1_centre_y - input$radius,
                                 input$embryo1_centre_y + input$radius,
                                 input$embryo1_centre_y + input$radius,
                                 
                                 input$embryo2_centre_y - input$radius,
                                 input$embryo2_centre_y - input$radius,
                                 input$embryo2_centre_y + input$radius,
                                 input$embryo2_centre_y + input$radius,
                                 
                                 input$embryo3_centre_y - input$radius,
                                 input$embryo3_centre_y - input$radius,
                                 input$embryo3_centre_y + input$radius,
                                 input$embryo3_centre_y + input$radius
            )
        )
        
        meta_sub = subset(meta, embryo %in% input$embryo_subset)
        
        validate(
            need(nrow(meta_sub) > 0,
                 "No cells are selected for subsetting, please tick at least one option for each of the embryo and z-slice categories.")
        )
        
        
        g = ggplot(meta_sub,
                   aes(x = x_global_affine,
                       y = -y_global_affine)) + 
            geom_point(size = 0.1, 
                       aes(colour = celltype_mapped_refined %in% input$celltype_outline_mRNA)) + 
            
            facet_wrap(~embryo, nrow = 1, labeller = labeller(embryo = embryolabeller)) + 
            theme_classic() + 
            coord_fixed() +
            xlab("") + 
            ylab("") +
            ggtitle("Regional selection") +
            scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey")) +
            theme(legend.position = "none") +
            theme(plot.title = element_text(size = 15, hjust = 0.5)) +
            NULL
        
        if (!input$mRNA_full) {
            g <- g + 
                geom_polygon(data = subset(polydf, embryo %in% input$embryo_subset),
                             aes(group = embryo),
                             colour = "red", fill = NA, size = 1) +
                NULL
        }
        
        return(g)
    }
    
    output$mRNARegionPlot <- renderPlot({
        
        g_mRNARegion_val = mRNARegionPlotGenerator()
        
        return(g_mRNARegion_val)
    })
    
    output$download_mRNAPlot <- downloadHandler(
        filename = function() {
            return("digital_in_situ.pdf")
        },
        content = function(file) {
            m1 = mRNARegionPlotGenerator()
            m2 = g_mRNA()
            pdf(file, paper="a4")
            if (!input$mRNA_full) {
                print(m1)
            }
            print(m2)
            dev.off()
        }
    )
    
    output$download_mRNAPlot_gg <- downloadHandler(
        filename = function() {
            return("digital_in_situ_gg.Rds")
        },
        content = function(file) {
            m2 = g_mRNA()
            saveRDS(m2, file)
        }
    )
    
    # for Virtual Dissection section
    VirtualDissectionPlotGenerator = function() {
        
        meta2 <- isolate(values$meta)
        
        meta_sub = subset(meta2,
                          meta2$embryo %in% input$embryo_subset & 
                              meta2$z %in% input$z_subset)
        
        validate(
            need(nrow(meta_sub) > 0,
                 "No cells are selected for subsetting, please tick at least one option for each of the embryo and z-slice categories.")
        )
        
        g = ggplot(meta_sub, 
                   aes(x = x_global_affine,
                       y = -y_global_affine,
                       tooltip = celltype_mapped_refined,
                       data_id = uniqueID)) + 
            geom_point_interactive(size = 0.1, aes(colour = selected)) + 
            facet_grid(~embryo, labeller = labeller(embryo = embryolabeller)) +
            scale_colour_manual(values = c("Group A" = "red",
                                           "Unselected" = "grey",
                                           "Group B" = "blue")) +
            theme_classic() + 
            coord_fixed() + 
            theme(legend.position = "bottom") + 
            theme(axis.line = element_blank()) + 
            theme(axis.text = element_blank()) + 
            theme(axis.title = element_blank()) +
            theme(axis.ticks = element_blank()) +
            theme(strip.background = element_rect(colour = "white")) +
            guides(colour = guide_legend(title = "",
                                         override.aes = list(size = 5))) +
            NULL
        
        gi = girafe(code = print(g),
                    options = list(
                        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                    ))
        gi
    }
    
    VirtualDissectionPlot = eventReactive(input$replotVirtualDissection,
                                          VirtualDissectionPlotGenerator()
    )
    
    output$virtualDissection <- renderGirafe({
        
        VirtualDissectionPlot = VirtualDissectionPlot()
        VirtualDissectionPlot
    })
    
    observeEvent(input$virtualDissection_reset, {
        session$sendCustomMessage(type = 'virtualDissection_set', message = character(0))
    })
    
    virtualDissection_selected_state <- reactive({
        input$virtualDissection_selected
    })
    
    addGenerator = function() {
        meta2 <- isolate(values$meta)
        meta2[meta2$uniqueID %in% virtualDissection_selected_state(),"selected"] <- "Group A"
        values$meta <<- meta2
        showNotification("Virtually dissected cells added to Group A! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    
    observeEvent(input$add, {
        addGenerator()
    })
    
    add2Generator = function() {
        meta2 <- isolate(values$meta)
        meta2[meta2$uniqueID %in% virtualDissection_selected_state(),"selected"] <- "Group B"
        values$meta <<- meta2
        showNotification("Virtually dissected cells added to Group B! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$add2, {
        add2Generator()
    })
    
    removeGenerator = function() {
        meta2 <- isolate(values$meta)
        meta2[meta2$uniqueID %in% virtualDissection_selected_state(),"selected"] <- "Unselected"
        values$meta <<- meta2
        showNotification("Virtually dissection cells removed from either group! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$remove, {
        removeGenerator()
    })
    
    removeAllGenerator = function() {
        meta2 <- isolate(values$meta)
        meta2[,"selected"] <- "Unselected"
        values$meta <<- meta2
        
        showNotification("All cells removed from selection! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$removeAll, {
        removeAllGenerator()
    })
    
    ###### CT start
    
    addCTGenerator = function() {
        
        meta2 <- isolate(values$meta)
        selectedUniqueIDs = subset(meta2,
                                   celltype_mapped_refined
                                   %in% input$virtualDissectionCelltypes & 
                                       embryo %in% input$embryo_subset & 
                                       z %in% input$z_subset)[,"uniqueID"]
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group A"
        values$meta <<- meta2
        showNotification("Selected cell type cells added to Group A! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$addCT, {
        addCTGenerator()
    })
    
    addCT2Generator = function() {
        meta2 <- isolate(values$meta)
        selectedUniqueIDs = subset(meta2,
                                   celltype_mapped_refined
                                   %in% input$virtualDissectionCelltypes & 
                                       embryo %in% input$embryo_subset & 
                                       z %in% input$z_subset)[,"uniqueID"]
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group B"
        values$meta <<- meta2
        showNotification("Selected cell type cells added to Group B! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$addCT2, {
        addCT2Generator()
    })
    
    removeCTGenerator = function() {
        meta2 <- isolate(values$meta)
        selectedUniqueIDs = subset(meta2,
                                   celltype_mapped_refined
                                   %in% input$virtualDissectionCelltypes & 
                                       embryo %in% input$embryo_subset & 
                                       z %in% input$z_subset)[,"uniqueID"]
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Unselected"
        values$meta <<- meta2
        showNotification("Selected cell type cells removed from either group! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$removeCT, {
        removeCTGenerator()
    })
    
    ######### CT end
    
    
    
    ###### Logical start
    
    addLogicalGenerator = function() {
        meta2 <- isolate(values$meta)
        selectedUniqueIDs = subsetCellsLogical(input$logical, input = input)
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group A"
        values$meta <<- meta2
        showNotification("Logical statement selected cells added to Group A! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$addLogical, {
        
        add_exprs()
        add_imp()
        
        addLogicalGenerator()
    })
    
    addLogical2Generator = function() {
        meta2 <- isolate(values$meta)
        selectedUniqueIDs = subsetCellsLogical(input$logical, input = input)
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group B"
        values$meta <<- meta2
        showNotification("Logical statement selected cells added to Group B! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$addLogical2, {
        
        add_exprs()
        add_imp()
        
        addLogical2Generator()
    })
    
    removeLogicalGenerator = function() {
        
        add_exprs()
        add_imp()
        
        meta2 <- isolate(values$meta)
        selectedUniqueIDs = subsetCellsLogical(input$logical, input = input)        
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Unselected"
        values$meta <<- meta2
        showNotification("Logical statement selected cells removed from either group! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$removeLogical, {
        removeLogicalGenerator()
    })
    
    ######### Logical end
    
    
    
    
    
    addPreselectedGenerator = function() {
        selectedUniqueIDs = do.call(c,
                                    sapply(input$preselected$datapath, 
                                           readRDS, simplify = FALSE)
        )
        meta2 <- isolate(values$meta)
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group A"
        values$meta <<- meta2
        showNotification("Preselected cells added to Group A! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$add_pre, {
        addPreselectedGenerator()
    })
    
    addPreselected2Generator = function() {
        selectedUniqueIDs = do.call(c,
                                    sapply(input$preselected$datapath, 
                                           readRDS, simplify = FALSE)
        )
        meta2 <- isolate(values$meta)
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group B"
        values$meta <<- meta2
        showNotification("Preselected cells added to Group B! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$add_pre2, {
        addPreselected2Generator()
    })
    
    removePreselectedGenerator = function() {
        selectedUniqueIDs = do.call(c,
                                    sapply(input$preselected$datapath, 
                                           readRDS, simplify = FALSE)
        )
        meta2 <- isolate(values$meta)
        meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Unselected"
        values$meta <<- meta2
        showNotification("Preselected cells removed from either group! Click Load/Reload plot to reload.")
        return(meta)
    }
    
    observeEvent(input$remove_pre, {
        removePreselectedGenerator()
    })
    
    output$virtualDissectionCellDownload <- downloadHandler(
        filename = function() {
            paste("cells_groupA.Rds")
        },
        content = function(file) {
            meta2 <- isolate(values$meta)
            meta_sub = subset(meta2, selected == "Group A" &
                                  embryo %in% input$embryo_subset & 
                                  z %in% input$z_subset)
            saveRDS(sort(as.character(unique(meta_sub[,"uniqueID"]))), file)
        }
    )
    
    output$virtualDissectionCellDownload2 <- downloadHandler(
        filename = function() {
            paste("cells_groupB.Rds")
        },
        content = function(file) {
            meta2 <- isolate(values$meta)
            meta_sub = subset(meta2, selected == "Group B" &
                                  embryo %in% input$embryo_subset & 
                                  z %in% input$z_subset)
            saveRDS(sort(as.character(unique(meta_sub[,"uniqueID"]))), file)
        }
    )
    
    VirtualDissectionBarPlotGenerator = function() {
        
        meta2 <- isolate(values$meta)
        
        meta_subList = sapply(c("Group A", "Group B"), function(sel) {
            subset(meta2, selected == sel &
                       embryo %in% input$embryo_subset &
                       z %in% input$z_subset)
        }, simplify = FALSE)
        
        meta_subList <- meta_subList[unlist(lapply(meta_subList, nrow)) > 0]
        
        gList = lapply(meta_subList, function(meta_sub) {
            
            if (nrow(meta_sub) == 0) return(NULL)
            
            meta_sub_sum = data.frame(
                count = tapply(meta_sub$celltype_mapped_refined, meta_sub$celltype_mapped_refined, length))
            meta_sub_sum$celltype_mapped_refined <- rownames(meta_sub_sum)
            meta_sub_sum = na.omit(meta_sub_sum)
            
            g = ggplot(meta_sub, 
                       aes(x=reorder(celltype_mapped_refined,celltype_mapped_refined,
                                     function(x)-length(x)))) +    
                geom_bar(aes(fill = celltype_mapped_refined)) + 
                geom_text(aes(y = count,
                              label = count), 
                          data = meta_sub_sum, 
                          vjust = "bottom",
                          size = 5) +
                scale_fill_manual(values = celltype_colours) + 
                theme_classic() + 
                theme(legend.position = "none") +
                theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 15)) +
                theme(axis.title.y = element_text(size = 15)) +
                theme(plot.title = element_text(size = 15)) +
                xlab("") +
                ylab("Number of cells selected") +
                ggtitle(paste0("Mapped cell types for ", 
                               sum(meta_sub_sum$count), 
                               " selected cells")) + 
                scale_x_discrete(expand = c(0.1, 0.15)) +
                scale_y_continuous(expand = c(0.1, 0.2)) +
                NULL
            return(g)
        })
        
        if (all(unlist(lapply(gList, is.null)))) return(NULL)
        
        gList_named <- sapply(names(gList), function(gname){
            gList[[gname]] + ylab(paste0("Number of ", gname, " cells"))
        }, simplify = FALSE)
        
        wrap_plots(gList_named, ncol = 1)
    }
    
    VirtualDissectionBarPlot = eventReactive(input$replotVirtualDissection,
                                             VirtualDissectionBarPlotGenerator()
    )
    
    output$VirtualDissectionBarPlot <- renderPlot({
        VirtualDissectionBarPlot = VirtualDissectionBarPlot()
        
        validate(
            need(!is.null(VirtualDissectionBarPlot),
                 "Need cells assigned to either Group A or Group B to be able to display barplot")
        )
        
        VirtualDissectionBarPlot
    })
    
    
    VirtualDissectionMAPlotGenerator = function() {
        
        
        meta2 <- isolate(values$meta)
        
        virtualDissectionCells = as.character(
            subset(meta2, 
                   selected == "Group A" &
                       embryo %in% input$embryo_subset & 
                       z %in% input$z_subset)[,"uniqueID"])
        
        validate(
            need(length(virtualDissectionCells) > 0,
                 "Need cells assigned to both Group A and Group B to be able to display MA-plot"
            )
        )
        
        
        
        virtualDissectionCells_not = as.character(
            subset(meta2, 
                   selected == "Group B" &
                       embryo %in% input$embryo_subset & 
                       z %in% input$z_subset)[,"uniqueID"])
        
        validate(
            need(length(virtualDissectionCells_not) > 0,
                 "Need cells assigned to both Group A and Group B to be able to display MA-plot")
        )
        
        add_exprs()
        
        
        
        MA_df_imp = NULL
        
        A_val = rowMeans(exprs[, c(virtualDissectionCells, virtualDissectionCells_not)])
        M_val = rowMeans(exprs[, virtualDissectionCells]) - rowMeans(exprs[, virtualDissectionCells_not])
        MA_df_seq = data.frame(
            A_val = A_val,
            M_val = M_val,
            gene = rownames(exprs),
            M_rank = rank(-M_val),
            M_rank_bottom = rank(M_val),
            type = "seqFISH"
        )
        
        MA_df = rbind(MA_df_seq, MA_df_imp)
        MA_df$type <- factor(MA_df$type, levels = c("seqFISH", "Imputed"))
        
        
        g = ggplot(MA_df, aes(x = A_val, y = M_val, label = gene, tooltip = gene, data_id = gene)) + 
            geom_point_interactive(colour = "grey") + 
            geom_hline(yintercept = 0) + 
            geom_text_interactive(aes(label = gene, colour = "red"),
                                  data = subset(MA_df, M_rank <= 5),
                                  size = 5,
                                  fontface = "italic") +
            geom_text_interactive(aes(label = gene, colour = "blue"),
                                  data = subset(MA_df, M_rank_bottom<= 5),
                                  size = 5,
                                  fontface = "italic") +
            theme_classic() +
            xlab("Mean expression of Group A and Group B cells") + 
            ylab("Difference in means of Group A versus Group B cells") +
            ggtitle("MA-Plot of Area-standardised gene expression\n(not batch-corrected)") +
            
            scale_color_identity(name = "",
                                 breaks = c("red", "blue"),
                                 labels = c("Higher in Group A", "Higher in Group B"),
                                 guide = "legend") +
            theme(legend.position = "none") +
            theme(legend.text = element_text(size = 15)) +
            
            guides(
                colour = guide_legend(
                    title = "",
                    override.aes = aes(label = "X",
                                       size = 10)
                )) +
            
            theme(axis.title.y = element_text(size = 15)) +
            theme(axis.title.x = element_text(size = 15)) +
            facet_grid(type~.) +
            NULL
        g
        
    }
    
    VirtualDissectionMAPlot = eventReactive(input$replotVirtualDissection,
                                            VirtualDissectionMAPlotGenerator()
    )
    
    # comment testing ggtips
    output$VirtualDissectionMAPlot <- renderGirafe({
        VirtualDissectionMAPlot = VirtualDissectionMAPlot()
        girafe(code = print(VirtualDissectionMAPlot),
               height_svg = 10,
               options = list(
                   opts_hover(css = "fill:#000000;stroke:black;cursor:pointer;", reactive = TRUE),
                   opts_selection(
                       type = "none", css = "fill:#000000;stroke:black;")
               ))
    })
    
    mRNARegionPlot_dblclickGenerator = function(coord = "mRNARegionPlot_dblclick") {
        
        dblclick_df = nearPoints(meta,
                                 coordinfo = input[[coord]],
                                 xvar = "x_global_affine",
                                 yvar = "y_global_affine",
                                 maxpoints = 1
        )
        uniq = dblclick_df[,"uniqueID"]
        
        observeEvent(input$mRNARegionPlot_dblclick, {
            newval_x = meta[uniq,"x_global_affine"]# - (input$embryo1_centre_x + radius)
            newval_y = meta[uniq,"y_global_affine"]# - (-input$embryo1_centre_y + radius)
            updateSliderInput(session,
                              paste0(meta[as.character(uniq), "embryo"],
                                     "_centre_x"),
                              value = newval_x,
                              min = embryo_coords_range_x[1],
                              max = embryo_coords_range_x[2]
            )
            updateSliderInput(session,
                              paste0(meta[as.character(uniq), "embryo"],
                                     "_centre_y"),
                              value = newval_y,
                              min = embryo_coords_range_y[1],
                              max = embryo_coords_range_y[2]
            )
        })
        
        
        return(uniq)
        
    }
    
    
    # this is here so that double click in digital in situ works
    output$info <- renderText({
        xy_str <- function(e) {
            if(is.null(e)) return("NULL\n")
            paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
        }
        
        paste0(
            "dblclick: ", xy_str(input$mRNARegionPlot_dblclick)
        )
        
        uniq = mRNARegionPlot_dblclickGenerator()
        paste0("")
    })
    
})
