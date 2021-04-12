library(shinyjs)
library(shinythemes)

shinyUI(fluidPage(
    useShinyjs(),
    theme=shinytheme("cerulean"),
    fluidRow(
        column(12,
               titlePanel("Highly multiplexed spatially resolved gene expression profiling of mouse organogenesis"))
    ),
    
    fluidRow(
        column(1,
               
               checkboxGroupInput("embryo_subset",
                                  "Subset by embryo:",
                                  choiceNames = c("Embryo 1", "Embryo 2", "Embryo 3"),
                                  choiceValues = c("embryo1", "embryo2", "embryo3"),
                                  selected = "embryo1"),
               
               HTML("z-slices are contiguous cell layers 12 um apart."),
               
               checkboxGroupInput("z_subset",
                                  "Subset by z-slice:",
                                  choiceNames = c("Slice 1", "Slice 2"),
                                  choiceValues = zvals,
                                  selected = 2),
               
               checkboxInput('show_segmentation',
                             'Display cell segmentation in Spatial plots',
                             value = FALSE)
               
        ),
        column(11,
               
               tabsetPanel(id = "tabs",
                           tabPanel("Landing page", fluid = TRUE,
                                    includeMarkdown("README.md")
                           ),
                           
                           tabPanel("Spatial plot", fluid = TRUE,
                                    
                                    column(2,
                                           
                                           selectInput("colour_by",
                                                       "Colour by:",
                                                       choices = c("Batch-corrected gene expression",
                                                                   "Mapped cell type"
                                                       ),
                                                       selected = "Batch-corrected gene expression",
                                                       multiple = FALSE),
                                           
                                           selectizeInput("gene_name",
                                                          "Gene name:",
                                                          choices = genes,
                                                          selected = "Ttn"),
                                           
                                           selectizeInput("celltype_subset",
                                                          "Subset by mapped cell type:",
                                                          choices = celltypes,
                                                          selected = setdiff(celltypes, "Low quality"),
                                                          multiple = TRUE),
                                           
                                    ),
                                    
                                    # Show a plot of the generated distribution
                                    column(8,
                                           mainPanel(
                                               
                                               HTML("<p style=\"text-align: center;\"><em>Plot of spatial expression or cell types, and distribution of expression per cell type (below). Check option on left panel for plot of cells' segmentation. Use the icon in the top-right corner to download each image.</em></p>"),
                                               
                                               girafeOutput("spatialPlot",
                                                            height = "900px",
                                                            width = "1200px"
                                               ) %>% withSpinner(color="#0dc5c1"),
                                               
                                               plotOutput("spatialPlot_leg",
                                                          height = "200px",
                                                          width = "1200px"),
                                               
                                               girafeOutput("violinPlot",
                                                            height = "600px",
                                                            width = "1200px"
                                               ) %>% withSpinner(color="#0dc5c1")
                                           )
                                    )
                           ),
                           tabPanel("Spatial plot (Imputed)", fluid = TRUE,
                                    
                                    column(2,
                                           
                                           textInput("gene_name_imp",
                                                     "Please type in a gene name:",
                                                     value = "Ttn"),
                                           
                                           textOutput("gene_name_imp_parse_status",
                                                      inline = TRUE),
                                           
                                           selectizeInput("celltype_subset_imp",
                                                          "Subset by mapped cell type:",
                                                          choices = celltypes,
                                                          selected = setdiff(celltypes, "Low quality"),
                                                          multiple = TRUE)
                                    ),
                                    
                                    column(8,
                                           mainPanel(
                                               
                                               HTML("<p style=\"text-align: center;\"><em>Plot of imputed spatial expression, and distribution of imputed expression per cell type (below). Check option on left panel for plot of cells' segmentation.</em></p>"),
                                               
                                               girafeOutput("spatialPlot_imp",
                                                            height = "900px",
                                                            width = "1200px"
                                               ) %>% withSpinner(color="#0dc5c1"),
                                               
                                               plotOutput("spatialPlot_leg_imp",
                                                          height = "200px",
                                                          width = "1200px"),
                                               
                                               girafeOutput("violinPlot_imp",
                                                            height = "600px",
                                                            width = "1200px"
                                               ) %>% withSpinner(color="#0dc5c1")
                                           )
                                    )
                           ),
                           tabPanel("Digital in situ", fluid = TRUE,
                                    
                                    column(2,
                                           
                                           
                                           selectizeInput("gene_names_mRNA",
                                                          "Select up to five genes:",
                                                          choices = genes,
                                                          selected = c("Cdh5","Dlk1","Postn","Tbx5"),
                                                          multiple = TRUE),
                                           
                                           textInput("colours_mRNA",
                                                     "mRNA dot colours (separated by spaces)",
                                                     value = "cyan green orange red",
                                                     placeholder = "e.g. \"yellow red blue\" (without the quotes)"
                                           ),
                                           
                                           actionButton("go","Plot digital in situ", width = "100%"),
                                           
                                           HTML("<i> </i>"),
                                           
                                           downloadButton("download_mRNAPlot", "Download PDF",
                                                          style="width:100%"),
                                           
                                           checkboxInput("mRNA_full",
                                                         "Plot whole embryo (ignores regional selection)",
                                                         value = FALSE),
                                           
                                           HTML("<i>FURTHER OPTIONS</i>"),
                                           
                                           selectizeInput('celltype_outline_mRNA',
                                                          'Cell types to outline',
                                                          choices = celltypes,
                                                          selected = "Cardiomyocytes",
                                                          multiple = TRUE),
                                           
                                           # this is here so that double click in digital in situ works
                                           verbatimTextOutput("info"),
                                           
                                           sliderInput('radius',
                                                       "Digital in situ area width (only used if \"Plot whole embryo\" is unchecked)",
                                                       min = 0.75,
                                                       max = 5,
                                                       value = 1.25,
                                                       round = -1),
                                           
                                           sliderInput('embryo1_centre_x',
                                                       "embryo1 centre value x",
                                                       min = embryo_coords_range_x[1],
                                                       max = embryo_coords_range_x[2],
                                                       value = -0.6,
                                                       round = -1),
                                           
                                           sliderInput('embryo1_centre_y',
                                                       "embryo1 centre value y",
                                                       min = embryo_coords_range_y[1],
                                                       max = embryo_coords_range_y[2],
                                                       value = -0.1,
                                                       round = -1),
                                           
                                           sliderInput('embryo2_centre_x',
                                                       "embryo2 centre value in x",
                                                       min = embryo_coords_range_x[1],
                                                       max = embryo_coords_range_x[2],
                                                       value = -1.3,
                                                       round = -1),
                                           
                                           sliderInput('embryo2_centre_y',
                                                       "embryo2 centre value in y",
                                                       min = embryo_coords_range_y[1],
                                                       max = embryo_coords_range_y[2],
                                                       value = -0.25,
                                                       round = -1),
                                           
                                           sliderInput('embryo3_centre_x',
                                                       "embryo3 centre value in x",
                                                       min = embryo_coords_range_x[1],
                                                       max = embryo_coords_range_x[2],
                                                       value = 0,
                                                       round = -1),
                                           
                                           sliderInput('embryo3_centre_y',
                                                       "embryo3 centre value in y",
                                                       min = embryo_coords_range_y[1],
                                                       max = embryo_coords_range_y[2],
                                                       value = -0.6,
                                                       round = -1),
                                           
                                           downloadButton("download_mRNAPlot_gg", "Download ggplot object",
                                                          style="width:100%")
                                           
                                    ),
                                    
                                    column(8,
                                           mainPanel(
                                               
                                               HTML("<p style=\"text-align: center;\"><em>Double click on each panel to re-centre the digital in situ regional selection, before clicking \"Plot digital in situ\".</em></p>"),
                                               
                                               plotOutput("mRNARegionPlot",
                                                          dblclick = "mRNARegionPlot_dblclick",
                                                          width = "1200px"
                                               ),
                                               
                                               plotOutput("mRNAPlot",
                                                          height = "600px",
                                                          width = "1200px"
                                               )  %>% withSpinner(color="#0dc5c1")
                                           )
                                    )
                                    
                           ),
                           tabPanel("Virtual dissection", fluid = TRUE,
                                    
                                    column(2,
                                           
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           actionButton("add", "Add selection to Group A",
                                                        width = "100%",
                                                        style="height:50px;color:#ff0000"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           actionButton("add2", "Add selection to Group B",
                                                        width = "100%",
                                                        style="height:50px;color:#0000ff"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           actionButton("remove", "Remove selection from either Group",
                                                        width = "100%",
                                                        style="height:50px"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           actionButton("replotVirtualDissection", HTML("<b>Load/Reload plot</b></br>(Click to begin)"),
                                                        width = "100%",
                                                        style="height:200px;font-size:25px"),
                                           
                                           HTML("<p style=\"margin-bottom:1cm;\"> </p>"),
                                           
                                           actionButton("removeAll", HTML("<i>Reset all cells to Unselected</i>"),width = "100%",
                                                        style="height:40px"),
                                           selectizeInput('virtualDissectionCelltypes',
                                                          'Mapped cell types',
                                                          choices = celltypes,
                                                          selected = NULL,
                                                          multiple = TRUE),
                                           
                                           actionButton("addCT", "Add listed cell types to Group A",
                                                        width = "100%",
                                                        style="height:50px;color:#ff0000"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           actionButton("addCT2", "Add listed cell types to Group B",
                                                        width = "100%",
                                                        style="height:50px;color:#0000ff"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           actionButton("removeCT", "Remove listed cell types from either Group",
                                                        width = "100%",
                                                        style="height:50px;font-size:12px"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           
                                           textInput("logical",
                                                     "Advanced use: logical statement to select cells",
                                                     value = "",
                                                     placeholder = "e.g. \"Shh > 0.5\", without the quotes"),
                                           
                                           actionButton("addLogical", "Add logically selected cells to Group A",
                                                        width = "100%",
                                                        style="height:50px;color:#ff0000"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           actionButton("addLogical2", "Add logically selected cells to Group B",
                                                        width = "100%",
                                                        style="height:50px;color:#0000ff"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           actionButton("removeLogical", "Remove logically selected from either Group",
                                                        width = "100%",
                                                        style="height:50px;font-size:13px"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           
                                           
                                           HTML("Only selected cells within the selected embryos and z-slices will be downloaded."),
                                           downloadButton("virtualDissectionCellDownload", "Download Group A cell names",
                                                          style="width:100%"),
                                           downloadButton("virtualDissectionCellDownload2", "Download Group B cell names",
                                                          style="width:100%"),
                                           HTML("<p style=\"margin-bottom:1cm;\"> </p>"),
                                           HTML("You can upload multiple cell pre-selection files. Once uploaded press \"Add pre-selected\" to either Group, and then \"Load/Reload plot\"."),
                                           fileInput("preselected", "Upload pre-selected cells",
                                                     multiple = TRUE,
                                                     accept = c(".Rds")),
                                           actionButton("add_pre", "Add pre-selected to Group A",
                                                        width = "100%", 
                                                        style = "font-size:13px;color:#ff0000"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           actionButton("add_pre2", "Add pre-selected to Group B",
                                                        width = "100%",
                                                        style = "font-size:13px;color:#0000ff"),
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           actionButton("remove_pre", "Remove pre-selected cells from either Group",
                                                        width = "100%",
                                                        style = "font-size:12px")
                                    ),
                                    
                                    # Show a plot of the generated distribution
                                    column(8,
                                           mainPanel(
                                               HTML("<p style=\"text-align: center;\"><em>Please press the \"Reload plot\" button to start and to update after cell additions or removals. Use the lasso selection icon (leftmost in the topright corner) to select cells, then click each button to assign them to a group. Once you have selected two groups, an MA-plot will appear below.</em></p>"),
                                               
                                               actionButton("virtualDissection_reset", label = "Clear lasso selection"),
                                               
                                               girafeOutput("virtualDissection",
                                                            height = "900px",
                                                            width = "1200px"
                                               ) %>% withSpinner(color="#0dc5c1"),
                                               
                                               splitLayout(cellWidths = c("600px", "600px"),
                                                           plotOutput("VirtualDissectionBarPlot",
                                                                      height = "900px") %>% withSpinner(color="#0dc5c1"),
                                                           girafeOutput("VirtualDissectionMAPlot", height = "900px") %>% withSpinner(color="#0dc5c1")
                                               )
                                               
                                           )
                                    )
                           )
                           
               )
        )
    )
)
)