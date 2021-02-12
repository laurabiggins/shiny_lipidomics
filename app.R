#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

ui <- fluidPage(
    tags$script("
                var myWidth = 0;
                $(document).on('shiny:connected', function(event) {
                  myWidth = $(window).width();
                  Shiny.onInputChange('shiny_width', myWidth);
                });
                $(window).resize(function(event) {
                   myWidth = $(window).width();
                   Shiny.onInputChange('shiny_width', myWidth);
                });
              "),
    tags$script("
                var myHeight = 0;
                $(document).on('shiny:connected', function(event) {
                  myHeight = $(window).height();
                  Shiny.onInputChange('shiny_height', myHeight);
                });
                $(window).resize(function(event) {
                   myHeight = $(window).height();
                   Shiny.onInputChange('shiny_height', myHeight);
                });
                "),
  titlePanel("Heatmap"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload file"),
      textInput("group_types", "group names", value = "WT, KO"),
      br(),
      textOutput("lipid_column"),
      br(),
      actionButton("plot_heatmap" , "Plot heatmap"),
      numericInput("plot_height", "plot height", 600),
      br(),
      br(),
      downloadButton("download_png", "png"),
      downloadButton("download_pdf", "pdf")
    ),
    mainPanel(
        tableOutput("head"),
        plotOutput("heatmap", height = "600px"),
        actionButton("browser", "browser")
    )
  )
)


server <- function(input, output) {
    
    rv <- reactiveValues()
    
    groups <- reactive({
        
        x <- stringr::str_split(input$group_types, pattern = ",", simplify = TRUE)
        stringr::str_trim(x)  
    })
    
    lipid_column <- reactive(colnames(raw_dataset())[1])
    
    raw_dataset <- reactive({
        
        req(input$file)
        ext <- tools::file_ext(input$file$name)
        switch(ext,
               csv = vroom::vroom(input$file$datapath, delim = ","),
               tsv = vroom::vroom(input$file$datapath, delim = "\t"),
               validate("Invalid file; Please upload a .csv or .tsv file")
        )
    })
    
    heatmap_data <- reactive({
        raw_dataset <- remove_0_variance(raw_dataset())
        tibble::column_to_rownames(raw_dataset(), lipid_column())
    })
    
    output$head <- renderTable({
        head(raw_dataset(), n = 3)
    })
    
    output$lipid_column <- renderText({
        paste0("Column containing lipid names is assumed to be ", lipid_column())
    })
    
    observeEvent(input$plot_heatmap, {
        
        df_groups <- sapply(colnames(heatmap_data()), get_group, groups = groups())

        # The lipid classes
        df_lipids <- parse_lipid_classes(raw_dataset(), lipid_column())

        # pheatmap requires a dataframe with sample names/lipid names as rownames
        rv$annot_col <- tibble::column_to_rownames(tibble::enframe(df_groups), "name")
        rv$annot_row <- tibble::column_to_rownames(df_lipids, "names")

        #===============================================================================
        # Lipid classes need further sorting to create custom row annotations where only
        # the class name is displayed rather than each individual lipid.

        lipid_summary <- get_lipid_summary(df_lipids)
        rv$lipid_labels <- create_lipid_class_labels(lipid_summary, nrow(df_lipids))

        lipid_colours <- create_lipid_colours(class_names = dplyr::pull(lipid_summary, class))
        group_colours <- stats::setNames(c("blue", "green"), groups())

        rv$colours <- list(
            value = group_colours,
            class = lipid_colours
        )
    })
    
    heatmap_obj <- eventReactive(input$plot_heatmap, {
        pheatmap::pheatmap(
            heatmap_data(),
            scale = "row",
            annotation_row = rv$annot_row,
            annotation_col = rv$annot_col,
            cluster_rows = FALSE,
            labels_row = rv$lipid_labels,
            fontsize_row = 8,
            annotation_colors = rv$colours,
            annotation_legend = FALSE,
            silent = TRUE
        )
    })    
   
    output$heatmap <- renderPlot({
        req(input$plot_height)
        plot(heatmap_obj()$gtable)
    }, 
    height = function(x) input$plot_height)
    
    
    output$download_png <- downloadHandler(
        filename = function() {
            paste0("heatmap.png")
        },
        content = function(file) {
            ggplot2::ggsave(
                file, 
                heatmap_obj(), 
                device = "png",
                width = input$shiny_width/4, ## 1pixel ~ 0.26mm at 96 dpi. it's ~0.35 at 72dpi
                height = input$shiny_height/4, 
                units = "mm"
            )
        }
    )
    

    observeEvent(input$browser, browser())
}

# Run the application 
shinyApp(ui = ui, server = server)



# #======================================================================================
# # sorting out the sample names and groups for colouring the row and column annotations
# 
# # The sample types eg KO/WT.

# 



