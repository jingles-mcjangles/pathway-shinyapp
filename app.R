library("pacman")
pacman::p_load("shiny", "shinydashboard", "shinycssloaders", 
               "ggplot2", 
               "MetaboAnalystR", "KEGGREST", "tidyverse", "ggfortify",
               "httr", "XML", "DT")

source("/Users/don/Documents/pathway-shinyapp/helper-func.R")

#https://stackoverflow.com/questions/40152857/how-to-dynamically-populate-dropdown-box-choices-in-shiny-dashboard
#https://stackoverflow.com/questions/34555548/adding-a-download-export-button-to-r-shinydashboard

dash.sidebar <- dashboardSidebar(width = 320,
    sidebarMenu(id="sidebarmenu", startExpanded=T,
        menuItem("Home", tabName="main", icon = icon("home")),
        fileInput("file1", "Upload CSV File",
                  accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
        sliderInput("fdr_num", "False Discovery Rate", min=0, max=1.0, step=0.01, value=0.05, width = '95%'),
        selectizeInput('grp_numerator', 'Select group1 (FC numerator)', ""),
        selectizeInput('grp_denominator', 'Select group2 (FC denominator)', ""),
        textInput("species_code", "KEGG 3-letter Species Code"),
        
        menuItem("Advanced Options", tabName="options", 
                 sliderInput("alpha", "Alpha for t-tests", min=0, max=1.0, step=0.01, value=0.01, width = '95%'),
                 sliderInput("fc", "Fold Change Threshold", min=0, max=1.0, step=0.05, value=0.2, width = '95%'),
                 checkboxInput("median_normalization","Use median normalization",TRUE),
                 checkboxInput("logtransform","Use log2 transform",TRUE),
                 checkboxInput("autoscale","Use Autoscaling",TRUE)
                 ),
        
        tags$hr(),
        menuItem("Help", tabName="help", icon = icon("info-circle")),
        menuItem("About", tabName="about", icon = icon("info-circle")),
        tags$hr(),
        downloadButton("downloadData", "Download Analysis Files")
        )
)

dash.body <- dashboardBody(
    tabItems(
        tabItem(tabName = "main",
            # ===== Analysis fluidRow =====
            fluidRow(
                tabBox(
                    width=8,
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "tabset_plots2", height = "580px",
                    tabPanel("iPath: Coverage", 
                             div(style = 'height:510px; overflow-y: scroll; overflow-x: scroll; font-size:75%', 
                                 htmlOutput("ipath_coverage") %>% withSpinner())),
                    tabPanel("iPath: Significant Metabolites", 
                             div(style = 'height:510px; overflow-y: scroll; overflow-x: scroll; font-size:75%', 
                                 htmlOutput("ipath_de_metabs") %>% withSpinner())),
                    tabPanel("Significantly Changing Metabolites", 
                             div(style = 'height:510px; overflow-y: scroll; overflow-x: scroll; font-size:75%', 
                                 dataTableOutput("de_tbl") %>% withSpinner()))
                    ),
                box(width=4,
                    title="Boxplot", height = "400px",
                    selectizeInput('selected_colname', 'Select Metabolite', ""),
                    div(style = 'height:330px; overflow-y: scroll; overflow-x: scroll', 
                        plotOutput("single_boxplot") %>% withSpinner())
                    )
                ),
            # ===== About Your Data fluidRow =====
            fluidRow(
                tabBox(
                    width=12,
                    title = "About Your Data",
                    id = "tabset_plots", height = "400px",
                    tabPanel("Input Data", 
                             div(style = 'height:330px; width:500px; overflow-y: scroll; overflow-x: scroll; font-size:65%', 
                                 dataTableOutput("input_table") %>% withSpinner())),
                    #tabPanel("PCA", 
                    #         div(style = 'height:330px; width:500px; overflow-y: scroll; overflow-x: scroll; font-size:65%', 
                    #             plotOutput("pca_plot") %>% withSpinner())),
                    tabPanel("Unidentified Compounds", 
                             div(style = 'height:330px; width:500px; overflow-y: scroll; overflow-x: scroll; font-size:65%', 
                                 dataTableOutput("unfound_table") %>% withSpinner()))
                    )
                )
        ),
        tabItem(tabName = "help", includeHTML("help.html")),
        tabItem(tabName = "about", includeHTML("about.html"))
    )
)

# Define UI
ui <- dashboardPage(
    dashboardHeader(titleWidth = 320, title="MS-MZ-Met-Metabo-Metabalyzer"),
    dash.sidebar,
    dash.body
)


# Define server logic
server <- function(input, output, session) {
    # ========== DEF REACTIVES =========
    # Reactive for DE metabs and ipath
    data_tbl <- reactive({
        req(input$file1)
        req(input$grp_numerator)
        req(input$grp_denominator)
        tbl0 <- read_csv(input$file1$datapath)
        #tbl0 <- clean_and_transform_tibble(tbl0, log_bool=TRUE, z_transform_bool=FALSE)
        tbl0 <- get_results_tibble(tbl0, 
                                   input_alpha=0.01, 
                                   grp_numerator = input$grp_numerator,
                                   grp_denominator = input$grp_denominator,
                                   input$fdr_num)
        tbl0
    })
    
    grp.name.vec <- reactive({
        req(input$file1)
        tbl0 <- read_csv(input$file1$datapath)
        as.vector(unlist(tbl0[2]))
    })
    
    outVar <- reactive({
        req(input$file1)
        tbl0 <- read_csv(input$file1$datapath)
        sort(colnames(tbl0[3:ncol(tbl0)]))
    })
    
    
    # ========== DEF OBSERVES =========
    # To update various selectizeInputs
    observe({
        updateSelectizeInput(session, "selected_colname",
                             choices = outVar())
    })
    observe({
        updateSelectizeInput(session, "grp_numerator",
                             choices = grp.name.vec(), 
                             selected = grp.name.vec()[1])
    })
    observe({
        updateSelectizeInput(session, "grp_denominator",
                             choices = grp.name.vec(), 
                             selected = grp.name.vec()[2])
    })
    
    # ========== DEF OUTPUTS =========
    # Individual metab box plots: 
    # dynamically update selectizeInput with sorted metab names
    output$single_boxplot <- renderPlot({
        req(input$file1)
        tbl0 <- read_csv(input$file1$datapath)
        colnames.to.keep <- c(colnames(tbl0[1:2]), input$selected_colname)
        tbl_selected <- tbl0 %>% select(colnames.to.keep)
        group.colname <- colnames(tbl0[2])
        ggplot(tbl_selected, aes(x=!!sym(group.colname), y=!!sym(input$selected_colname))) + 
            geom_boxplot() + 
            geom_dotplot(binaxis='y')
    })
    # Dynamically update selectizeInputs with group names
    
    
    # Display untransformed user input
    output$input_table <- renderDataTable({
        req(input$file1)
        DT::datatable(read_csv(input$file1$datapath), options = list(ordering = F, searching = T, pageLength=20))
    })
    
    # KEGGless compounds
    output$unfound_table <- renderDataTable({
        req(input$file1)
        tbl0 <- read_csv(input$file1$datapath)
        
        metab.names.ls <- colnames(tbl0)[3:ncol(tbl0)]
        mSet<-InitDataObjects("NA", "utils", FALSE)
        mSet<-Setup.MapData(mSet, metab.names.ls);
        mSet<-CrossReferencing(mSet, "name", hmdb=F, pubchem=F, chebi=F, kegg=T, metlin=F);
        mSet<-CreateMappingResultTable(mSet)
        tmp <- as_tibble(mSet$dataSet$map.table) %>% mutate_each(funs(replace(., .=='', "NA")))
        unfound <- tmp %>% filter(KEGG=="NA") %>% select(Query)
        DT::datatable(unfound)
    })
    
    # PCA plot
    #output$pca_plot <- renderPlot({
    #    req(input$file1)
    #    tbl0 <- read_csv(input$file1$datapath)
    #    autoplot(prcomp(tbl0[,3:ncol(tbl0)]), 
    #             data=tbl0, 
    #             colour = !!sym(colnames(tbl0[2])), 
    #             label=F)
    #})
    
    # DE metabolites
    output$de_tbl <- renderDataTable({
        tbl1 <- data_tbl() %>% filter(adj_p_val < input$fdr_num)
        DT::datatable(tbl1, options = list(ordering = T, searching = T, pageLength=20))
    })
    
    # ipath
    output$ipath_coverage <- renderText({
        req(input$file1)
        req(input$grp_numerator)
        req(input$grp_denominator)
        req(input$species_code)
        
        tbl1 <- data_tbl() %>% filter(KEGG != "undef")
        selection.str <- get_ipath_selection_str(tbl1, "Sample", "KEGG", "fc_colour", "W20")
        
        ## PARAMS
        highlight_path_width <- 10
        highlight_path_opacity <- 0.3
        module_ellipse_radius <- 5
        
        #palette <- brewer.pal(length(significant.pws), name="Set3")
        url <- "https://pathways.embl.de/mapping.cgi"
        
        # POST request to ipath where `whole_modules`=1
        # This is the main file being modified
        body <- list(selection = selection.str,
                     export_type="svg", 
                     default_opacity="0.7",
                     default_width="1",
                     default_radius="2",
                     whole_modules="0", 
                     tax_filter=input$species_code)
        #print("Making POST request where whole_module == 1...")
        r <- POST(url, body = body, encode = "form")
        print(http_status(r)$message)
        xml_main <- content(r, "text")
        #cat("check\n")
        #outfile <- strsplit(xml_main, "\n")[[1]]
        x0_text <- strsplit(xml_main, "\n")[[1]]
        # Edit SVG total height and width
        x0_text[3] <- gsub('2250', '450', x0_text[3])
        x0_text[3] <- gsub('3774', '754.8', x0_text[3])
        xml_svg_string <- paste0(x0_text, collapse="\n")        
        c(xml_svg_string)
    })
    
    output$ipath_de_metabs <- renderText({
        req(input$file1)
        req(input$grp_numerator)
        req(input$grp_denominator)
        req(input$species_code)
        
        tbl0 <- read_csv(input$file1$datapath)
        tbl2 <- get_results_tibble(tbl0, 
                                   input_alpha=0.01, 
                                   grp_numerator = input$grp_numerator,
                                   grp_denominator = input$grp_denominator,
                                   input$fdr_num)
        
        tbl2 <- tbl2 %>% filter(KEGG != "undef") %>% filter(adj_p_val < input$fdr_num)
        selection.str2 <- get_ipath_selection_str(tbl2, "Sample", "KEGG", "fc_colour", "W20")
        
        ## PARAMS
        highlight_path_width <- 10
        highlight_path_opacity <- 0.3
        module_ellipse_radius <- 5
        
        url <- "https://pathways.embl.de/mapping.cgi"
        
        # POST request to ipath where `whole_modules`=1
        # This is the main file being modified
        body <- list(selection = selection.str2,
                     export_type="svg", 
                     default_opacity="0.7",
                     default_width="1",
                     default_radius="2",
                     whole_modules="0", 
                     tax_filter=input$species_code)
        #print("Making POST request where whole_module == 1...")
        r <- POST(url, body = body, encode = "form")
        print(http_status(r)$message)
        xml_main <- content(r, "text")
        #cat("check\n")
        #outfile <- strsplit(xml_main, "\n")[[1]]
        x0_text <- strsplit(xml_main, "\n")[[1]]
        # Edit SVG total height and width
        x0_text[3] <- gsub('2250', '450', x0_text[3])
        x0_text[3] <- gsub('3774', '754.8', x0_text[3])
        xml_svg_string <- paste0(x0_text, collapse="\n")        
        c(xml_svg_string)
        
    })
}

# Run the app
shinyApp(ui, server)

