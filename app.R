library("pacman")
pacman::p_load("shiny", "shinydashboard", "shinycssloaders", 
               "ggplot2", 
               "MetaboAnalystR", "KEGGREST", "tidyverse", "ggfortify",
               "httr", "XML", "DT")

source("/Users/don/Documents/pathway-shinyapp/helper-func.R")

#https://stackoverflow.com/questions/40152857/how-to-dynamically-populate-dropdown-box-choices-in-shiny-dashboard
#https://stackoverflow.com/questions/34555548/adding-a-download-export-button-to-r-shinydashboard

dash.sidebar <- dashboardSidebar(
    sidebarMenu(id="sidebarmenu", startExpanded=T,
        menuItem("Home", tabName="main", icon = icon("home")),
        fileInput("file1", "Upload CSV File",
                  accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
        sliderInput("fdr_num", "False Discovery Rate", min=0, max=1.0, step=0.01, value=0.3, width = '95%'),
        menuItem("Advanced Options", tabName="options", 
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
                    width=12,
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "tabset_plots2", height = "580px",
                    tabPanel("iPath: Coverage", 
                             div(style = 'height:510px; width:800px; overflow-y: scroll; overflow-x: scroll; font-size:65%', 
                                 htmlOutput("ipath") %>% withSpinner())),
                    tabPanel("Significantly Changing Metabolites", 
                             div(style = 'height:510px; width:800px; overflow-y: scroll; overflow-x: scroll; font-size:65%', 
                                 dataTableOutput("de_tbl") %>% withSpinner()))
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
                    tabPanel("PCA", 
                             div(style = 'height:330px; width:500px; overflow-y: scroll; overflow-x: scroll; font-size:65%', 
                                 plotOutput("pca_plot") %>% withSpinner())),
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
    dashboardHeader(title="MS-MZ-Met-Metabo-Metabalyzer"),
    dash.sidebar,
    dash.body
)


# Define server logic
server <- function(input, output, session) {

    # Read user-uploaded csv as DT
    output$input_table <- renderDataTable({
        req(input$file1)
        DT::datatable(read_csv(input$file1$datapath), options = list(ordering = F, searching = T, pageLength=20))
    })
    
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
    
    # Show PCA
    output$pca_plot <- renderPlot({
        req(input$file1)
        tbl0 <- read_csv(input$file1$datapath)
        #DT::datatable(read_csv(input$file1$datapath), options = list(ordering = F, searching = T, pageLength=20))
        
        autoplot(prcomp(tbl0[,3:ncol(tbl0)]), data=tbl0, colour = "Group", label=T, label.size=3)
    })
    
    # DE metabolites
    output$de_tbl <- renderDataTable({
        
        req(input$file1)
        # Load data, specify some inputs
        tbl0 <- read_csv(input$file1$datapath)
        
        group.name.col <- names(tbl0[,2])
        groupnames <- as.vector(unlist(unique(tbl0[[group.name.col]])))
        control.group.name <- "Plus" # name of denominator group. 
        alpha <- input$fdr_num
        #kegg_species_id <- "dme"
        
        #  ==================== START ==================== 
        cnames <- colnames(tbl0)[3:length(colnames(tbl0))]
        # Init vec of t.stats p-vals, and fold change of averages
        # Compute t-stats
        p.vals.ls <- vector(mode="numeric", length = length(cnames))
        fc.ls <- vector(mode="numeric", length = length(cnames))
        for (i in 1:length(cnames)) {
            g1.vec <- as.vector(unlist(tbl0 %>% filter(!!sym(group.name.col)==groupnames[1]) %>% select(cnames[i])))
            g2.vec <- as.vector(unlist(tbl0 %>% filter(!!sym(group.name.col)==groupnames[2]) %>% select(cnames[i])))
            x <- t.test(g1.vec, g2.vec)
            p.vals.ls[i] <- x$p.value
            
            mu1 <- mean(g1.vec)
            mu2 <- mean(g2.vec)
            fc.ls[i] <- mu1/mu2
        }
        names(p.vals.ls) <- colnames(tbl0)[3:length(colnames(tbl0))]
        names(fc.ls) <- colnames(tbl0)[3:length(colnames(tbl0))]
        
        # adjust: BH correction
        # pick out names where p < alpha
        p.vals.ls <- p.adjust(p.vals.ls, method = "hochberg", n = length(p.vals.ls))
        significant.metabs <- c()
        for (nm in names(p.vals.ls)) {
            if (p.vals.ls[[nm]] < alpha) {
                significant.metabs <- c(significant.metabs, nm)
            }
        }
        print(paste0(length(significant.metabs), " significant metabolites found at alpha = ", alpha))
        
        # Filter for significant metabs
        p.vals.ls2 <- p.vals.ls[significant.metabs]
        fc.ls2 <- fc.ls[significant.metabs]
        
        # Colour by FC
        th.lower <- 0.8
        th.upper <- 1.2
        fc.colour.ls <- rep("#000000", length(significant.metabs))
        names(fc.colour.ls) <- significant.metabs
        for (i in 1:length(significant.metabs)) {
            nm <- significant.metabs[i]
            if (fc.ls[nm] <= th.lower) {
                fc.colour.ls[i] <- "red" #crimson
            }
            if (fc.ls[nm] > th.upper) {
                fc.colour.ls[i] <- "green" #forest green
            }
        }
        fc.tbl <- tibble::enframe(fc.ls2) %>% rename("fc"=value, "Sample"=name)
        p.val.tbl <- tibble::enframe(p.vals.ls2) %>% rename("FDR"=value, "Sample"=name)
        
        #  ==================== GET CPD IDS ==================== 
        kegg.id.tbl <- lookup_chem_id(names(fc.ls2))
        kegg.id.vec <- as.vector(unlist(kegg.id.tbl %>% select("KEGG")))
        
        tbl1 <- inner_join(fc.tbl, kegg.id.tbl, by="Sample")
        tbl1 <- inner_join(tbl1, p.val.tbl, by="Sample") %>% select(Sample, fc, KEGG, FDR)
        
        # Display
        DT::datatable(tbl1, options = list(ordering = F, searching = T, pageLength=20))
    })
    
    # PE Analysis
    output$ipath <- renderText({
        req(input$file1)
        tbl0 <- read_csv(input$file1$datapath)
        
        group.name.col <- names(tbl0[,2])
        groupnames <- as.vector(unlist(unique(tbl0[[group.name.col]])))
        control.group.name <- "Plus" # name of denominator group. 
        alpha <- input$fdr_num
        kegg_species_id <- "dme"
        
        #  ==================== START ==================== 
        cnames <- colnames(tbl0)[3:length(colnames(tbl0))]
        # Init vec of t.stats p-vals, and fold change of averages
        # Compute t-stats
        p.vals.ls <- vector(mode="numeric", length = length(cnames))
        fc.ls <- vector(mode="numeric", length = length(cnames))
        for (i in 1:length(cnames)) {
            g1.vec <- as.vector(unlist(tbl0 %>% filter(!!sym(group.name.col)==groupnames[1]) %>% select(cnames[i])))
            g2.vec <- as.vector(unlist(tbl0 %>% filter(!!sym(group.name.col)==groupnames[2]) %>% select(cnames[i])))
            x <- t.test(g1.vec, g2.vec)
            p.vals.ls[i] <- x$p.value
            
            mu1 <- mean(g1.vec)
            mu2 <- mean(g2.vec)
            fc.ls[i] <- mu1/mu2
        }
        names(p.vals.ls) <- colnames(tbl0)[3:length(colnames(tbl0))]
        names(fc.ls) <- colnames(tbl0)[3:length(colnames(tbl0))]
        
        # adjust: BH correction
        # pick out names where p < alpha
        p.vals.ls <- p.adjust(p.vals.ls, method = "hochberg", n = length(p.vals.ls))
        significant.metabs <- c()
        for (nm in names(p.vals.ls)) {
            if (p.vals.ls[[nm]] < alpha) {
                significant.metabs <- c(significant.metabs, nm)
            }
        }
        print(paste0(length(significant.metabs), " significant metabolites found at alpha = ", alpha))
        
        # Filter for significant metabs
        p.vals.ls2 <- p.vals.ls[significant.metabs]
        fc.ls2 <- fc.ls[significant.metabs]
        
        # Colour by FC
        th.lower <- 0.8
        th.upper <- 1.2
        fc.colour.ls <- rep("#000000", length(significant.metabs))
        names(fc.colour.ls) <- significant.metabs
        for (i in 1:length(significant.metabs)) {
            nm <- significant.metabs[i]
            if (fc.ls[nm] <= th.lower) {
                fc.colour.ls[i] <- "#e41a1c" #red
            }
            if (fc.ls[nm] > th.upper) {
                fc.colour.ls[i] <- "#4daf4a" #green
            }
        }
        print("check")
        fc.tbl <- tibble::enframe(fc.ls2) %>% rename("fc"=value, "Sample"=name)
        fc.colour.tbl <- tibble::enframe(fc.colour.ls) %>% rename("fc_colour"=value, "Sample"=name)
        fc.tbl <- inner_join(fc.tbl, fc.colour.tbl, by="Sample")
        p.val.tbl <- tibble::enframe(p.vals.ls2) %>% rename("FDR"=value, "Sample"=name)
        
        #  ==================== GET CPD IDS ==================== 
        kegg.id.tbl <- lookup_chem_id(names(fc.ls2))
        kegg.id.vec <- as.vector(unlist(kegg.id.tbl %>% select("KEGG")))
        
        tbl1 <- inner_join(fc.tbl, kegg.id.tbl, by="Sample")
        tbl1 <- inner_join(tbl1, p.val.tbl, by="Sample") %>% select(Sample, fc, fc_colour, KEGG, FDR)
        
        #  ======================================== POST Calls to iPath API ========================================
        # Make selection str
        kegg.vec <- as.vector(unlist(kegg.id.tbl["KEGG"]))
        names(kegg.vec) <- as.vector(unlist(kegg.id.tbl["Sample"]))
        selection.str.ls <- c()
        for (nm in names(kegg.vec)) {
            selection.str.ls <- c(selection.str.ls, paste0(kegg.vec[nm], " W20 ", fc.colour.ls[nm]))
        }
        selection.str <- paste0(selection.str.ls, collapse = "\n")
        
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
                     default_radius="5",
                     whole_modules="0")
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
