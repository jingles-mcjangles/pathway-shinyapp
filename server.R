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