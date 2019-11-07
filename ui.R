dash.sidebar <- dashboardSidebar(width = 320,
                                 # HTML tag for margin
                                 tags$style(HTML("#sidebar_margin {margin: 10px}")),
                                 
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
                                             fluidRow(id = "sidebar_margin", downloadButton("downloadData", "Download Significantly Changed Metabolites"))
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

ui <- dashboardPage(
    dashboardHeader(titleWidth = 320, title="MS-MZ-Met-Metabo-Metabalyzer"),
    dash.sidebar,
    dash.body
)
