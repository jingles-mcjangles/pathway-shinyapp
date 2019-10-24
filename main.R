library("pacman")
pacman::p_load("shiny", "shinydashboard", "shinycssloaders", 
               "ggplot2", 
               "MetaboAnalystR", "KEGGREST", "tidyverse", "ggfortify",
               "httr", "XML", "DT")

source("helper-func.R")
source("ui.R", local = TRUE)
source("server.R", local = TRUE)

# Run the app
shinyApp(ui = ui, server= server)
