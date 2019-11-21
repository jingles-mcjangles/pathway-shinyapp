library("pacman")
pacman::p_load(
  "shiny", "shinydashboard", "shinycssloaders", "ggplot2",
  "MetaboAnalystR", "KEGGREST", "tidyverse", "ggfortify",
  "httr", "XML", "DT"
)

# set working directory
# setwd('./pathway-shinyapp/')

# source functions, ui.R and server.R files
source("helper-func.R")
source("ui.R")
source("server.R")

runApp()
