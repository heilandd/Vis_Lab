Vis_Lab <- function(folder=c("/VIS_Lab")){
  
  library(shinydashboard)
  library(shiny)
  ################## Sources of shiny parts of the App

  ##--------------------------------------------------------------------------##
  ## Load server and UI functions.
  ##--------------------------------------------------------------------------##
  
  
  setwd(folder)
  source(paste0(folder,"/App/ALL_UI.R"))
  source(paste0(folder,"/App/ALL_server.R"))
  ##--------------------------------------------------------------------------##
  ## Launch app.
  ##--------------------------------------------------------------------------##
  
  shinyApp(ui = ui, server = server)
  
}






















