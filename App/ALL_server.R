##----------------------------------------------------------------------------##
##  Server Functions 
##----------------------------------------------------------------------------##

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=1000*1024^2)
  
  ##--------------------------------------------------------------------------##
  ## Colors
  ##--------------------------------------------------------------------------##
  
  
  ##--------------------------------------------------------------------------##
  ## Central parameters.
  ##--------------------------------------------------------------------------##
  scatter_plot_dot_size <- list(
    min = 1,
    max = 20,
    step = 1,
    default = 5
  )
  
  scatter_plot_alpha <- list(
    min = 0.1,
    max = 1.0,
    step = 0.1,
    default = 1.0
  )
  
  preferences <- reactiveValues(use_webgl = TRUE)
  
  ##--------------------------------------------------------------------------##
  ## Sidebar menu.
  ##--------------------------------------------------------------------------##
  output[["sidebar_menu"]] <- renderMenu({
    sidebarMenu(id = "sidebar",
                menuItem(
                  "Load data", tabName = "loadData",
                  icon = icon("spinner"), selected = TRUE
                ),
                menuItem(
                  "DE-Analysis", tabName = "DEAnalysis",
                  icon = icon("brain"), selected = F
                ),
                menuItem(
                  "GSEA", tabName = "GSEA",
                  icon = icon("chart-line"), selected = F
                )
    )
  })

  
  ##--------------------------------------------------------------------------##
  ## Sample data
  ##--------------------------------------------------------------------------##
  sample_data <- reactive({
    if ( is.null(input[["input_file"]]) || is.na(input[["input_file"]]) ) {
      sample_data <- readRDS("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/VIS_Lab/sample_data.RDS")
      print("We use non selected DF")
      return(sample_data)
    } else {
      req(input[["input_file"]])
      sample_data <- readRDS(input[["input_file"]]$datapath)
      print(sample_data@fdata)
      print(sample_data@data@counts[1:4,1:4])
      return(sample_data)
    }
  })
  
  
  ##--------------------------------------------------------------------------##
  ## Tabs.
  ##--------------------------------------------------------------------------##
  source(paste0(folder,"/Test/server.R"), local = T)
  source(paste0(folder,"/DE/server.R"), local = T)
  source(paste0(folder,"/GSEA/server.R"), local = T)
}

