
##----------------------------------------------------------------------------##
## UI-Load Data
##----------------------------------------------------------------------------##

tab_load_data <- tabItem(
  tabName = "loadData",
  fluidRow(
    column(12,
           titlePanel("Load data"),
           fileInput(
             inputId = "input_file",
             label = "Select input data (.milo_rds)...",
             multiple = FALSE,
             accept = c(".milo_rds"),
             width = '350px',
             buttonLabel = "Browse...",
             placeholder = "No file selected"
           )
    )
  ),
  fluidRow(
    valueBoxOutput("load_data_number_of_cells"),
    valueBoxOutput("load_data_number_of_Groups"),
    valueBoxOutput("load_data_number_of_Batches")
  )
)
