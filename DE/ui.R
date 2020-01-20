
##----------------------------------------------------------------------------##
## UI-DE
##----------------------------------------------------------------------------##

tab_DEAnalysis <- tabItem(
  tabName = "DEAnalysis",
  fluidRow(valueBoxOutput("DE"),
           valueBoxOutput("FDR"),
           valueBoxOutput("TotalGenes")),
  fluidRow(
    
    
    
    
    #Selection Box
    box(
      width = 3, status = "warning", solidHeader = TRUE,
      title = "Select conditions for Differential expressed Genes:",
      
    #### Checkbox for Groups
    uiOutput("Choice"),
    
    "Please select two conditions to Analyze:",
    
    uiOutput("C1"),
    
    uiOutput("C2"),
                
    #### Select Gene
    
    uiOutput("Gene"),
    
    
    ),
    
    
    #Plot TabBox
    tabBox(
      title = "", 
      id="Plots2",
      width = 8,
      
      tabPanel(title="Volcano Plot",
              downloadLink("downloadVolcano_plot", "PDF"),
              plotOutput("Volcano_plot")),
      
      tabPanel(title="Heatmap",
               downloadLink("downloadHeatmap", "PDF"),
               uiOutput("Number_of_genes"),
               plotOutput("Heatmap"))
      
    ),
    
    tabBox(
      title = "Plots", 
      id="Plots1",
      width = 6,
      tabPanel(title="Violine Plot",
               downloadLink("downloadplot_violin", "PDF"),
               plotOutput("plot_violin")),
      tabPanel(title="Barplot",
               downloadLink("downloadplot_bar", "PDF"),
               plotOutput("Barplot"))
      
    ),
    
    #Put Table Here
    box(
      width = 12, status = "warning", solidHeader = TRUE,
      title = "Table of Differential Expressed Genes:",
      DT::DTOutput("DE_Tab")
      
      
    ),
    
    ),
  )

