
##----------------------------------------------------------------------------##
## UI-GSEA
##----------------------------------------------------------------------------##

tab_GSEA <- tabItem(
  tabName = "GSEA",

  fluidRow(

    #Selection Box
    box(
      width = 4, status = "warning", solidHeader = TRUE,
      title = "Select conditions for GSEA:",
      
    #### Checkbox for Groups
    uiOutput("Choice_GSEA"),
    uiOutput("Condition_GSEA"),
    uiOutput("PW"),

    ),
    
    
    #Plot TabBox
    tabBox(
      title = "", 
      id="Plots_GSEA_1",
      width = 8,
      
      tabPanel(title="GSEA-Plot",
              downloadLink("downloadplot_GSEA_Plot", "PDF"),
              plotOutput("GSEA_Plot")),
      
      tabPanel(title="Sorted GSEA: ",
               downloadLink("downloadSort_GSEA_Plot", "PDF"),
               uiOutput("PW_X1"),
               
               plotOutput("Sort_GSEA_Plot"))
      
    ),
    

    
    #Put Table Here
    box(
      width = 12, status = "warning", solidHeader = TRUE,
      title = "GSEA Table:",
      DT::DTOutput("GSEA_out_tab", width = "100%")
      
      
    ),
    
    ),
  )


