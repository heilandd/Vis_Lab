##----------------------------------------------------------------------------##
## Custom functions.
##----------------------------------------------------------------------------##
APPBox <- function(title, content) {
  box(
    title = title,
    status = "primary",
    solidHeader = TRUE,
    width = 12,
    collapsible = TRUE,
    content
  )
}
APPInfoButton <- function(id) {
  actionButton(
    inputId = id,
    label = "info",
    icon = NULL,
    class = "btn-xs",
    title = "Show additional information for this panel."
  )
}
boxTitle <- function(title) {
  p(title, style = "padding-right: 5px; display: inline")
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##

source(paste0(folder,"/Test/ui.R"),local = T)
source(paste0(folder,"/DE/ui.R"),local = T)
source(paste0(folder,"/GSEA/ui.R"),local = T)

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
ui <- dashboardPage(
  dashboardHeader(
    title = span("VIS_Lab ", style = "color: white; font-size: 28px; font-weight: bold")
  ),
  dashboardSidebar(
    tags$head(tags$style(HTML(".content-wrapper {overflow-x: scroll;}"))),
    sidebarMenu(
      sidebarMenuOutput("sidebar_menu")
    )
  ),
  dashboardBody(
    tags$script(HTML('$("body").addClass("fixed");')),
    tabItems(
      tab_load_data,
      tab_DEAnalysis,
      tab_GSEA
    )
  )
)

