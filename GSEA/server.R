##----------------------------------------------------------------------------##
## Server-GSEA
##----------------------------------------------------------------------------##


##----------------------------------------------------------------------------##
## uiOutput
##----------------------------------------------------------------------------##

output$Choice_GSEA=renderUI({
  selectInput(inputId="Choice_GSEA", 
              label="Select Condition Factor:",
              choices = c(sample_data()@condition_Analysis),
              selected= 1
  )
})
output$Condition_GSEA=renderUI({
  selectInput(inputId="Condition_GSEA", 
              label="Condition Factor",
              choices = as.character(na.omit(c(sample_data()@GSEA_sample_list)))
  )
})
output$PW=renderUI({
  selectInput(inputId="PW", 
              label="Pathway",
              choices = rownames(sample_data()@GSEA[[1]][[1]]@result)
  )
})
output$PW_X1=renderUI({
  selectInput(inputId="PW_X1", 
              label="Pathway1",
              choices = rownames(sample_data()@GSEA[[1]][[1]]@result),
              selected= rownames(sample_data()@GSEA[[1]][[1]]@result)[1],
              multiple = TRUE
  )
})






##----------------------------------------------------------------------------##
## Functions for Tab
##----------------------------------------------------------------------------##
require(dplyr)
require(DESeq2)

gseaplot_DHH=function(object,geneSetID,main="",add_line=T, color="black", cex=2, Add=F){
  
  require(enrichplot)
  gseaScores <- getFromNamespace("gseaScores", "DOSE")
  gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList
    
    if (is.numeric(geneSetID))
      geneSetID <- object@result[geneSetID, "ID"]
    
    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin=0
    df$ymax=0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
  }
  gsdata <- gsInfo(object, geneSetID)
  
  if(Add==T){
    points(x=gsdata[gsdata$position==1, ]$x, y=gsdata[gsdata$position==1, ]$runningScore,pch=16, lwd=1, cex=cex,bty="n", 
           ylim=c(-1,1), main=main, xlab="", ylab="Enrichment Score",
           xaxt='n', col=color)
    if(add_line==T){points(x=gsdata$x, y=gsdata$runningScore, type = "l",lty=1, col=color)}
    
  }
  if(Add==F){
    par(mar=c(5,5,5,5), las=2)
    plot(x=gsdata[gsdata$position==1, ]$x, y=gsdata[gsdata$position==1, ]$runningScore,pch=16, lwd=1, cex=cex,bty="n", 
         ylim=c(-1,1), main=main, xlab="", ylab="Enrichment Score",
         xaxt='n',col=color, xlim=c(head(gsdata$x,1),tail(gsdata$x,1)))
    if(add_line==T){points(x=gsdata$x, y=gsdata$runningScore, type = "l",lty=1,col=color)}
    abline(h=0, lty=1, lwd=1)
  }
}
GSEA_sample=function(x){
  #x=sample_data
  selection=input$Choice_GSEA
  condition=input$Condition_GSEA
  
  GSEA_input=x@GSEA[[selection]][[condition]]
}

##----------------------------------------------------------------------------##
##  Tab Outputs
##----------------------------------------------------------------------------##
output$GSEA_Plot <- renderPlot({

  gseaplot_DHH(GSEA_sample(sample_data()), input$PW)

 
})
output$Sort_GSEA_Plot <- renderPlot({
  
  col=brewer.pal(9,"Set1")
  print(input$PW_X1)
  for(i in 1:length(input$PW_X1)){
    if(i==1){ gseaplot_DHH(GSEA_sample(sample_data()), input$PW_X1[i], color = col[i]) }
    else{ gseaplot_DHH(GSEA_sample(sample_data()), input$PW_X1[i], color = col[i], Add = T) }
  }


  
  
})
output$GSEA_out_tab <- DT::renderDT({
  
  GSEA_sample(sample_data())@result[,c(4,5,6,7)]
  
  },
  options = list(
    autoWidth = T,scrollX=T,
    columnDefs = list(list(width = '10px', targets = "_all"))
    )
)




# Downloads of Plots  
output$downloadplot_GSEA_Plot <- downloadHandler(
  filename = function() { paste(input$PW, '.pdf', sep='') },
  content = function(file) {pdf(file, useDingbats = F)
    gseaplot_DHH(GSEA_sample(sample_data()), input$PW)
    dev.off()},
  contentType = "application/pdf"
)

output$downloadSort_GSEA_Plot <- downloadHandler(
  filename = function() { paste(input$PW, '.pdf', sep='') },
  content = function(file) {
    pdf(file, useDingbats = F)
    col=brewer.pal(9,"Set1")
    print(input$PW_X1)
    for(i in 1:length(input$PW_X1)){
      if(i==1){ gseaplot_DHH(GSEA_sample(sample_data()), input$PW_X1[i], color = col[i]) }
      else{ gseaplot_DHH(GSEA_sample(sample_data()), input$PW_X1[i], color = col[i], Add = T) }
    }
    dev.off()},
  contentType = "application/pdf"
)


