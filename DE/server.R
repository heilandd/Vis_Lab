##----------------------------------------------------------------------------##
## Server-DE
##----------------------------------------------------------------------------##


##----------------------------------------------------------------------------##
## UiOutputs
##----------------------------------------------------------------------------##
output$Choice <- renderUI({
selectInput(inputId="Choice", 
            label="Select Condition Factor:",
            choices = c(sample_data()@condition_Analysis),
            selected= c(sample_data()@condition_Analysis)[1]
)
})

output$C1 <- renderUI({
  selectInput(inputId="C1", 
              label="Condition Factor 1:",
              choices = c(levels(sample_data()@fdata$group),levels(sample_data()@fdata$treatment)),
              selected= c(levels(sample_data()@fdata$group),levels(sample_data()@fdata$treatment))[1]
  )
})


output$C2 <- renderUI({
  selectInput(inputId="C2", 
              label="Condition Factor 2:",
              choices = c(levels(sample_data()@fdata$group),levels(sample_data()@fdata$treatment)),
              selected= c(levels(sample_data()@fdata$group),levels(sample_data()@fdata$treatment))[2]
  )
})

output$Gene <- renderUI({
  selectInput(inputId="Gene", 
              label="Select Gene:",
              choices = c(rownames(sample_data()@data@counts)),
              selected= NULL
  )
})

output$Number_of_genes <- renderUI({
  selectInput(inputId="Number_of_genes", 
              label="Number of genes for Heatmap:",
              choices = c(10,20,30,40,50,60,70,80,90,100),
              selected= 20
  )
})




##----------------------------------------------------------------------------##
## Functions for Tab
##----------------------------------------------------------------------------##
require(dplyr)
require(DESeq2)

diff_gene <- reactive({
  
  #call DE from input
  select=input$Choice
  print(select)
  diff_gene=results(sample_data()@DE[[select]], contrast = c(select, input$C1, input$C2))
  print(diff_gene)
  return(diff_gene)
  
})
plot_DE=function(diff_gene){
  plot(x=diff_gene$log2FoldChange, y=c(-log10(diff_gene$padj)), pch=19, bty="o", xlab="Differential Expression", ylab="-log10 FDR")
  abline(v=1, lty=2);abline(v=-1, lty=2);abline(h=-log10(0.05), lty=2)
  diff_gene[is.na(diff_gene$padj), ]$padj=1
  points(x=diff_gene[diff_gene$log2FoldChange>1 & diff_gene$padj<0.05 | diff_gene$log2FoldChange<c(-1) & diff_gene$padj<0.05 , ]$log2FoldChange, 
         y=c(-log10(diff_gene[diff_gene$log2FoldChange>1 & diff_gene$padj<0.05 | diff_gene$log2FoldChange<c(-1) & diff_gene$padj<0.05 , ]$padj)), 
         pch=19, col="red")
  diff_gene=diff_gene[order(diff_gene$padj, decreasing = F), ]
  text(head(diff_gene, 10)$log2FoldChange, -log10(head(diff_gene, 10)$padj)+0.4, labels = rownames(head(diff_gene, 10)), cex=0.7) 
}
get_matrix=function(x){
  counts=sample_data()@data@norm_Exp
  feat=sample_data()@fdata
  
  #counts=sample_data@data@counts
  #feat=sample_data@fdata
  
  print(counts[1:4,1:4])
  print(dim(counts))
  
  
  
  #Selected condtions from the Check
  G1=levels(feat$group)
  
  print(c(G1))
  
  ##Select Samples to test
  feat=feat[feat$group %in% c(G1), ]
  counts=counts[,rownames(feat)]
  
  print(feat)
  print(dim(counts))
  
  gene=input$Gene
  print(gene)
  
  
  mat=matrix(NA, nrow(feat),  length(unique(feat$group)))
  print(mat)
  
  lv=unique(feat$group)
  
  for(i in 1:ncol(mat)){
    mat[,i]=c(as.numeric(counts[gene, as.character(feat[feat$group==lv[i], ]$Sample)]), rep(NA,nrow(mat)-length(as.numeric(counts[gene, as.character(feat[feat$group==lv[i], ]$Sample)]))))  
  }
  
  colnames(mat)=paste0("Group_", lv)
  
  print(mat)
  return(mat)
}
vioplotDHH=function (x, range = 1.5,mar=c(6,4,4,4), h = NULL, ylim = NULL, names = colnames(x), 
                     horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                     lwd = 1, rectCol = "white", colMed = "white", pchMed = 19, 
                     at, add = FALSE, wex = 1, drawRect = TRUE, main="main", ylab="") 
{
  par(mar=mar,las=1)
  datas = lapply(1:ncol(x), function(i){as.numeric(na.omit(x[,i]))})
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    require(sm)
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    
    if (!horizontal) {
      if (!add) {
        plot(NA,xlim = xlim, ylim = ylim, main=main, xaxt="n", bty="n", xlab="", ylab=ylab)
        axis(2)
        
        text(1:length(names),par("usr")[3]-0.3, srt = 60, adj= 1, xpd = TRUE,labels = names, cex=0.8)
      }
      
      for (i in 1:n) {
        polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
                lty = lty, lwd = lwd)
        if (drawRect) {
          lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                lty = lty)
          rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
               q3[i], col = rectCol)
          points(at[i], med[i], pch = pchMed, col = colMed)
        }
      }
    }
  else {
    if (!add) {
      plot(NA,xlim = ylim, ylim = xlim,main=main, bty="n",xlab="", ylab="",xaxt="n")
      axis(1)
      
      text(1:length(names),par("usr")[3]-0.3, srt = 60, adj= 1, xpd = TRUE,labels = names, cex=0.8)
    }
    
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}
Barplot_SC=function(input,min_x=0,Text_cex=1, points=T, 
                    main="main",
                    ylab="",
                    ylim=c(0,c(max(input,na.rm = TRUE)+3)),
                    mar=c(6,4,4,4), 
                    col_Schema=F,
                    col_Self=c("red", "grey"),
                    p_values=T,
                    AbstandxAxe=0.2,
                    AbstandBalken=0.1,
                    AbstandPWert=0.1){
  
  
  
  
  if(class(input)=="matrix"){len=ncol(input)}else{print("Input needs to be matrix")}
  par(mar=mar,xaxs = "i",yaxs = "i")
  maxx=0
  for(i in 1:ncol(input)){maxx=c(maxx,max(na.omit(input[,i])))}
  
  plot(x=c(1,len+1), y=c(c(min_x-c(max(maxx)/29)),max(maxx)), bty="n", type="n",ylab=ylab, xlab="",xaxt="n", main=main,ylim=ylim)
  
  if (col_Schema==T){
    require(RColorBrewer); rf <- colorRampPalette(rev(brewer.pal(9,'Set1')))
    r <- rf(len)
    col=sample(r) 
  } else {col=col_Self}
  
  #plot
  for(i in 1: len){
    val=na.omit(input[,i])
    xx=c(i,i,i+0.8,i+0.8 )
    yy=c(min_x,mean(val), mean(val), min_x)
    polygon(x=xx,y=yy, col=col[i], border="black", cex=1)
  }
  for(i in 1: len){
    val=na.omit(input[,i])
    xx=c(i+0.4,i+0.4 )
    yy=c(mean(val), mean(val)+sd(val))
    polygon(x=xx,y=yy, col="black",cex=1)
    
    xx=c(i+0.3,i+0.5 )
    yy=c(mean(val)+sd(val), mean(val)+sd(val))
    polygon(x=xx,y=yy, col="black")
    
    
  }
  if(points==T){
    for(i in 1: len){
      val=na.omit(input[,i])
      points(x=rnorm(mean=i+0.4, sd=0.05, n=length(val)), y=as.numeric(val))
    }
  }
  
  if (p_values==T){
    if(ncol(input)>2){
      for(i in 1:c(ncol(input)-1)){
        p=t.test(as.numeric(na.omit(input[,i])),as.numeric(na.omit(input[,i+1])))$p.value
        text(x=c(i+0.8) ,y=c(max(input,na.rm = TRUE)+i/AbstandPWert), labels = paste("p=",round(p, digits = 3),cex=Text_cex, sep=""))
        polygon(x=c(i+0.4, i+1.4),y=c(max(input,na.rm = TRUE)+(i/AbstandBalken),max(input,na.rm = TRUE)+(i/AbstandBalken) ), col="black")
      }
      
    }
    else{
      p=t.test(as.numeric(na.omit(input[,1])),as.numeric(na.omit(input[,2])))$p.value
      text(x=1.8 ,y=c(max(input,na.rm = TRUE)+AbstandPWert), labels = paste("p=",round(p, digits = 3),cex=Text_cex, sep=""))
      polygon(x=c(1.4, 2.4),y=c(max(input,na.rm = TRUE)+AbstandBalken,max(input,na.rm = TRUE)+AbstandBalken ), col="black")
    }
  }
  
  
  #put Axisi
  polygon(x=c(0,len),y=c(0,0), col="black",cex=1)
  for(i in 1:len){
    polygon(y=c(c(-max(maxx)/35),max(maxx)/35),x=c(i+0.4,i+0.4), col="black",cex=1)
  }
  
  #input names
  
  text(seq(1,len,by=1)+0.4, par("usr")[3]-AbstandxAxe, srt = 60, adj= 1, xpd = TRUE,labels = colnames(input), cex=Text_cex)
}

heatmapDHH_DE_tab=function(x,DE){
  
  #Define input
  select=input$Choice
  #Select samples that were used in Analnalysis

  groups=levels(x@fdata[,select])[levels(x@fdata[,select]) %in%   c(input$C1, input$C2) ]
  samples_in_use=levels(x@fdata$Sample)[x@fdata[,select] %in% groups]
  
  print(samples_in_use)
  
  #Expressionfile
  exp=x@data@norm_Exp[,samples_in_use]
  
  #DE_Genes
  DE=data.frame(DE)
  DE$Genes=rownames(DE)
  
  #Number of genes that were printed in heatmap
  nr=as.numeric(input$Number_of_genes)
  
  
  
  #select Genes by mosly diff expressed genes
  DE=DE[DE$pvalue<0.05, ]
  print(head(DE))
  selected_Genes=c(head(DE[order(DE$log2FoldChange, decreasing = T), ]$Genes, nr/2),
                   head(DE[order(DE$log2FoldChange, decreasing = F), ]$Genes, nr/2) )
  print(selected_Genes)
  require(pheatmap);require(viridis)
  print(nr)
  
  pheatmap(exp[rev(selected_Genes),samples_in_use ],
           cluster_row = T,
           cluster_cols = F,
           scale="row",
           color = viridis(50),
           cutree_rows = 2,
           cutree_cols = 2)
}




##----------------------------------------------------------------------------##
##  Tab Outputs
##----------------------------------------------------------------------------##
output$Volcano_plot <- renderPlot({

  #diff_gene=as.data.frame(diff_gene())
  #print(head(diff_gene))
  #Print volcano plot
  plot_DE(diff_gene())

 
})
output$plot_violin <- renderPlot({
  ### Get Groups that need to be tested
  par(mar=c(5,5,5,5), las=1)
  vioplotDHH(get_matrix(sample_data()), 
             col = RColorBrewer::brewer.pal(9,"Set1"), 
             main=input$Gene,
             ylab="Normalized Geneexpression")
  })
output$Barplot <- renderPlot({
  ### Get Groups that need to be tested
  mat=get_matrix(sample_data())
  mat=mat/mat[1,1]
  par(mar=c(5,5,5,5), las=1)
  Barplot_SC(mat, col_Schema = T, main=input$Gene, ylab="Normalized Geneexpression")
  
  
})
output$DE_Tab <- DT::renderDT({
  out=data.frame(diff_gene())
  out=data.frame(Genes=rownames(out),FC= out$log2FoldChange, Adj_p_value=out$padj)
  
  },
  options = list(
    autoWidth = T,scrollX=T,
    columnDefs = list(list(width = '30px', targets = "_all"))
  ),rownames=F
  )


output$Heatmap <- renderPlot({
  heatmapDHH_DE_tab(x=sample_data(), DE=diff_gene())
})



#Render Outputs on Top
output$DE <- renderValueBox({
  valueBox(
    value = nrow(data.frame(diff_gene()) %>% filter(pvalue<0.05)),
    subtitle = "Number of Significant DE",
    color = "light-blue"
  )})
output$FDR <- renderValueBox({
  valueBox(
    value = nrow(data.frame(diff_gene()) %>% filter(padj<0.05)),
    subtitle = "FDR corrected",
    color = "light-blue"
  )})
output$TotalGenes <- renderValueBox({
  valueBox(
    value = nrow(data.frame(diff_gene())),
    subtitle = "Number of Genes",
    color = "light-blue"
  )})


# Downloads of Plots  
output$downloadplot_violin <- downloadHandler(
  filename = function() { paste(gene, '.pdf', sep='') },
  content = function(file) {pdf(file, useDingbats = F)
    ar(mar=c(5,5,5,5), las=1)
    vioplotDHH(get_matrix(sample_data()), 
               col = RColorBrewer::brewer.pal(9,"Set1"), 
               main=input$Gene,
               ylab="Normalized Geneexpression")
    dev.off()},
  contentType = "application/pdf"
)

output$downloadVolcano_plot <- downloadHandler(
  filename = function() { paste(c(select, input$C1, input$C2), '.pdf', sep='') },
  content = function(file) {
    pdf(file, useDingbats = F)
    plot_DE(diff_gene())
    dev.off()},
  contentType = "application/pdf"
)

output$downloadHeatmap <- downloadHandler(
  filename = function() { paste(c(select, input$C1, input$C2), '.pdf', sep='') },
  content = function(file) {
    pdf(file, useDingbats = F)
    heatmapDHH_DE_tab(x=sample_data(), DE=diff_gene())
    dev.off()},
  contentType = "application/pdf"
)


output$downloadplot_bar <- downloadHandler(
  filename = function() { paste(gene, '.pdf', sep='') },
  content = function(file) {
    pdf(file, useDingbats = F)
    mat=get_matrix(sample_data())
    mat=mat/mat[1,1]
    par(mar=c(5,5,5,5), las=1)
    Barplot_SC(mat, col_Schema = T, main=input$Gene, ylab="Normalized Geneexpression")
    dev.off()},
  contentType = "application/pdf"
)


