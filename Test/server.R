##----------------------------------------------------------------------------##
## Server-Load Data
##----------------------------------------------------------------------------##
# number of samples
output$load_data_number_of_cells <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) formatC(nrow(sample_data()@fdata), format = "f", big.mark = ",", digits = 0) else 0,
    subtitle = "Number of samples",
    color = "light-blue"
  )
})

output$load_data_number_of_Groups <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) formatC(length(unique(sample_data()@fdata$group)), format = "f", big.mark = ",", digits = 0) else 0,
    subtitle = "Number of groups",
    color = "light-blue"
  )
})

output$load_data_number_of_Batches <- renderValueBox({
  valueBox(
    value = if ( !is.null(sample_data()) ) formatC(length(unique(sample_data()@fdata$Batch)), format = "f", big.mark = ",", digits = 0) else 0,
    subtitle = "Number of batches",
    color = "light-blue"
  )
})




