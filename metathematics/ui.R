library(shiny)
library(plotly)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  titlePanel('Metagenomic Thematic Structure: Gevers Data'),
  
  plotOutput('est',height='100px'),
  plotlyOutput('tax'),
  plotlyOutput('fxn'),
  tableOutput('tbl')
  
))
