#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
source('smth.R')

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Smoothing Span Selection"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("band",
                        label = "choose a band",
                        choices = list("Sentinel-1 VH",
                                       "Sentinel-1 VV",
                                       "Sentinel-2 Blue",
                                       "Sentinel-2 Green",
                                       "Sentinel-2 Red",
                                       "Sentinel-2 NIR",
                                       "Sentinel-2 Red Edege 1",
                                       "Sentinel-2 Red Edege 2",
                                       "Sentinel-2 Red Edege 3",
                                       "Sentinel-2 Red Edege 4",
                                       "Sentinel-2 SWIR1",
                                       "Sentinel-2 SWIR2"),
                        selected = "Sentinel-2 Red"),
            sliderInput("span1",
                        "smoothing span 1:",
                        min = 0.1,
                        max = 0.9,
                        value = 0.25),
            sliderInput("span2",
                        "smoothing span 2:",
                        min = 0.1,
                        max = 0.9,
                        value = 0.5),
            sliderInput("span3",
                        "smoothing span 3:",
                        min = 0.1,
                        max = 0.9,
                        value = 0.75)
        ),
            
        

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("smthPlot"),
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  dat <- reactive({
    if (input$band == "Sentinel-1 VH"){
      file_name <- 'data/nmds_VH.csv'
      band_name <- 'VH'
    }else if(input$band == "Sentinel-1 VV"){
      file_name <- 'data/nmds_VV.csv'
      band_name <- 'VV'
    }else if(input$band == "Sentinel-2 Blue"){
      file_name <- 'data/nmds_BLUE.csv'
      band_name <- 'Blue'
    }else if(input$band == "Sentinel-2 Green"){
      file_name <- 'data/nmds_GREEN.csv'
      band_name <- 'Green'
    }else if(input$band == "Sentinel-2 Red"){
      file_name <- 'data/nmds_RED.csv'
      band_name <- 'Red'
    }else if(input$band == "Sentinel-2 NIR"){
      file_name <- 'data/nmds_NIR.csv'
      band_name <- 'Nir'
    }else if(input$band == "Sentinel-2 Red Edege 1"){
      file_name <- 'data/nmds_RE1.csv'
      band_name <- 'RE1'
    }else if(input$band == "Sentinel-2 Red Edege 2"){
      file_name <- 'data/nmds_RE2.csv'
      band_name <- 'RE2'
    }
    else if(input$band == "Sentinel-2 Red Edege 3"){
      file_name <- 'data/nmds_RE3.csv'
      band_name <- 'RE3'
    }else if(input$band == "Sentiel-2 Red Edege 4"){
      file_name <- 'data/nmds_RE4.csv'
      band_name <- 'RE4'
    }
    else if(input$band == "Sentiel-2 SWIR1"){
      file_name <- 'data/nmds_SWIR1.csv'
      band_name <- 'SWIR1'
    }else if(input$band == "Sentinel-2 SWIR2"){
      file_name <- 'data/nmds_SWIR2.csv'
      band_name <- 'SWIR2'
    }
    
    span1 <- input$span1
    span2 <- input$span2
    span3 <- input$span3
    
    df <- smth_df(file_name, band_name, span1, span2, span3)
    df
 
  })
  


    output$smthPlot <- renderPlot({
      if (input$band == "Sentinel-1 VH"){
        band_name <- 'VH'
      }else if(input$band == "Sentinel-1 VV"){
        band_name <- 'VV'
      }else if(input$band == "Sentinel-2 Blue"){
        band_name <- 'Blue'
      }else if(input$band == "Sentinel-2 Green"){
        band_name <- 'Green'
      }else if(input$band == "Sentinel-2 Red"){
        band_name <- 'Red'
      }else if(input$band == "Sentinel-2 NIR"){
        band_name <- 'Nir'
      }else if(input$band == "Sentinel-2 Red Edege 1"){
        band_name <- 'RE1'
      }else if(input$band == "Sentinel-2 Red Edege 2"){
        band_name <- 'RE2'
      }
      else if(input$band == "Sentinel-2 Red Edege 3"){
        band_name <- 'RE3'
      }else if(input$band == "Sentiel-2 Red Edege 4"){
        band_name <- 'RE4'
      }
      else if(input$band == "Sentiel-2 SWIR1"){
        band_name <- 'SWIR1'
      }else if(input$band == "Sentinel-2 SWIR2"){
        band_name <- 'SWIR2'
      }
      eval(parse(text=paste('ggplot(data = dat(), aes(x = DOY, y = ', band_name, ', group = span, color = span)) + 
                            geom_line()', sep = '')))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
