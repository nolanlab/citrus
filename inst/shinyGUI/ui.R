library(shiny)

actionButton <- function(inputId, label) {
  tagList(
    singleton(tags$head(tags$script(src = "js/actionbutton.js"))),
    tags$button(id=inputId, type="button", class="btn action-button", label)
  )
}

eventButton <- function(inputId, value) {
  tagList(
    singleton(tags$head(tags$script(src = "js/eventbutton.js"))),
    tags$button(id = inputId,
                class = "eventbutton btn",
                type = "button",
                as.character(value))
  )
}


shinyUI(pageWithSidebar(
  
  headerPanel("Citrus GUI v0.01"),
  
  sidebarPanel(
    h4("Runtime Summary:"),
    tableOutput("sampleGroupsTable")
  ),
  
  mainPanel(
    tabsetPanel(
      
      tabPanel("Sample Group Setup",
               numericInput(inputId="numberOfGroups",label="Number Of Sample Groups",value=2,min=2),
               uiOutput("groupNameInput")
               ), 
      
      tabPanel("Sample Group Assignment", 
               uiOutput("sampleGroupSelector")
               ),
      
      tabPanel("Clustering Configuration",
               numericInput(inputId="fileSampleSize","Events Sampled Per File",min=1,value=1000),
               uiOutput("clusterCols"),
               uiOutput("transformCols")
               ),
      
      tabPanel("Cluster Characterization",
                sliderInput("minClusterSize",label="Minimum Cluster Size as a percentage of aggregate data",min=0.01,max=1,step=0.01,value=0.05),
                wellPanel(
                tags$b("Computed Cluster Features:"),
                tags$hr(),
                checkboxInput(inputId="calcDensityFeatures",label="Calculate Cluster Densities"),
                checkboxInput(inputId="calcMedianFeatures",label="Calculate Cluster Medians"),
                uiOutput("medianCols")
                )
               )
      
        
      
    )
  )
))


