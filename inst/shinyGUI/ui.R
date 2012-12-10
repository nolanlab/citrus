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
    tags$head(tags$link(rel="stylesheet",type="text/css",href="citrus.css")),
    h4("Citrus Runtime Summary:"),
    tags$hr(),
    tags$em("Working Directory:"),
    uiOutput("workingDirectorySummary"),
    tags$em("Groups Summary:"),
    uiOutput("groupSummary"),
    tags$em("Clustering Summary:"),
    uiOutput("clusteringSummary"),
    tags$em("Cluster Characterization Summary:"),
    uiOutput("featureSummary"),
    tags$em("Classification Summary"),
    tags$hr(),
    tableOutput("sampleGroupsTable")
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Sample Group Setup",
               tags$table(class="sampleGroupTable",
                 tags$tr(tagList(
                   tags$td(tagList(
                     numericInput(inputId="numberOfGroups",label="Number Of Sample Groups",value=2,min=2),
                     uiOutput("groupNameInput"))
                   ),
                   tags$td(
                     uiOutput("sampleGroupSelector")
                   )
                 ))
                )
      ),
               
      tabPanel("Clustering Setup",
               numericInput(inputId="fileSampleSize","Events Sampled Per File",min=1,value=1000),
               uiOutput("clusterCols"),
               uiOutput("transformCols")
               ),
      
      tabPanel("Cluster Characterization",
                sliderInput("minClusterSize",label="Minimum Cluster Size as a percentage of aggregate data",min=0.01,max=1,step=0.01,value=0.05),
                tags$hr(),
                checkboxGroupInput(inputId="computedFeatures",label="Computed Cluster Features:",choices=c("Cluster Densities","Cluster Medians"),selected="Cluster Densities"),
                uiOutput("medianCols")
               ),
      
      tabPanel("Classification Setup",
               uiOutput("crossValidationRange"),
               checkboxInput(inputId="buildPAMRModel",label="Construct PAMR Model",value=T),
               checkboxInput(inputId="buildGLMModel",label="Construct GLMNet Model",value=T)
               )
      )
      
    # End Main Panel
    )
  
))


