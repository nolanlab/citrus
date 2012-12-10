library(shiny)


disabledCheckbox = function(inputId,label){
  tagList(
    tags$input(type="checkbox",id=inputId,disabled="disabled"),
    tags$span(label,class="disabled"),
    tags$br()
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
    tags$em("Classification Summary:"),
    uiOutput("classificationSummary"),
    tags$hr(),
    tags$em("Sample Summary:"),
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
                tags$div(sliderInput("minClusterSize",label="Minimum Cluster Size: A percentage of aggregate data size",min=1,max=100,step=1,value=5),style="width:300px;"),
                tags$hr(),
                checkboxGroupInput(inputId="computedFeatures",label="Computed Cluster Features:",choices=c("Cluster Densities","Cluster Medians"),selected="Cluster Densities"),
                uiOutput("medianCols")
               ),
      
      tabPanel("Classification Setup",
               uiOutput("crossValidationRange"),
               checkboxGroupInput(inputId="classificationModelTypes",label="Classification Models:",choices=c("PAMR","GLMNET"),selected=c("PAMR","GLMNET"))
               ),
      
      tabPanel("Run!",
               disabledCheckbox(inputId="multithread",label="Run Multithreaded"),
               disabledCheckbox(inputId="exportClusters",label="Export Identified Clusters"),
               disabledCheckbox(inputId="optimisticMode",label="Run in Naive Mode"),
               uiOutput("quitAndRun")
               )
      )
      
    # End Main Panel
    )
  
))


