library(shiny)

shinyUI(pageWithSidebar(
  
  headerPanel(paste("Citrus UI v",citrus.version(),sep=""),windowTitle="Citrus UI"),
  
  sidebarPanel(
    tags$head(tags$link(rel="stylesheet",type="text/css",href="citrus.css")),
    h4("Citrus Runtime Summary:"),
    tags$hr(),
    tags$em("Working Directory:"),
    uiOutput("workingDirectorySummary"),
    tags$em("Groups Summary:"),
    uiOutput("groupSummary"),
    tags$em("Condition Summary:"),
    uiOutput("conditionSummary"),
    tags$em("Clustering Summary:"),
    uiOutput("clusteringSummary"),
    tags$em("Cluster Characterization Summary:"),
    uiOutput("featureSummary"),
    tags$em("Class Summary:"),
    uiOutput("classificationSummary"),
    tags$hr(),
    tags$em("Sample Summary:"),
    tableOutput("sampleGroupsTable")
  ),
  
  mainPanel(
    tabsetPanel(
      
        tabPanel("Sample Group Setup",
                 if (preload){
                    disableInput(numericInput(inputId="numberOfGroups",label="Number Of Sample Groups",value=length(unique(keyFile[,labelCol]))))
                 } else {
                   numericInput(inputId="numberOfGroups",label="Number Of Sample Groups",value=2,min=2)
                   #disableInput(numericInput(inputId="numberOfGroups",label="Number Of Sample Groups",value=2))
                 },
                 tags$table(class="sampleGroupTable",
                            tagList(
                              tags$tr(uiOutput("groupNameInput")),
                              tags$tr(uiOutput("sampleGroupSelector"))
                            )),
                 if (preload){
                   tagList(tags$hr(),tags$label("Condition Comparisons"),uiOutput("conditionComparaMatrixInput"))
                 } else {
                   tags$br()
                 }
                 
        ),
               
      tabPanel("Clustering Setup",
               numericInput(inputId="fileSampleSize","Events Sampled Per File",min=1,value=1000),
               # Scale cols not attached to package yet. Uncomment when fixed.
               #tags$table(tags$tr(tags$td(uiOutput("clusterCols")),tags$td(uiOutput("transformCols")),tags$td(uiOutput("scaleCols"))))
               tags$table(
                 tags$tr(
                   tags$td(
                     tagList(
                       uiOutput("clusterCols"),
                       actionButton("selectAllCluster",label="Select All/None")
                      )
                    ),
                    tags$td(
                      tagList(
                        uiOutput("transformCols"),
                        actionButton("selectAllTransform",label="Select All/None")
                      )
                    )
                  )
               )
            ),
      
      tabPanel("Cluster Characterization",
                tags$div(numericInput("minimumClusterSizePercent",label="Minimum Cluster Size: A percentage of aggregate data size",min=0.1,max=100,step=0.1,value=5),style="width:300px;"),
                tags$hr(),
                uiOutput("calculatedFeatures"),
                uiOutput("medianCols"),
                uiOutput("emdCols")
               ),
      
      tabPanel("Two-Class Setup",
               uiOutput("crossValidationRange"),
               tags$hr(),
               uiOutput("classificationModels")
               ),
      
      tabPanel("Run!",
               
               radioButtons("citrusRunAction", "Citrus Execution Options:",
                            list("Write runCitrus.R file to data directory only" = "wrc",
                                 "Quit GUI and run Citrus in R" = "qar")),
               tags$hr(),
               uiOutput("run"),
               tags$hr(),
               tags$em("TBD Runtime Options:",class="control-label"),
               tags$br(),
               disableInput(checkboxInput(inputId="exportClusters",label="Export Identified Clusters"))
               )
      )
      
    # End Main Panel
    )
  
))


