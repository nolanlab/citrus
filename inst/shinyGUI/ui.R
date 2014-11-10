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
               numericInput(inputId="fileSampleSize","Events Sampled Per File:",min=1,value=1000),
               uiOutput("estimatedClusteredEvents"),
               tags$hr(),
               tags$div(numericInput("minimumClusterSizePercent",label="Minimum Cluster Size: A percentage of aggregate data size",min=0.1,max=100,step=0.1,value=5),style="width:300px;"),
               uiOutput("estimatedClusterSize"),
               tags$hr(),
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
                    ),
                    tags$td(
                     tagList(
                       uiOutput("scaleCols"),
                       actionButton("selectAllScale",label="Select All/None")
                     )
                    ),
                    tags$td(
                      uiOutput("transformCofactor")
                    )
                  )
               )
            ),
      
      tabPanel("Cluster Characterization",
                radioButtons(inputId="featureType",label="Feature Types:",choices=citrus.featureTypes(),selected="abundances"),
                uiOutput("medianCols")
               ),
      
      tabPanel("Association Model Configuration",
               uiOutput("crossValidationRange"),
               tags$hr(),
               uiOutput("classificationModels")
               ),
      
      tabPanel("Run!",
               
               radioButtons("citrusRunAction", "Citrus Execution Options:",
                            list("Quit GUI and run Citrus in R" = "qar","Write runCitrus.R file to data directory only" = "wrc")
                            ),
               tags$hr(),
               uiOutput("run"),
               tags$hr(),
               tags$em("Multithreading Options:",class="control-label"),
               checkboxInput(inputId="coreLimit",label="Limit Clustering Multicore Usage"),
               numericInput(inputId="analysisCores",label="Number of Cores:",value=1,min=1,max=16)
               
               )
      )
      
    # End Main Panel
    )
  
))


