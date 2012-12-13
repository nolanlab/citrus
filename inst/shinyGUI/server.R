shinyServer(function(input, output) {
  
  output$groupNameInput = reactiveUI(function() {
    return(tagList(lapply(1:input$numberOfGroups,serialGroupNameInput)))
  })
  
  output$sampleGroupSelector = reactiveUI(function(){
    return(tagList(lapply(getGroupNames(input),serialGroupSelectors,fileList=fileList)))
  })
  
  
  output$sampleGroupsTable = reactiveTable(function(){
    return(getAssignmentsTable(input,fileList))
  })
  
  output$clusterCols = reactiveUI(function(){
    choices = getParameterIntersections(input,fileList,fileCols);
    if (is.null(choices)){
      return(tagList(tags$b("Assign files to groups to enable selection of clustering parameters.")))
    } else {
      return(selectInput("clusterCols",label="Clustering Parameters",choices=choices,multiple=T))  
    }
  })
  
  output$transformCols = reactiveUI(function(){
    choices = getParameterIntersections(input,fileList,fileCols);
    if (is.null(choices)){
      return(tagList(tags$b("Assign samples to groups to enable selection of transform parameters.")))
    } else {
      return(selectInput("transformCols",label="Transform Parameters",choices=choices,multiple=T))  
    }
  })
  
  output$calculatedFeatures = reactiveUI(function(){
    return(tagList(
                    tags$span("Calculated Cluster Features:"),
                    tags$br(),
                    checkboxInput(inputId="computeDensities",label="Cluster Densities",value=T),
                    checkboxInput(inputId="computeMedians",label="Cluster Medians",value=F)
    ))
  })
  
  output$medianCols = reactiveUI(function(){
    if ((!is.null(input$computeMedians))&&(input$computeMedians)){
      choices = getParameterIntersections(input,fileList,fileCols);
      if (is.null(choices)){
        return(tagList(tags$b("Assign samples to groups to enable selection of median parameters.")))
      } else {
        return(selectInput("medianCols",label="Cluster Median Parameters",choices=choices,multiple=T))  
      }
    } else {
      return(tags$span(""))  
    }
  })
  
  output$crossValidationRange = reactiveUI(function(){
    selectedFiles = getSelectedFiles(input)
    nFiles = length(unlist(selectedFiles))
    
    if ((length(names(selectedFiles))<2)||(!all(unlist(lapply(selectedFiles,length))>1))){
      return(tagList(tags$b("Assign samples to groups to enable specification of cross-validation folds")))
    } else {
      return(tagList(numericInput(inputId="crossValidationFolds",label="Cross Validation Folds",value=min(5,nFiles),min=2,max=nFiles)))
    }
  })
  
  output$groupSummary = reactiveUI(function(){
    selectedFiles=getSelectedFiles(input)
    groupNames = getGroupNames(input)
    return(
      tagList(
        tags$ul(lapply(groupNames,serialGroupSummary,selectedFiles=selectedFiles))
      )
    )
  })
  
  output$clusteringSummary = reactiveUI(function(){
    if (is.null(input$clusterCols)){
      ccTag = tags$li(tagList(tags$span("Clustering Parameters:"),tags$span("None Selected",class="red-error")))
    } else {
      ccTag = tags$li(tagList(tags$span("Clustering Parameters:"),tags$span(paste(input$clusterCols,collapse=", "))))
    }
    if (is.null(input$transformCols)){
      tcTag = tags$li(tagList(tags$span("Transform Parameters:"),tags$span("None Selected")))
    } else {
      tcTag = tags$li(tagList(tags$span("Transform Parameters:"),tags$span(paste(input$transformCols,collapse=", "))))
    }
    return(tags$ul(tagList(tags$li(paste("Events sampled per file:",input$fileSampleSize)),ccTag,tcTag)));  
  })
  
  output$workingDirectorySummary = reactiveUI(function(){
    return(tagList(tags$span(dataDir),tags$br(),tags$br()))
  })
  
  output$featureSummary = reactiveUI(function(){
    featureSetTags = tags$span("None",class="red-error")
    features=list();
    if (!is.null(input$computeDensities)&&input$computeDensities){
      features[["Densities"]] = tags$li("Cluster Densities")
    }
    if (!is.null(input$computeMedians)&&input$computeMedians){
        medianCols = input$medianCols
        medianVals = tags$span("No Median Parameters Selected",class="red-error");
        if (length(input$medianCols)>0){
          medianVals = tags$span(paste(input$medianCols,collapse=", "))
        }
        features[["Medians"]] = tags$li(tagList(tags$span("Cluster Medians:"),medianVals))
    }
    if (length(features)>0){
      featureSetTags = tags$ul(tagList(features))
    } 
    featureSetTags = tagList(tags$span("Computed Cluster Features:"),featureSetTags)
    
    return(
      tags$ul(
        tagList(
          tags$li(paste("Minimum Cluster Size: ",input$minClusterSize,"%",sep="")),
          tags$li(featureSetTags)
          )
        )
    )
  })
  
  output$classificationSummary = reactiveUI(function(){
    if (is.null(input$crossValidationFolds)){
      cvTag = tagList(tags$span("Cross Validation Folds:"),tags$span("None",class="red-error"))
    } else {
      cvTag = tagList(tags$span("Cross Validation Folds:"),tags$span(input$crossValidationFolds))
    }
    if (is.null(input$classificationModelTypes)){
      mTag = tagList(tags$span("Classification Models:"),tags$span("None",class="red-error"))
    } else {
      mTag = tags$span(paste("Classification Models:",paste(input$classificationModelTypes,collapse=", ")))
    }
    return(tags$ul(tagList(tags$li(cvTag),tags$li(mTag))))
    
  })
  
  output$run = reactiveUI(function(){
    if ((!is.null(input$runCitrus))&&(input$runCitrus)){
      writeRunCitrusFile(input)
      stop(simpleWarning("GUI Setup Complete"))
    }
    errors = errorCheckInput(input)
    if (length(errors)==0){
      return(tagList(checkboxInput(inputId="runCitrus",label="<- Check to run Citrus"),tags$em("Checkbox to be replaced by button when input-specific reactivity becomes avaialble.")))  
    } else {
      return(tagList(tags$em("The following problems must be corrected before running citrus:"),tags$ul(tagList(lapply(errors,tags$li)))))
    }
  })
  
})

