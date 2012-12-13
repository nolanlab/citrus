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

##############################
# UI OUTPUT CONTROLS 
##############################
serialGroupSummary = function(groupName,selectedFiles){
  if (is.null(groupName)){
    return(tags$b(""))
  }
  countTag=tags$span("0",class="red-error")
  if ((length(selectedFiles)>0) && (groupName %in% names(selectedFiles))){
    if (length(selectedFiles[[groupName]])>1){
      countTag = tags$span(length(selectedFiles[[groupName]]))  
    } else {
      countTag = tags$span(length(selectedFiles[[groupName]]),class="red-error")
    }
    
  }
  return(tags$li(tagList(tags$span(paste(groupName,"Samples: ")),countTag)))
}

serialGroupNameInput = function(x){
  tags$td(textInput(inputId=paste("Group",x,"name",sep="_"),label=paste("Group",x,"name",sep=" "),value=paste("Group",x,sep=" ")))
}

serialGroupSelectors = function(groupName,fileList){
  tags$td(selectInput(paste(groupName,"files",sep=""),label=paste(groupName,"samples"),selected=fileList[grep(groupName,fileList,ignore.case=T)],choices=fileList,multiple=T))
}

writeRunCitrusFile = function(input,templateFile=NULL){
  templateData = as.list(input)
  templateData[["dataDir"]]=dataDir
  outputDir = file.path(dataDir,"citrusOutput")
  if (!file.exists(outputDir)){
    dir.create(file.path(dataDir,"citrusOutput"),showWarnings=F)
  }
  brew(
    file=file.path(system.file(package="citrus"),"shinyGUI","runCitrus.template"),
    output=file.path(outputDir,"runCitrus.R")
  )
}

convertToLabeledFileList = function(x){
  df = data.frame(labels=c(),defaultCondition=c())
  for (group in names(x)){
    df = rbind(df,data.frame(labels=rep(group,length(x[[group]])),defaultCondition=x[[group]]))
  }
  df[,"labels"]=as.factor(df[,"labels"])
  return(df)
}

convertColToDefinition = function(colname,df){
  paste(colname,"=c(",paste(sapply(df[,colname],stringQuote),collapse=","),")",sep="")
}

getAssignmentsTable = function(input,fileList){
  fileGroupAssignments = rep("",length(fileList))
  for (groupName in getGroupNames(input)){
    fileGroupAssignments[fileList %in% as.list(input)[[paste(groupName,"files",sep="")]]]=groupName
  }
  return(data.frame("File"=fileList,"Group"=fileGroupAssignments));
}

getGroupNames = function(input){
  inputList = as.list(input)
  vals = c()
  for (i in 1:input$numberOfGroups){
    name = paste("Group",i,"name",sep="_");
    if (name %in% names(inputList)){
      vals[i]=inputList[[name]]  
    } else {
      vals[i] = "EmptyGroup";
    }
    
  }
  names(vals) = vals
  return(vals)
}

getSelectedFiles = function(input){
  sf = list();
  for (groupName in getGroupNames(input)){
    sf[[groupName]] = as.list(input)[[paste(groupName,"files",sep="")]]
  }
  return(sf)
}

getParameterIntersections = function(input,fileList,fileCols){
  selectedFiles = unlist(getSelectedFiles(input))
  return(Reduce(intersect,fileCols[which(fileList %in% selectedFiles)]))
}

stringQuote = function(x){
  return(paste("\"",x,"\"",sep=""));  
}

getComputedFeatures = function(input){
  features = list();
  if (!is.null(input$computeDensities)&&input$computeDensities){
    features[["densities"]] = TRUE
  } else {
    features[["densities"]] = FALSE
  }
  if (!is.null(input$computeMedians)&&input$computeMedians){
    features[["medians"]] = TRUE
  } else {
    features[["medians"]] = FALSE
  }
  return(features)
}

errorCheckInput = function(input){
  errors = c();
  if (is.null(input$clusterCols)){
    errors = c(errors,"No clustering parameters selected");
  }
  if (is.null(input$crossValidationFolds)){
    errors = c(errors,"Number of cross validation folds not specified");
  }
  computedFeatures = getComputedFeatures(input)
  if (all(!unlist(computedFeatures))){
    errors = c(errors,"No computed cluster features selected")
  } else {
    if (computedFeatures[["medians"]]&&(length(input$medianCols)==0)){
      errors = c(errors,"No cluster median parameters selected")
    }
  }
  selectedFiles = getSelectedFiles(input)
  counts = unlist(lapply(selectedFiles,length))
  if ((length(counts)<2)||any(counts<2)){
    errors = c(errors,"2 or more samples must be assigned to each group")
  }
  
  return(errors);
}
