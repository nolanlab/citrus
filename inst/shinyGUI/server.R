shinyServer(function(input, output) {
  
  
  #currentGroupNames = reactive({})
  #currentClusterParameters = reactive({getParameterIntersections(input,fileList,fileCols)})
  #currentSelectedFiles = reactive({getSelectedFiles(input)})
  
  output$groupNameInput = renderUI({
    return(tagList(lapply(1:input$numberOfGroups,serialGroupNameInput)))
  })
  
  output$sampleGroupSelector = renderUI({
    return(tagList(lapply(getGroupNames(input),serialGroupSelectors,fileList=fileList)))
  })
  
  output$sampleGroupsTable = renderTable({
    if (preload){
      return(keyFile)
    } else {
      return(getAssignmentsTable(input,fileList))  
    }
  })
  
  # Estimates the number of events to be clustered. May be
  # less if files have fewer events than specified sample size
  output$estimatedClusteredEvents = renderUI({
    #selectedFiles = currentSelectedFiles()
    selectedFiles = getSelectedFiles(input)
    eventEstimate = format(input$fileSampleSize*length(unlist(selectedFiles)),big.mark=",",scientific=F)
    return(tags$div(paste0("Estimated maximum number of events to be clustered: ",eventEstimate)))
  })
  
  output$estimatedClusterSize = renderUI({
    numFiles = length(unlist(getSelectedFiles(input)))
    #numFiles = length(unlist(currentSelectedFiles()))
    eventEstimate = input$fileSampleSize*numFiles
    minClusterSize = format(floor(eventEstimate*(input$minimumClusterSizePercent/100)),big.mark=",",scientific=F)
    return(tags$div(paste0("Estimated minimum cluster size: ",minClusterSize," cells")))
  })
  
  # Transform cofactor 
  output$transformCofactor = renderUI({
    if (is.null(input$transformCols)||(length(input$transformCols)==0)){
      tags$input(type="text",id="transformCofactor",value="NULL",style="display:none")
    } else {
      numericInput(inputId="transformCofactor",label="Transform Cofactor",value=5,min=1)  
    }
  })
  
  output$conditionComparaMatrixInput = renderUI({
    if (preload){
      conditions = c("index",colnames(keyFile[,-labelCol]))
      rowList = list()
      for (row in conditions){
        colList = list()
        for (col in conditions){
          if ((col=="index")&&(row=="index")){
            colList[[col]] = tags$td(tags$b("Condition",class="blank"))  
          } else if (col=="index"){
            colList[[col]] = tags$th(row)  
          } else if (row=="index"){
            colList[[col]] = tags$th(col)  
          } else {
            if (col==row){
              colList[[col]] = tags$td(checkboxInput(inputId=paste(row,col,sep="_vs_"),label="",value=T))
            } else {
              colList[[col]] = tags$td(checkboxInput(inputId=paste(row,col,sep="_vs_"),label="",value=F))  
            }
            
          }
        }
        rowList[[row]] = tags$tr(tagList(colList))
      }
      tagList(rowList)
      return(tags$table(tagList(rowList),class="comparaTable"))  
    }
    return("SHOULD NOT BE USED");
  })
  
  output$clusterCols = renderUI({
    choices = isolate(getParameterIntersections(input,fileList,fileCols))
    #choices = isolate(currentClusterParameters())
    if (is.null(choices)){
      return(tagList(tags$b("Assign files to groups to enable selection of clustering parameters.")))
    } else {
      if ((!is.null(input$selectAllCluster))&&(input$selectAllCluster%%2 == 0)){
        checked = checkboxGroupInput("clusterCols",label="Clustering Parameters",choices=choices,selected=choices)
      } else {
        checked = checkboxGroupInput("clusterCols",label="Clustering Parameters",choices=choices)
      }
      return(checked)
    }
  })
  
  output$transformCols = renderUI({
    choices = isolate(getParameterIntersections(input,fileList,fileCols))
    #choices = isolate(currentClusterParameters())
    if (is.null(choices)){
      return(tagList(tags$b("Assign samples to groups to enable selection of transform parameters.")))
    } else {
      if ((!is.null(input$selectAllTransform))&&(input$selectAllTransform%%2 == 0)){
        checked = checkboxGroupInput("transformCols",label="Transform Parameters",choices=choices,selected=choices)
      } else {
        checked = checkboxGroupInput("transformCols",label="Transform Parameters",choices=choices)
      }
      return(checked)
      
    }
  })
  
  output$scaleCols = renderUI({
    choices = isolate(getParameterIntersections(input,fileList,fileCols))
    #choices = isolate(currentClusterParameters())
    if (is.null(choices)){
      return(tagList(tags$b("Assign samples to groups to enable selection of scaling parameters.")))
    } else {
      if ((!is.null(input$selectAllScale))&&(input$selectAllScale%%2 == 0)){
        checked = checkboxGroupInput("scaleCols",label="Scale Parameters",choices=choices,selected=choices)
      } else {
        checked = checkboxGroupInput("scaleCols",label="Scale Parameters",choices=choices)
      }
      return(checked)
    }
  })
  
  
  output$calculatedFeatures = renderUI({
    return(tagList(
                    tags$span("Calculated Cluster Features:",class="control-label"),
                    tags$br(),
                    tagList(lapply(citrus.featureTypes(),serialFeaturesInput))
    ))
  })
  
  output$medianCols = renderUI({
    if ((!is.null(input$featureType))&&(input$featureType=="medians")){
      choices = isolate(getParameterIntersections(input,fileList,fileCols))
      #choices = isolate(currentClusterParameters())
      if (is.null(choices)){
        return(tagList(tags$b("Assign samples to groups to enable selection of median parameters.")))
      } else {
        return(tagList(tags$hr(),checkboxGroupInput("medianColumns",label="Cluster Median Parameters:",choices=choices)))  
      }
    } else {
      return(tags$span(""))  
    }
  })
  
  
  output$crossValidationRange = renderUI({
    selectedFiles = getSelectedFiles(input)
    #selectedFiles = currentSelectedFiles()
    nFiles = length(unlist(selectedFiles))
    
    if ((length(names(selectedFiles))<2)||(!all(unlist(lapply(selectedFiles,length))>1))){
      return(tagList(tags$b("Assign samples to groups to enable specification of cross-validation folds")))
    } else {
      print(paste("FILES:",nFiles))
      return(tagList(numericInput(inputId="crossValidationFolds",label="Cross Validation Folds",value=1,min=1,max=nFiles)))
    }
  })
  
  output$groupSummary = renderUI({
    #selectedFiles=currentSelectedFiles()
    selectedFiles = getSelectedFiles(input)
    groupNames = getGroupNames(input)
    #groupNames = currentGroupNames()
    return(
      tagList(
        tags$ul(lapply(groupNames,serialGroupSummary,selectedFiles=selectedFiles))
      )
    )
  })
  
  output$conditionSummary = renderUI({
    if (preload){
      comparaConditions = getComparaConditions(input,conditions=colnames(keyFile[,-labelCol]))
      if (length(comparaConditions)==0){
        return(tags$ul(tags$li(tags$span("None Conditions Selected",class="red-error"))))
      }
      return(tags$ul(lapply(comparaConditions,tags$li)))
    } else {
      return(tags$ul(tags$li("Default Condition")))
    }
  })
  
  output$clusteringSummary = renderUI({
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
  
  output$workingDirectorySummary = renderUI({
    return(tagList(tags$span(dataDir),tags$br(),tags$br()))
  })
  
  output$featureSummary = renderUI({
    featureSetTags = tags$span("None",class="red-error")
    features=list();
    if (!is.null(input$featureType)&&(input$featureType=="abundances")){
      features[["Abundances"]] = tags$li("Cluster Abundances")
    } else if (!is.null(input$featureType)&&(input$featureType=="medians")){
      medianColumns = input$medianColumns
        medianVals = tags$span("No Median Parameters Selected",class="red-error");
        if (length(input$medianColumns)>0){
          medianVals = tags$span(paste(input$medianColumns,collapse=", "))
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
          tags$li(paste("Minimum Cluster Size: ",input$minimumClusterSizePercent,"%",sep="")),
          tags$li(featureSetTags)
          )
        )
    )
  })
  
  output$classificationSummary = renderUI({
    if (is.null(input$crossValidationFolds)){
      cvTag = tagList(tags$span("Cross Validation Folds:"),tags$span("None",class="red-error"))
    } else {
      cvTag = tagList(tags$span("Cross Validation Folds:"),tags$span(input$crossValidationFolds))
    }
    if (sum(getSelectedModels(input))==0){
      mTag = tagList(tags$span("Two-Class Models:"),tags$span("None",class="red-error"))
    } else {
      mTag = tags$span(paste("Two-Class Models:",paste(citrus.modelTypes()[getSelectedModels(input)],collapse=", ")))
    }
    return(tags$ul(tagList(tags$li(cvTag),tags$li(mTag))))
    
  })
  
  output$classificationModels = renderUI({
    # This is a hack to disable glmnet until lasso multinomial regression is supported. 
    modelTypes = citrus.modelTypes()
    if (input$numberOfGroups>2){
      modelTypes = modelTypes[modelTypes!="glmnet"]
    }
    tagList(tags$span("Two-Class Models:"),lapply(modelTypes,serialClassificationModel))  
  })
  
  output$run = renderUI({
    errors = errorCheckInput(input)
    if ((!is.null(input$runCitrus))&&(input$runCitrus)){
      if (length(errors)==0){
        writeRunCitrusFile(input)  
      } else {
        stop(simpleError("Can't write runCitrus.R with errors."))
      }
      if (input$citrusRunAction %in% c("wrc","qar")){
        runCitrus<<-FALSE
        if (input$citrusRunAction=="qar"){
          runCitrus<<-TRUE
        } 
        stop(simpleWarning("GUI Setup Complete"))    
      } else {
        # Do something fancy
      }
      
    } else {
      if (length(errors)==0){
        return(actionButton(inputId="runCitrus",label="Run Citrus"))  
      } else {
        return(tagList(tags$em("The following problems must be corrected before running citrus:"),tags$ul(tagList(lapply(errors,tags$li)))))
      }  
    } 
    
  })
  
})

##############################
# UI OUTPUT CONTROLS 
##############################
serialClassificationModel = function(modelName){
  checkboxInput(modelName,modelName,value=F)
}

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
  
  if (preload){
    inputTag = textInput(inputId=paste("Group",x,"name",sep="_"),value=unique(keyFile[,labelCol])[x],label=paste("Group",x,sep=" "))
    inputTag = disableInput(inputTag)
  } else {
    inputTag = textInput(inputId=paste("Group",x,"name",sep="_"),value=paste("Group",x,sep=" "),label=paste("Group",x,sep=" "))
  }
  tags$td(inputTag)
}


serialGroupSelectors = function(groupName,fileList){
  
  if (preload){
    inputTag = selectInput(paste(groupName,"files",sep=""),label=paste(groupName,"samples"),selected=fileList[fileGroupAssignments==groupName],choices=fileList,multiple=T,selectize=F)
    inputTag = disableInput(inputTag)
  } else {
    inputTag = selectInput(paste(groupName,"files",sep=""),label=paste(groupName,"samples"),selected=fileList[grep(groupName,fileList,ignore.case=T)],choices=fileList,multiple=T,selectize=F)
  }
  tags$td(inputTag)
}

#############################################
#
# TEMPLATE CONTROLS
#
#############################################
writeRunCitrusFile = function(input,templateFile=NULL){
  templateData = reactiveValuesToList(input)
  templateData[["minimumClusterSizePercent"]] = templateData[["minimumClusterSizePercent"]]/100;
  templateData[["citrusVersion"]] = citrus.version();
  templateData[["preload"]]=preload
  templateData[["dataDir"]]=dataDir
  templateData[["computedFeatures"]] = names(getComputedFeatures(input))[unlist(getComputedFeatures(input))]
  templateData[["classificationModels"]] = citrus.modelTypes()[getSelectedModels(input)]
  if (preload){
    templateData[["keyFile"]]=keyFile
    templateData[["conditionComparaMatrix"]]=getConditionComparaMatrix(input,conditions=colnames(keyFile[,-labelCol]))
    templateData[["conditions"]]=colnames(keyFile[,-labelCol])
  }
  outputDir = file.path(dataDir,"citrusOutput")
  if (!file.exists(outputDir)){
    dir.create(file.path(dataDir,"citrusOutput"),showWarnings=F)
  }
  
  runCitrusTemplateFilePath = file.path(system.file(package="citrus"),"shinyGUI","runCitrus.template")
  if (exists("debugTemplate")&&(debugTemplate)){
    runCitrusTemplateFilePath = "/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI/runCitrus.template" 
  }
  
  brew(file=runCitrusTemplateFilePath,output=file.path(outputDir,"runCitrus.R"))
  
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
    fileGroupAssignments[fileList %in% reactiveValuesToList(input)[[paste(groupName,"files",sep="")]]]=groupName
  }
  return(data.frame("File"=fileList,"Group"=fileGroupAssignments));
}

getGroupNames = function(input){
  
  if (preload){
    return(unique(keyFile[,labelCol]))
  }

  vals = c()
  for (i in 1:input$numberOfGroups){
    name = paste("Group",i,"name",sep="_");
    if (name %in% names(input)){
      vals[i]=input[[name]]  
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
    sf[[groupName]] = reactiveValuesToList(input)[[paste(groupName,"files",sep="")]]
  }
  return(sf)
}

getParameterIntersections = function(input,fileList,fileCols){
  selectedFiles = unlist(getSelectedFiles(input))
  params = Reduce(intersect,fileCols[which(fileList %in% selectedFiles)])
  names(params) = names(fileCols[[which(fileList %in% selectedFiles)[1]]])
  return(params)
}

stringQuote = function(x){
  return(paste("\"",x,"\"",sep=""));  
}



getComputedFeatures = function(input){
  features = list();
  featureSelections = reactiveValuesToList(input)
  for (type in citrus.featureTypes()){
    if ((paste("compute",type,sep="") %in% names(featureSelections))&&(featureSelections[[paste("compute",type,sep="")]])){
      features[[type]]=T
    } else {
      features[[type]]=F
    }
  }
  return(features)
}

getSelectedModels = function(input){
  selectedModels = rep(FALSE,length(citrus.modelTypes()))
  names(selectedModels) = citrus.modelTypes();
  input = reactiveValuesToList(input)
  for (modelType in citrus.modelTypes()){
    if ((modelType %in% names(input))&&(input[[modelType]])){
      selectedModels[[modelType]]=T
    }
  }
  return(selectedModels)
}

errorCheckInput = function(input){
  errors = c();
  if (preload){
    if (length(getComparaConditions(input,conditions=colnames(keyFile[,-labelCol])))==0){
      errors = c(errors,"No conditions selected for analysis")
    }
  }
  if (is.null(input$clusterCols)){
    errors = c(errors,"No clustering parameters selected");
  }
  if (is.null(input$crossValidationFolds)){
    errors = c(errors,"Number of cross validation folds not specified");
  }
  if (is.null(input$featureType)){
    errors = c(errors,"No computed cluster features selected")
  } else {
    if ((input$featureType=="medians")&&(length(input$medianColumns)==0)){
      errors = c(errors,"No cluster median parameters selected")
    }
  }
  selectedFiles = getSelectedFiles(input)
  counts = unlist(lapply(selectedFiles,length))
  if ((length(counts)<2)||any(counts<2)){
    errors = c(errors,"2 or more samples must be assigned to each group")
  }
  if (!any(getSelectedModels(input))){
    errors = c(errors,"At least one differential model must be selected")
  }
  
  
  
  return(errors);
}

getComparaConditions = function(input,conditions){
  input = reactiveValuesToList(input)
  comparaConditions = c()
  for (condition1 in conditions){
    for (condition2 in conditions){
      inputName = paste(condition1,condition2,sep="_vs_")
      if (inputName %in% names(input)&&(input[[inputName]])){
        if (condition1==condition2){
          comparaConditions = c(comparaConditions,condition1)
        } else {
          comparaConditions = c(comparaConditions,paste(condition2,condition1,sep=" vs. "))  
        }
        
      }
    }
  }
  return(comparaConditions)
}

getConditionComparaMatrix = function(input,conditions){
  input = reactiveValuesToList(input)
  comparaMatrix = matrix(F,nrow=length(conditions),ncol=length(conditions),dimnames=list(conditions,conditions))
  for (condition1 in conditions){
    for (condition2 in conditions){
      inputName = paste(condition1,condition2,sep="_vs_")
      if (inputName %in% names(input)&&(input[[inputName]])){
        comparaMatrix[condition1,condition2]=T 
      }
    }
  }
  return(comparaMatrix)
}