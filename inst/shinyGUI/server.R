library(shiny)

serialGroupSummary = function(groupName,selectedFiles){
  if (is.null(groupName)){
    return(tags$b(""))
  }
  countTag=tags$span("0",class="red")
  if ((length(selectedFiles)>0) && (groupName %in% names(selectedFiles))){
    if (length(selectedFiles[[groupName]])>1){
      countTag = tags$span(length(selectedFiles[[groupName]]))  
    } else {
      countTag = tags$span(length(selectedFiles[[groupName]]),class="red")
    }
    
  }
  return(tags$li(tagList(tags$span(paste(groupName,"Samples: ")),countTag)))
}

serialGroupNameInput = function(x){
  textInput(inputId=paste("Group",x,"name",sep="_"),label=paste("Group",x,"name",sep=" "),value=paste("Group",x,sep=" "))  
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

serialGroupSelectors = function(groupName,fileList){
  selectInput(paste(groupName,"files",sep=""),label=paste(groupName,"samples"),selected=fileList[grep(groupName,fileList,ignore.case=T)],choices=fileList,multiple=T)
}

getAssignmentsTable = function(input,fileList){
  fileGroupAssignments = rep("",length(fileList))
  for (groupName in getGroupNames(input)){
    fileGroupAssignments[fileList %in% as.list(input)[[paste(groupName,"files",sep="")]]]=groupName
  }
  return(data.frame("File"=fileList,"Group"=fileGroupAssignments));
}

getClusterCols = function(fileName,dataDir){
  colnames(read.FCS(paste(dataDir,fileName,sep=""),which.lines=1))
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

dataDir = "/Users/rbruggner/Desktop/work/citrus/data/syntheticData/train/unstim/"
fileList = list.files(dataDir)
fileGroupAssignments = rep("",length(fileList))  
fileCols = lapply(fileList,getClusterCols,dataDir=dataDir)

shinyServer(function(input, output) {
  
  output$groupNameInput = reactiveUI(function() {
    return( 
        tagList(
          sapply(1:input$numberOfGroups,serialGroupNameInput)
        )
    )
  })
  
  output$sampleGroupSelector = reactiveUI(function(){
    return(
      tagList(sapply(getGroupNames(input),serialGroupSelectors,fileList=fileList))
    )
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
      return(tagList(tags$b("Assign files to groups to enable selection of transform parameters.")))
    } else {
      return(selectInput("transformCols",label="Transform Parameters",choices=choices,multiple=T))  
    }
  })
  
  output$medianCols = reactiveUI(function(){
    cat(paste(input$computedFeatures,"\n"));
    if ("Cluster Medians" %in% input$computedFeatures){
      choices = getParameterIntersections(input,fileList,fileCols);
      if (is.null(choices)){
        return(tagList(tags$b("Assign files to groups to enable selection of median parameters.")))
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
      return(tagList(tags$b("Need at least 2 samples assigned to each group")))
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
      ccTag = tags$li(tagList(tags$span("Clustering Parameters:"),tags$span("None Selected",class="red")))
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
    return(tags$ul(tagList(tags$li(paste("Minimum Cluster Size: ",as.numeric(input$minClusterSize)*100,"%",sep="")))))
  })
  
})

