library(shiny)



serialGroupNameInput = function(x){
  textInput(inputId=paste("Group",x,"name",sep="_"),label=paste("Group",x,"name",sep=" "),value=paste("Group",x,sep=" "))  
}

getGroupNames = function(input){
  inputList = as.list(input)
  vals = c()
  for (i in 1:input$numberOfGroups){
    vals[i]=inputList[[paste("Group",i,"name",sep="_")]]
  }
  names(vals) = vals
  return(vals)
}

serialGroupSelectors = function(groupName,fileList){
  selectInput(paste(groupName,"files",sep=""),label=groupName,selected=fileList[grep(groupName,fileList,ignore.case=T)],choices=fileList,multiple=T)
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
      return(
        tagList(
          tags$b("No files assigned to groups.")
        ) 
      )
    } else {
      return(selectInput("clusterCols",label="Clustering Parameters",choices=choices,multiple=T))  
    }
  })
  
  output$transformCols = reactiveUI(function(){
    choices = getParameterIntersections(input,fileList,fileCols);
    if (is.null(choices)){
      return(
        tagList(tags$b("No files assigned to groups.")) 
      )
    } else {
      return(selectInput("transformCols",label="Transform Parameters",choices=choices,multiple=T))  
    }
  })
  
  output$medianCols = reactiveUI(function(){
    if (input$calcMedianFeatures){
      choices = getParameterIntersections(input,fileList,fileCols);
      if (is.null(choices)){
        return(
          tagList(tags$b("No files assigned to groups.")) 
        )
      } else {
        return(selectInput("medianCols",label="Cluster Median Parameters",choices=choices,multiple=T))  
      }
    } else {
      return(tags$b(""))  
    }
    cat(paste("M:",input$calcMedianFeatures,"\n"))
    
  })
  
})

