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

checkMissingInput = function(input){
  errors = c();
  if (is.null(input$clusterCols)){
    errors = c(errors,"No clustering parameters selected");
  }
  if (is.null(input$crossValidationFolds)){
    errors = c(errors,"Number of cross validation folds not specified");
  }
  return(errors);
}
