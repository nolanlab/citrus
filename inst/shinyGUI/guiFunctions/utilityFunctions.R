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
