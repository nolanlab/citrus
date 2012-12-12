
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
