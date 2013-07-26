rm(list=ls(all=T))
library('flowCore')

inputDir = "~/Desktop/work/citrus/inst/extdata/example2/"
outputDir = "~/Desktop/work/citrus/inst/extdata/example2.1/"
for (fcsFileName in list.files(inputDir,pattern=".fcs")){
  fcsFile = read.FCS(file.path(system.file(package="citrus"),"extdata","example2",fcsFileName))
  fileData = exprs(fcsFile)
  for (popId in unique(fileData[,"trueLabel"])){
    popData = fileData[fileData[,"trueLabel"]==popId,]
    newFCS = flowFrame(exprs=popData)
    write.FCS(newFCS,filename=file.path(outputDir,sub(".fcs",paste("_pop",popId,".fcs",sep=""),fcsFileName)))
  }
}