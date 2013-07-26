rm(list=ls(all=T))
library('flowCore')

inputDir = "~/Desktop/work/citrus/inst/extdata/example1/"
outputDir = "~/Desktop/work/citrus/inst/extdata/example4.1/"
for (fcsFileName in list.files(inputDir,pattern=".fcs")){
  fcsFile = read.FCS(file.path(system.file(package="citrus"),"extdata","example1",fcsFileName))
  fileData = exprs(fcsFile)
  for (popId in 1:3){
    popData = fileData[fileData[,"TrueLabel"]==popId,]
    newFCS = flowFrame(exprs=popData)
    write.FCS(newFCS,filename=file.path(outputDir,sub(".fcs",paste("_pop",popId,".fcs",sep=""),fcsFileName)))
  }
}