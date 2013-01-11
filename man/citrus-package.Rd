\name{citrus-package}
\alias{citrus-package}
\alias{citrus}
\docType{package}
\title{
Citrus: Identification of Stratifying Subpopulations in Flow Cytometry Data
}
\description{
Citrus atuomatically identifies populations of cells in multi-dimensional flow cytometry data whose behavior differs between sample groups. 
}
\details{
\tabular{ll}{
Package: \tab citrus\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2012-12-04\cr
License: \tab Use by author's permission only.\cr
}
Citrus identifies subpopulations of cells in multidimensional flow cytometry data by clustering cells using heiarchcial clustering, automaitcally characterizing them using a number of features, and using regularized supervised learning algorithms to identify the features that are best predictors of user-specified sample groups. A citrus analysis may be configured and launched via the \code{citrus.launchUI()} command, or by using the \code{citrus.full()} driver function.
}
\author{
Robert Bruggner

Maintainer: Robert Bruggner <rbruggner@gmail.com>
}
\references{
https://github.com/nolanlab/citrus
}
\keyword{ package }
\examples{
# Analyze package example data
dataDir=system.file("extdata","example1",package="citrus")

# Create output directory
outputDir = file.path(dataDir,"output")
dir.create(outputDir)

# Read in files to be analyzed and assign them to groups
fileList = cbind(labels=c(rep("healthy",10),rep("diseased",10)),
                 condition1=list.files(dataDir,pattern=".fcs"))

# Run a ctrus analysis
citrus.full(dataDir,
            outputDir=outputDir,
            clusterCols=c("Red","Blue"),
            fileSampleSize=500,
            fileList=fileList,
            nFolds=5)

# Inform user of output location
paste(paste("Output is in",outputDir))
}
