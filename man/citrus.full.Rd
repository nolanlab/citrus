\name{citrus.full}
\alias{citrus.full}
\title{Driver to run a full citrus analysis}
\usage{
  citrus.full(dataDir, outputDir, clusterCols,
    fileSampleSize, fileList, nFolds,
    modelTypes = c("pamr", "glmnet"),
    featureTypes = c("densities"),
    minimumClusterSizePercent = 0.05, transformCols = NULL,
    ...)
}
\arguments{
  \item{dataDir}{The path to the data directory containing
  the FCS files to be analyzed.}

  \item{outputDir}{Path to a directory where the citrus
  output will be placed.}

  \item{clusterCols}{A vector of integers or parameter
  names that should be used for clustering of the data}

  \item{fileSampleSize}{Number of events to be sampled from
  each analyzed file. Files with fewer events contribute
  all of their events.}

  \item{fileList}{A matrix containing file names,
  conditions, and group assignments. See details.}

  \item{nFolds}{Number of cross validation folds used to
  assess model accuracy}

  \item{modelTypes}{A vector of classification models to
  construct. Valid arguments are \code{pamr} and
  \code{glmnet}.}

  \item{featureTypes}{A vector of descriptive feature types
  to be calculated for each cluster. Valid arguments are
  \code{densities} and \code{medians}. See details.}

  \item{minimumClusterSizePercent}{Specifies the minimum
  cluster size to be analyzed as a percentage of the total
  aggregate datasize. A value etween \code{0} and
  \code{1}.}

  \item{transformCols}{A vector of integer or parameter
  names to be transformed before analysis.}

  \item{...}{Further arguments to be passed to Citrus
  subcomponents.}
}
\description{
  \code{citrus.full} runs a full Citrus analysis, including
  clustering, feature extraction and regularized
  classification. Plots are generated that show
  classification cross validation error rates, informative
  classifier features, and phenotype plots for informative
  clusters.
}
\details{
  Details about the cluster conditions matrix, fold
  features, etc.
}
\author{
  Robert Bruggner
}
\references{
  http://github.com/nolanlab/citrus/
}

