\name{citrus.buildFoldFeatureDifferences}
\alias{citrus.buildFoldFeatureDifferences}
\title{Calculate differences in feature values between two conditions}
\usage{
  citrus.buildFoldFeatureDifferences(index, features,
    conditions, citrus.dataArray, folds,
    calculateLeaveoutData = F)
}
\arguments{
  \item{index}{The fold index of files to be subtracted}

  \item{features}{An array of features from two conditions}

  \item{conditions}{A vector of 2 conditions to calculate
  feature differences between.}

  \item{citrus.daataArray}{A citrus.dataArray object
  containing the data used to calculate he features}

  \item{folds}{A list of fold indices}

  \item{calculateLeaveoutData}{A logical value indicating
  whether or not the difference should be calulated for
  files in the fold. If TRUE, then calculates differences
  for files in the fold. If FALSE, caluclates differences
  for files not in the fold.}
}
\description{
  \code{citrus.buildFoldFeatureDifferences} takes feature
  values from two conditions and returns the difference.
  Essentially it is simple array subtraction using names
  for indices.
}
\author{
  Robert Bruggner
}
\references{
  http://github.com/nolanlab/citrus/
}

