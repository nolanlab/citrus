rm(list=ls())

library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("ggplot2")
#library("spade",lib.loc="~/Desktop/work/R/spade_with_nn/lib")
library("citrus")

dataDir="~/Desktop/work/citrus/inst/extdata/example1/"
runApp("inst/shinyGUI")