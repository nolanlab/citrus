rm(list=ls())

library("citrus")
library("shiny")
library("brew")

dataDir="~/Desktop/work/citrus/inst/extdata/example1/"

debugTemplate=F

runApp("inst/shinyGUI")