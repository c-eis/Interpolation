rm(list = ls())
check = 'check'
Library = '/home-nfs/ollie/cluettig/R/x86_64-redhat-linux-gnu-library/3.4’
library(doParallel, lib.loc = Library)
library(gstat, lib.loc = Library)
library(sp, lib.loc = Library)
library(rgdal, lib.loc = Library)
library(raster, lib.loc = Library)
setwd("/work/ollie/cluettig/github/Interpolation/Python")
source("fitmodel.R")


######### INPUT #############################################################################\
######################

args <- commandArgs(TRUE)
InputFileName <- args[1]
RadiusFileName <- args[2]
OutputFileName <- args[3]
Directory <- args[4]
#InputFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/vx_input.tif"
#RadiusFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/kriging_radiu\
s.tif"
#Directory <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test"
#OutputFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/vx_int.tif"

######### FLAGS #############################################################################\
######################

RandomKrigingSubset = FALSE
BoundaryKrigingSubset = TRUE
VariogramSubset = TRUE

######### READ ##############################################################################\
#####################

InputRaster <- raster(InputFileName)
RadiusRaster <- raster(RadiusFileName)

######### FORMAT DATA #######################################################################\
#####################
GivenPoints <- as.data.frame(InputRaster, xy = TRUE, na.rm = TRUE) # values and coordinates o\
f given positions
Radius <- as.data.frame(RadiusRaster)
AllPoints <- as.data.frame(InputRaster, xy = TRUE)
InputLayerName <- strsplit(tail(strsplit(InputFileName,"/")[[1]], n = 1),"\\.")[[1]][1]
SearchedPoints <- AllPoints[(is.na(AllPoints[InputLayerName])) & Radius$kriging_radius != 0,]
SearchedPoints[InputLayerName] <- NULL
SearchedPoints <- cbind(SearchedPoints, "z" = Radius