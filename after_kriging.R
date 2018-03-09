rm(list = ls())
library(rgdal)
library(raster)
setwd("~/Documents/Daten/Antarctica/Recovery/vel/Interpolation")
load("image.RData")

KrigingResults <- readRDS("_rslurm_kriging_2/results_0.RDS")
Result <- do.call("rbind", KrigingResults)
rows <-0
for (i in 1:168){
    KrigingResults <- readRDS(paste("_rslurm_kriging_2/results_",i,".RDS",sep=""))
    Result_i <- do.call("rbind", KrigingResults)
    Result <- rbind(Result, Result_i)
}

# Addition of substracted trend model:
if (length(Model$coefficients) == 16){
 Result$z <- Result$z + (Model$coefficients[1]
   + Model$coefficients[2] * Result$x
   + Model$coefficients[3] * Result$x^2
   + Model$coefficients[4] * Result$x^3
   + Model$coefficients[5] * Result$y
   + Model$coefficients[6] * Result$y^2
   + Model$coefficients[7] * Result$y^3
   + Model$coefficients[8] * Result$x * Result$y
   + Model$coefficients[9] * Result$x^2 * Result$y
   + Model$coefficients[10] * Result$x^3 * Result$y
   + Model$coefficients[11] * Result$x * Result$y^2
   + Model$coefficients[12] * Result$x^2 * Result$y^2
   + Model$coefficients[13] * Result$x^3 * Result$y^2
   + Model$coefficients[14] * Result$x * Result$y^3
   + Model$coefficients[15] * Result$x^2 * Result$y^3
   + Model$coefficients[16] * Result$x^3 * Result$y^3)
} else if (length(Model$coefficients) == 9) {
 Result$z <- Result$z + (Model$coefficients[1]
   + Model$coefficients[2] * Result$x
   + Model$coefficients[3] * Result$x^2
   + Model$coefficients[4] * Result$y
   + Model$coefficients[5] * Result$y^2
   + Model$coefficients[6] * Result$x * Result$y
   + Model$coefficients[7] * Result$x^2 * Result$y
   + Model$coefficients[8] * Result$x * Result$y^2
   + Model$coefficients[9] * Result$x^2 * Result$y^2)
} else if (length(Model$coefficients) == 4) {
 Result$z <- Result$z + (Model$coefficients[1]
   + Model$coefficients[2] * Result$x
   + Model$coefficients[3] * Result$y
   + Model$coefficients[4] * Result$x * Result$y)
} else {
print("Error while addition of the model")
}

Result <- rbind(Result, GivenPoints)
coordinates(Result) <- ~ x + y
gridded(Result) <- TRUE
ResultRaster <- raster(Result)
writeRaster(ResultRaster, filename = "result.tif", format = "GTiff", overwrite = TRUE)
