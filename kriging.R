rm(list = ls())
check = "check"
Library = '/home-nfs/ollie/cluettig/R/x86_64-redhat-linux-gnu-library/3.4’
library(doParallel, lib.loc = Library)
library(gstat, lib.loc = Library)
library(sp, lib.loc = Library)
library(rgdal, lib.loc = Library)
library(raster, lib.loc = Library)
setwd("/work/ollie/cluettig/github/Interpolation/Python")
source("fitmodel.R")


######### INPUT ###################################################################################################

args <- commandArgs(TRUE)
InputFileName <- args[1]
RadiusFileName <- args[2]
OutputFileName <- args[3]
Directory <- args[4]
#InputFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/vx_input.tif"
#RadiusFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/kriging_radius.tif"
#Directory <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test"
#OutputFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/vx_int.tif"

######### FLAGS ###################################################################################################

RandomKrigingSubset = FALSE
BoundaryKrigingSubset = TRUE
VariogramSubset = TRUE

######### READ ###################################################################################################

InputRaster <- raster(InputFileName)
RadiusRaster <- raster(RadiusFileName)

######### FORMAT DATA ############################################################################################
GivenPoints <- as.data.frame(InputRaster, xy = TRUE, na.rm = TRUE) # values and coordinates of given positions
Radius <- as.data.frame(RadiusRaster)
AllPoints <- as.data.frame(InputRaster, xy = TRUE)
InputLayerName <- strsplit(tail(strsplit(InputFileName,"/")[[1]], n = 1),"\\.")[[1]][1]
SearchedPoints <- AllPoints[(is.na(AllPoints[InputLayerName])) & Radius$kriging_radius != 0,]
SearchedPoints[InputLayerName] <- NULL 
SearchedPoints <- cbind(SearchedPoints, "z" = Radius[(is.na(AllPoints[InputLayerName])) & Radius$kriging_radius != 0,])  #Coordinates and kriging radius of searched positions

if (BoundaryKrigingSubset || RandomKrigingSubset) {
    # Possibility to take only boundary points for kriging:
  if (BoundaryKrigingSubset) {
    BoundaryRaster <- boundaries(clump(InputRaster))
    HelpKrigingSample <- as.data.frame(BoundaryRaster)
    HelpKrigingSample[is.na(HelpKrigingSample)] <- FALSE
    KrigingSample <- AllPoints[HelpKrigingSample == 1,]
  } else {
    KrigingSample <- GivenPoints
  }
    # Possibility to take only a subsample for kriging:
  if (RandomKrigingSubset) {
    SampleSize <- 1000
    if (SampleSize < nrow(KrigingSample)) {
      KrigingSample <- KrigingSample[sample(nrow(KrigingSample), SampleSize),]
    }
  }
} else {
    KrigingSample <- GivenPoints
}

# Possibility to take only a subsample for variogram estimation:
if (VariogramSubset) {
  SampleSize <- 1000
  if (SampleSize < nrow(KrigingSample)) {
    VariogramSample <- KrigingSample[sample(nrow(KrigingSample), SampleSize),]
  } else {
    VariogramSample <- KrigingSample
  }
} else {
  VariogramSample <- KrigingSample
}
names(VariogramSample)[names(VariogramSample) == InputLayerName] <- "z"
names(GivenPoints)[names(GivenPoints) == InputLayerName] <- "z"
names(KrigingSample)[names(KrigingSample) == InputLayerName] <- "z"

rm(Radius, AllPoints, RadiusRaster, RadiusFileName, InputFileName, HelpKrigingSample, BoundaryRaster, InputRaster)


################ TREND REMOVAL ###########################################################################################
# print("trend removal")
png(paste(Directory, "/correlation.png", sep = ""))
plot(VariogramSample)
dev.off()

Model <- fitmodel(unlist(VariogramSample["x"]),unlist(VariogramSample["y"]), unlist(VariogramSample["z"]), Directory)
Model$coefficients[is.na(Model$coefficients)] <- 0
if (length(Model$coefficients) == 16) {
  print("cubic trend model")
  VariogramSample["z"] <- VariogramSample["z"] - (Model$coefficients[1] 
    + Model$coefficients[2] * VariogramSample["x"]
    + Model$coefficients[3] * VariogramSample["x"]^2 
    + Model$coefficients[4] * VariogramSample["x"]^3 
    + Model$coefficients[5] * VariogramSample["y"] 
    + Model$coefficients[6] * VariogramSample["y"]^2
    + Model$coefficients[7] * VariogramSample["y"]^3 
    + Model$coefficients[8] * VariogramSample["x"] * VariogramSample["y"] 
    + Model$coefficients[9] * VariogramSample["x"]^2 * VariogramSample["y"]
    + Model$coefficients[10] * VariogramSample["x"]^3 * VariogramSample["y"] 
    + Model$coefficients[11] * VariogramSample["x"] * VariogramSample["y"]^2 
    + Model$coefficients[12] * VariogramSample["x"]^2 * VariogramSample["y"]^2
    + Model$coefficients[13] * VariogramSample["x"]^3 * VariogramSample["y"]^2 
    + Model$coefficients[14] * VariogramSample["x"] * VariogramSample["y"]^3 
    + Model$coefficients[15] * VariogramSample["x"]^2 * VariogramSample["y"]^3
    + Model$coefficients[16] * VariogramSample["x"]^3 * VariogramSample["y"]^3)
  KrigingSample["z"] <- KrigingSample["z"] - (Model$coefficients[1]
    + Model$coefficients[2] * KrigingSample["x"] 
    + Model$coefficients[3] * KrigingSample["x"]^2
    + Model$coefficients[4] * KrigingSample["x"]^3 
    + Model$coefficients[5] * KrigingSample["y"]
    + Model$coefficients[6] * KrigingSample["y"]^2
    + Model$coefficients[7] * KrigingSample["y"]^3 
    + Model$coefficients[8] * KrigingSample["x"] * KrigingSample["y"] 
    + Model$coefficients[9] * KrigingSample["x"]^2 * KrigingSample["y"]
    + Model$coefficients[10] * KrigingSample["x"]^3 * KrigingSample["y"]
    + Model$coefficients[11] * KrigingSample["x"] * KrigingSample["y"]^2 
    + Model$coefficients[12] * KrigingSample["x"]^2 * KrigingSample["y"]^2
    + Model$coefficients[13] * KrigingSample["x"]^3 * KrigingSample["y"]^2 
    + Model$coefficients[14] * KrigingSample["x"] * KrigingSample["y"]^3 
    + Model$coefficients[15] * KrigingSample["x"]^2 * KrigingSample["y"]^3
    + Model$coefficients[16] * KrigingSample["x"]^3 * KrigingSample["y"]^3)
} else if (length(Model$coefficients) == 9) {
  print("quadratic trend model")
  VariogramSample["z"] <- VariogramSample["z"] - (Model$coefficients[1] 
    + Model$coefficients[2] * VariogramSample["x"] 
    + Model$coefficients[3] * VariogramSample["x"]^2
    + Model$coefficients[4] * VariogramSample["y"] 
    + Model$coefficients[5] * VariogramSample["y"]^2 
    + Model$coefficients[6] * VariogramSample["x"]*VariogramSample["y"]
    + Model$coefficients[7] * VariogramSample["x"]^2*VariogramSample["y"] 
    + Model$coefficients[8] * VariogramSample["x"]*VariogramSample["y"]^2
    + Model$coefficients[9] * VariogramSample["x"]^2*VariogramSample["y"]^2)
  KrigingSample["z"] <- KrigingSample["z"] - (Model$coefficients[1] 
    + Model$coefficients[2] * KrigingSample["x"] 
    + Model$coefficients[3] * KrigingSample["x"]^2
    + Model$coefficients[4] * KrigingSample["y"] 
    + Model$coefficients[5] * KrigingSample["y"]^2 
    + Model$coefficients[6] * KrigingSample["x"] * KrigingSample["y"]
    + Model$coefficients[7] * KrigingSample["x"]^2 * KrigingSample["y"]
    + Model$coefficients[8] * KrigingSample["x"] * KrigingSample["y"]^2 
    + Model$coefficients[9] * KrigingSample["x"]^2 * KrigingSample["y"]^2)
} else if (length(Model$coefficients) == 4) {
  print("linear trend model")
  VariogramSample["z"] <- VariogramSample["z"] - (Model$coefficients[1] 
    + Model$coefficients[2] * VariogramSample["x"]
    + Model$coefficients[3] * VariogramSample["y"] 
    + Model$coefficients[4] * VariogramSample["x"] * VariogramSample["y"])
  KrigingSample["z"] <- KrigingSample["z"] - (model$coefficients[1] 
    + model$coefficients[2] * KrigingSample["x"]
    + model$coefficients[3] * KrigingSample["y"]
    + model$coefficients[4] * KrigingSample["x"] * KrigingSample["y"])
} else {
  print("Error while substracting the trend model")
}
# Plot of variables after removing the trend:
png(paste(Directory, "/correlation_trend_removal.png", sep = ""))
plot(VariogramSample)
dev.off()

########## VARIOGRAM ESTIMATION ###################################################################################################
coordinates(VariogramSample) <- ~x + y
# d <- data.frame(x,y,z)
# coordinates(d) = ~x+y
G <- gstat(formula = z~x+y,data = VariogramSample)
MeasuredVariogram <- variogram(G, width = 10, cutoff = 10000000)
png(paste(Directory, "/variogram.png", sep = ""))
plot(MeasuredVariogram)
dev.off()
EstimatedVariogram <- vgm(psill = max(MeasuredVariogram$gamma)/2, "Sph", range = max(MeasuredVariogram$dist)/3) # ,anis = c(45,0.3))
FittedVariogram <- fit.variogram(MeasuredVariogram, EstimatedVariogram)
#Plot of the fitted variogram function:
png(paste(Directory, "/variogram_fit.png", sep = ""))
plot(MeasuredVariogram, FittedVariogram)
dev.off()
rm(MeasuredVariogram, EstimatedVariogram, G, VariogramSample)

########## KRIGING #############################################################################################################
coordinates(KrigingSample) <- ~x+y
kriging <- function(i){
    if (unique(SearchedPoints$z)[i] < 25000) {max_kriging_points = 3000}else{max_kriging_points = 500}
    ActualRadius <- unique(SearchedPoints$z)[i]
    ActualPoints <- SearchedPoints[SearchedPoints$z == ActualRadius,]
    coordinates(ActualPoints) <- ~x + y
    KrigingValues <- krige(formula = z~1, locations = KrigingSample, newdata = ActualPoints, model = FittedVariogram, maxdist = ActualRadius, nmax = 1000)
    ActualPoints@data$z <- KrigingValues@data$var1.pred
    ActualPoints <- data.frame(ActualPoints)[,1:3]
}
print("number of cores:")
NumberOfCores <- detectCores()
print(NumberOfCores)
Cluster <- makeCluster(NumberOfCores)
registerDoParallel(Cluster)
MaxGapLength <- length(SearchedPoints$z)/NumberOfCores
NumberOfIterations <- length(unique(SearchedPoints$z))
if (MaxGapLength > 1) {
  # divide gaps with to many points in parts with a maximal length of all searched points divided by number of cores
  for (j in unique(SearchedPoints$z)) {
    if (length(SearchedPoints[SearchedPoints$z == j,]) > MaxGapLength) {
      print("maxgaplength überschritten")
      print(length(SearchedPoints[SearchedPoints$z == j,]))
      print(MaxGapLength)
      NumberOfParts <- MaxGapLength/length(SearchedPoints[SearchedPoints$z == j,])
      NumberOfIterations <- NumberOfIterations + (NumberOfParts - 1)
      print(NumberOfParts)
      for (k in 1:floor(NumberOfParts)) {  
        SearchedPoints[SearchedPoints$z == j,]$z[(k - 1) * MaxGapLength + 1:MaxGapLength] <- SearchedPoints[SearchedPoints$z == j,]$z[(k - 1) * MaxGapLength + 1:MaxGapLength] + k * 0.01
      }
    }
  }
}
print("number of iterarions:")
print(NumberOfIterations)
KrigingResults <- foreach(Iteration = 1:NumberOfIterations,
                           .combine = list,
                           .multicombine = TRUE,
                           .packages = c("gstat", "sp"))  %dopar%
  kriging(Iteration)

stopImplicitCluster()
rm(SearchedPoints)
Result <- KrigingResults[[1]]
for (i in 2:length(KrigingResults)) {
  Result <- rbind(Result,KrigingResults[[i]])
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

############# WRITE ##################################################################################################
# d_z_n <- data.frame(x_n,y_n,z_n)
# d_z_n <- rename(d_z_n, c("x_n" = "x", "y_n" = "y", "z_n" = "z"))
Result <- rbind(Result, GivenPoints)
coordinates(Result) <- ~ x + y
gridded(Result) <- TRUE
ResultRaster <- raster(Result)
writeRaster(ResultRaster, filename = OutputFileName, format = "GTiff", overwrite = TRUE)
