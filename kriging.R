rm(list = ls())
check = 'check'
Library = '/home-nfs/ollie/cluettig/R/x86_64-redhat-linux-gnu-library/3.3'
library('doParallel', lib.loc = Library)
library('gstat', lib.loc = Library)
library('sp', lib.loc = Library)
Library = '/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3'
library('rgdal', lib.loc = Library)
library('raster', lib.loc = Library)
library('rslurm', lib.loc = Library)
rm(Library)
setwd("/work/ollie/cluettig/github/Interpolation/Python")
source("fitmodel.R")


######### INPUT ###################################################################################################

args <- commandArgs(TRUE)
InputFileName <- args[1]
RadiusFileName <- args[2]
PriorFileName <- args[3]
OutputFileName <- args[4]
Directory <- args[5]
rm(args)
#InputFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/vx_input.tif"
#RadiusFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/kriging_radius.tif"
#Directory <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test"
#OutputFileName <- "/Users/cluettig/Documents/Daten/Greenland/data/Kriging_test/vx_int.tif"

######### FLAGS ###################################################################################################

RandomKrigingSubset = FALSE
BoundaryKrigingSubset = FALSE
VariogramSubset = TRUE

######### READ ###################################################################################################

InputRaster <- raster(InputFileName)
PriorRaster <- raster(PriorFileName)
RadiusRaster <- raster(RadiusFileName)

######### FORMAT DATA ############################################################################################
check
GivenPoints <- as.data.frame(InputRaster, xy = TRUE, na.rm = TRUE) # values and coordinates of given positions
check
ObsGivenPoints <- GivenPoints
InputLayerName <- strsplit(tail(strsplit(InputFileName,"/")[[1]], n = 1),"\\.")[[1]][1]
check
PriorPoints <- as.data.frame(PriorRaster, xy = TRUE, na.rm = TRUE)
check
names(PriorPoints)[names(PriorPoints) == "rignot_vx_rand_5000"] <- InputLayerName
check
print(names(PriorPoints))
print(names(GivenPoints))
GivenPoints <- rbind(GivenPoints, PriorPoints)
check
Radius <- as.data.frame(RadiusRaster)
rm(RadiusRaster)
AllPoints <- as.data.frame(InputRaster, xy = TRUE)

SearchedPoints <- AllPoints[(is.na(AllPoints[InputLayerName])) & Radius$kriging_radius != 0,]
SearchedPoints[InputLayerName] <- NULL 
SearchedPoints <- cbind(SearchedPoints, "z" = Radius[(is.na(AllPoints[InputLayerName])) & Radius$kriging_radius != 0,]) 
rm(Radius)
 #Coordinates and kriging radius of searched positions

if (BoundaryKrigingSubset || RandomKrigingSubset) {
    # Possibility to take only boundary points for kriging:
  if (BoundaryKrigingSubset) {
    BoundaryRaster <- boundaries(clump(InputRaster))
    HelpKrigingSample <- as.data.frame(BoundaryRaster)
    rm(BoundaryRaster)
    HelpKrigingSample[is.na(HelpKrigingSample)] <- FALSE
    KrigingSample <- AllPoints[HelpKrigingSample == 1,]
    rm(HelpKrigingSample)
  } else {
    KrigingSample <- GivenPoints
  }
  rm(BoundaryKrigingSubset)
    # Possibility to take only a subsample for kriging:
  if (RandomKrigingSubset) {
    SampleSize <- 100
    if (SampleSize < nrow(KrigingSample)) {
      KrigingSample <- KrigingSample[sample(nrow(KrigingSample), SampleSize),]
    }
  }
} else {
    KrigingSample <- GivenPoints
}
rm(InputRaster, AllPoints, RandomKrigingSubset)

# Possibility to take only a subsample for variogram estimation:
if (VariogramSubset) {
  SampleSize <- 1000
  if (SampleSize < nrow(GivenPoints)) {
    VariogramSample <- GivenPoints[sample(nrow(GivenPoints), SampleSize),]
  } else {
    VariogramSample <- GivenPoints
  }
} else {
  VariogramSample <- GivenPoints
}
rm(SampleSize, VariogramSubset)
print(nrow(VariogramSample))
print(nrow(KrigingSample))
names(VariogramSample)[names(VariogramSample) == InputLayerName] <- "z"
names(GivenPoints)[names(GivenPoints) == InputLayerName] <- "z"
names(KrigingSample)[names(KrigingSample) == InputLayerName] <- "z"

#rm(RadiusFileName, InputFileName, InputLayerName)


################ TREND REMOVAL ###########################################################################################
 print("trend removal")
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
  print(nrow(VariogramSample))
  print(nrow(KrigingSample))
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
  print(nrow(VariogramSample))
  print(nrow(KrigingSample))
# Plot of variables after removing the trend:
png(paste(Directory, "/correlation_trend_removal.png", sep = ""))
plot(VariogramSample)
dev.off()

########## VARIOGRAM ESTIMATION ###################################################################################################
coordinates(VariogramSample) <- ~x + y
# d <- data.frame(x,y,z)
# coordinates(d) = ~x+y
check
G <- gstat(formula = z~x+y,data = VariogramSample)
MeasuredVariogram <- variogram(G, width = 800, cutoff = 200000)
png(paste(Directory, "/variogram.png", sep = ""))
plot(MeasuredVariogram)
dev.off()
print(max(MeasuredVariogram$dist)/3)
EstimatedVariogram <- vgm(psill = max(MeasuredVariogram$gamma)/2, "Sph", range = max(MeasuredVariogram$dist)/3) # ,anis = c(45,0.3))
FittedVariogram <- fit.variogram(MeasuredVariogram, EstimatedVariogram)
#Plot of the fitted variogram function:
png(paste(Directory, "/variogram_fit.png", sep = ""))
plot(MeasuredVariogram, FittedVariogram)
dev.off()
rm(MeasuredVariogram, EstimatedVariogram, G, VariogramSample, Directory)
check
########## KRIGING #############################################################################################################
coordinates(KrigingSample) <- ~x+y
check
kriging <- function(i){
    if (unique(SearchedPoints$z)[i] > 100000) {max_kriging_points = 2000}else{max_kriging_points = 100}
    ActualRadius <- unique(SearchedPoints$z)[i]
    ActualPoints <- SearchedPoints[SearchedPoints$z == ActualRadius,]
    print(nrow(ActualPoints))
    print(max_kriging_points)
    coordinates(ActualPoints) <- ~x + y
    KrigingValues <- krige(formula = z~1, locations = KrigingSample, newdata = ActualPoints, model = FittedVariogram, maxdist = ActualRadius, nmax = max_kriging_points)
    #rm(FittedVariogram)
    ActualPoints@data$z <- KrigingValues@data$var1.pred
    ActualPoints <- data.frame(ActualPoints)[,1:3]
    print("check")
    ActualPoints
}
check
#print("number of cores:")
NumberOfCores <- 100
check
#print(NumberOfCores)
#Cluster <- makeCluster(NumberOfCores)
#registerDoParallel(Cluster)
MaxGapLength <- length(SearchedPoints$z)/NumberOfCores
rm(NumberOfCores)

check
ls()
u <- unique(SearchedPoints$z)
if (MaxGapLength > 1) {
  # divide gaps with to many points in parts with a maximal length of all searched points divided by number of cores
  for (j in u) {
    if (nrow(SearchedPoints[SearchedPoints$z == j,]) > MaxGapLength) {
      print("maxgaplength überschritten")
      print(nrow(SearchedPoints[SearchedPoints$z == j,]))
      print(MaxGapLength)
      NumberOfParts <- floor(nrow(SearchedPoints[SearchedPoints$z == j,])/MaxGapLength)    
      print(NumberOfParts)
      NewGroups <- list(numeric(nrow(SearchedPoints[SearchedPoints$z == j,]))+j)
      for (k in 1:NumberOfParts) {  
        NewGroups[[1]][((k - 1) * MaxGapLength + 1):(k * MaxGapLength )] <- SearchedPoints[SearchedPoints$z == j,]$z[((k - 1) * MaxGapLength + 1):(k*MaxGapLength)] + k * 0.01
      }
      SearchedPoints[SearchedPoints$z == j,][,3] <- as.data.frame(NewGroups)
    }
  }
}
rm(MaxGapLength)
check
NumberOfIterations <- length(unique(SearchedPoints$z))
print("number of iterarions:")
print(NumberOfIterations)

check
 sopt <- list(time = '12:00:00', mem = '5000M')
 check
save.image("image.RData")
sjob <- slurm_apply(kriging, data.frame(i = 1:NumberOfIterations), jobname = 'kriging', nodes = NumberOfIterations, cpus_per_node = 1, submit = FALSE, slurm_options = sopt, add_objects = c('SearchedPoints', 'FittedVariogram', 'KrigingSample'))

