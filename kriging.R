#!/usr/bin/env Rscript
rm(list = ls())
library("foreach",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-redhat-linux-gnu-library/3.3")
library("doParallel",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("gstat",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("sp",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("nlme",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("ncdf4",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("ncdf4.helpers",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("lattice",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("raster",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
library("plyr",lib.loc="/home-nfs/ollie/cluettig/R/x86_64-pc-linux-gnu-library/3.3")
source("fitmodel.R")

######### INPUT ###################################################################################################
args <- commandArgs(TRUE)
fname <- args[1]
fname2 <- args[2]
fnewname <- args[3]
dir <- args[4]
#setwd("~/Documents/Daten/Interpolation/test4_Kombination")
######### READ ###################################################################################################
nc_data <- nc_open(fname)       # input
z_data <- ncvar_get(nc_data,varid = "z")
x_coor <- ncvar_get(nc_data,varid = "x")
y_coor <- ncvar_get(nc_data,varid = "y")
nc_close(nc_data)

nc_data <- nc_open(fname2)      # contains polygon diameters
krig_data <- ncvar_get(nc_data,varid = "Band1")
nc_close(nc_data)
######### FORMAT DATA ############################################################################################
# raster to vectors:
z_a <- as.vector(z_data)
y_a <- rep(x_coor,times = length(y_coor))
x_a <- rep(y_coor,each = length(x_coor))
krig_a <- as.vector(krig_data)

# NA positions
z_a[z_a > 1e+30] <- NA    # convert to NaN
krig_nan <- krig_a[krig_a != 0 & is.na(z_a)]
x_nan <- x_a[krig_a != 0 & is.na(z_a)]
y_nan <- y_a[krig_a != 0 & is.na(z_a)]
z_nan <- z_a[krig_a != 0 & is.na(z_a)]


# Values in given positions:
x_notnan <- x_a[!is.na(z_a)]
y_notnan <- y_a[!is.na(z_a)]
z_notnan <- z_a[!is.na(z_a)]

# Possibility to take only a subsample for kriging:
sub1 = length(z_notnan)/1 # number of data points in subsample
if (length(z_notnan) > sub1)
{
  s <- 1:length(z_notnan)
  IX <- sample(s,sub1)
  x <- x_notnan[IX]
  y <- y_notnan[IX]
  z <- z_notnan[IX]
} else
{
  x <- x_notnan
  y <- y_notnan
  z <- z_notnan
}

# Possibility to take only outer points for kriging:
outer_p <- 1
print("reaches outer")
if (outer_p == 1)
{
  z_outer <- z_data
  nan_mat <- is.na(z_data)
  for (i in 1:length(x_coor))
  {
      for (j in 1:length(y_coor))
      {
	if (nan_mat[i,j]==0)
	{
	  if (i-1 >= 1 && i+1 <= length(x_coor) && j-1 >= 1 && j+1 <= length(y_coor))
	  {
		n_nan_neighbours <- nan_mat[i-1,j] + nan_mat[i+1,j] + nan_mat[i,j-1] + nan_mat[i,j+1]
		if (n_nan_neighbours == 0)
		{
			z_outer[i,j] <- NA
		}
	  }
	}
      }
  }  
  z_o <- as.vector(z_outer)
  y_o <- rep(x_coor,times = length(y_coor))
  x_o <- rep(y_coor,each = length(x_coor))
  x <- x_o[!is.na(z_o)]
  y <- y_o[!is.na(z_o)]
  z <- z_o[!is.na(z_o)]
} else
{
  x <- x_notnan
  y <- y_notnan
  z <- z_notnan 
}

# Possibility to take only a subsample for variogram estimation:
sub2 = 1000 # number of data points in subsample
if (length(z) > sub2)
{
  s <- 1:length(z)
  IX <- sample(s,sub2)
  x_v <- x[IX]
  y_v <- y[IX]
  z_v <- z[IX]
} else {
  x_v <- x
  y_v <- y
  z_v <- z
}
rm(x_o,y_o,z_o,z_outer,nan_mat,z_data,s,IX)

################ TREND REMOVAL ###########################################################################################
print("trend removal")
model <- fitmodel(x_v,y_v,z_v, dir)
model$coefficients[is.na(model$coefficients)] <- 0
if (length(model$coefficients) == 16)
{
  print("cubic trend model")
  z_v <- z_v - (model$coefficients[1] + model$coefficients[2]*x_v + model$coefficients[3]*x_v^2 
               + model$coefficients[4]*x_v^3 + model$coefficients[5]*y_v + model$coefficients[6]*y_v^2 
               + model$coefficients[7]*y_v^3 + model$coefficients[8]*x_v*y_v + model$coefficients[9]*x_v^2*y_v 
               + model$coefficients[10]*x_v^3*y_v + model$coefficients[11]*x_v*y_v^2 + model$coefficients[12]*x_v^2*y_v^2 
               + model$coefficients[13]*x_v^3*y_v^2 + model$coefficients[14]*x_v*y_v^3 + model$coefficients[15]*x_v^2*y_v^3 
               + model$coefficients[16]*x_v^3*y_v^3)
  z <- z - (model$coefficients[1] + model$coefficients[2]*x + model$coefficients[3]*x^2 
            + model$coefficients[4]*x^3 + model$coefficients[5]*y + model$coefficients[6]*y^2 
            + model$coefficients[7]*y^3 + model$coefficients[8]*x*y + model$coefficients[9]*x^2*y 
            + model$coefficients[10]*x^3*y + model$coefficients[11]*x*y^2 + model$coefficients[12]*x^2*y^2 
            + model$coefficients[13]*x^3*y^2 + model$coefficients[14]*x*y^3 + model$coefficients[15]*x^2*y^3 
            + model$coefficients[16]*x^3*y^3)
}else if (length(model$coefficients) == 9)
{
  print("quadratic trend model")
  z_v <- z_v - (model$coefficients[1] + model$coefficients[2]*x_v + model$coefficients[3]*x_v^2 
               + model$coefficients[4]*y_v + model$coefficients[5]*y_v^2 + model$coefficients[6]*x_v*y_v 
               + model$coefficients[7]*x_v^2*y_v + model$coefficients[8]*x_v*y_v^2 
               + model$coefficients[9]*x_v^2*y_v^2)
  z <- z - (model$coefficients[1] + model$coefficients[2]*x + model$coefficients[3]*x^2 
            + model$coefficients[4]*y + model$coefficients[5]*y^2 + model$coefficients[6]*x*y 
            + model$coefficients[7]*x^2*y + model$coefficients[8]*x*y^2 + model$coefficients[9]*x^2*y^2)
}else if (length(model$coefficients) == 4)
{
  print("linear trend model")
  z_v <- z_v - (model$coefficients[1] + model$coefficients[2]*x_v 
               + model$coefficients[3]*y_v + model$coefficients[4]*x_v*y_v)
  z <- z - (model$coefficients[1] + model$coefficients[2]*x 
            + model$coefficients[3]*y + model$coefficients[4]*x*y)
}else{
  print("Error while substracting the trend model")
}
d_v <- data.frame(x_v,y_v,z_v)
# Plot of variables after removing the trend:
png(paste(dir,"corr_model_substr.png"))
plot(d_v)
dev.off()

########## VARIOGRAM ESTIMATION ###################################################################################################
print("variogram")
coordinates(d_v) = ~x_v+y_v
d <- data.frame(x,y,z)
coordinates(d) = ~x+y
g <- gstat(id = "test", formula = z_v~1,data = d_v)
v <- variogram(g, alpha = c(0,45,90,135),width = 10,cutoff = 10000000)
v2 <- variogram(g, width = 100,cutoff = 10000000)
# Plot of the variogram:
png(paste(dir,"variogram.png"))
plot(v2)
dev.off()
vm <- vgm(psill = max(v2$gamma)/2,"Sph",range = max(v2$dist)/3)# ,anis = c(45,0.3))
fit <- fit.variogram(v,vm)
#Plot of the fitted variogram function:
png(paste(dir,"variogram_fit.png"))
plot(v, fit)
dev.off()
rm(x_v,y_v,z_v,d_v,v,v2,g,v,v2,vm)

t3 <- Sys.time()

########## KRIGING #############################################################################################################
print("kriging")
x_notnan <- x_a[krig_a == 0 | (krig_a != 0 & !is.na(z_a))]
y_notnan <- y_a[krig_a == 0 | (krig_a != 0 & !is.na(z_a))]
z_notnan <- z_a[krig_a == 0 | (krig_a != 0 & !is.na(z_a))]
d_notnan <- data.frame(x_notnan,y_notnan,z_notnan)
d_notnan <- rename(d_notnan, c("x_notnan" = "x", "y_notnan" = "y","z_notnan" = "z"))
rm(x_notnan,y_notnan,z_notnan,x_a,y_a,z_a)
kriging <- function(i){
    print("kriging run")
    print(i)
    if (unique(krig_nan)[i] > 25000) {max_kriging_points = 3000}else{max_kriging_points = 500}
    x_nan_i <- x_nan[krig_nan == unique(krig_nan)[i]]
    y_nan_i <- y_nan[krig_nan == unique(krig_nan)[i]]
    d_n <- data.frame(x_nan_i,y_nan_i)
    coordinates(d_n) = ~x_nan_i+y_nan_i
    cat(i, 'of', n_iterations, " at time", as.character(Sys.time()), "for",length(x_nan_i),"points and with maximal", max_kriging_points, "kriging points \n")
    zk <- krige(formula = z~1, locations = d, newdata = d_n, model = fit, maxdist = 5*unique(krig_nan)[i], nmax=1000)
    #maxdist = unique(krig_nan)[i]
    zw <- zk$var1.pred
    xw <- x_nan_i
    yw <- y_nan_i
    dw <- data.frame(xw,yw,zw)
    dw <- rename(dw, c("xw"="x", "yw"="y", "zw"="z"))
}
print("number of cores:")
n_cores <- detectCores()
print(n_cores)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
max_length_gap=length(krig_nan/n_cores)
for (j in unique(krig_nan))
{
if (length(krig_nan[krig_nan==j])>max_length_gap)
   {
   part=max_length_gap/length(krig_nan[krig_nan==j])
   for (k in 1:floor(part))
       {  
		krig_nan[krig_nan==j][(k-1)*max_length_gap+1:k*max_length_gap]=krig_nan[(k-1)*max_length_gap+1:k*max_length_gap]+k*0.01
       }
   }
} 

n_iterations <- length(unique(krig_nan))
print("number of iterarions:")
print(n_iterations)
kriging_results <- foreach(test_i = 1:n_iterations,
                           .combine = list,
                            .multicombine = TRUE,
                            .packages = c("gstat", "sp", "nlme", "ncdf4", "lattice", "plyr"))  %dopar%
    kriging(test_i)

stopImplicitCluster()
d_result <- kriging_results[[1]]
for (i in 2:length(kriging_results)) {
    d_result <- rbind(d_result,kriging_results[[i]])
}
z_n <- d_result$z
x_n <- d_result$x
y_n <- d_result$y
t4 <- Sys.time()
cat("Kriging runtime:", t4 - t3)
# Addition of substracted trend model:
if (length(model$coefficients) == 16)
{
   z_n <- z_n + (model$coefficients[1] + model$coefficients[2]*x_n + model$coefficients[3]*x_n^2
             + model$coefficients[4]*x_n^3 + model$coefficients[5]*y_n + model$coefficients[6]*y_n^2
             + model$coefficients[7]*y_n^3 + model$coefficients[8]*x_n*y_n + model$coefficients[9]*x_n^2*y_n
             + model$coefficients[10]*x_n^3*y_n + model$coefficients[11]*x_n*y_n^2 + model$coefficients[12]*x_n^2*y_n^2
             + model$coefficients[13]*x_n^3*y_n^2 + model$coefficients[14]*x_n*y_n^3 + model$coefficients[15]*x_n^2*y_n^3
             + model$coefficients[16]*x_n^3*y_n^3)
}else if (length(model$coefficients) == 9)
{
  z_n <- z_n + (model$coefficients[1] + model$coefficients[2]*x_n + model$coefficients[3]*x_n^2
             + model$coefficients[4]*y_n + model$coefficients[5]*y_n^2 + model$coefficients[6]*x_n*y_n
             + model$coefficients[7]*x_n^2*y_n + model$coefficients[8]*x_n*y_n^2 + model$coefficients[9]*x_n^2*y_n^2)
}else if (length(model$coefficients) == 4)
{
   z_n <- z_n + (model$coefficients[1] + model$coefficients[2]*x_n
             + model$coefficients[3]*y_n + model$coefficients[4]*x_n*y_n)
}else{
   print("Error while addition of the model")
}

############# WRITE ##################################################################################################
print("writing")
# kriged values:
d_z_n <- data.frame(x_n,y_n,z_n)
d_z_n <- rename(d_z_n, c("x_n" = "x", "y_n" = "y", "z_n" = "z"))
total <- rbind(d_z_n,d_notnan)
total <- total[with(total, order(x, y)),]
z_new <- as.matrix(total[c("z")], nrow = length(x_coor), ncol = length(y_coor))
nc_data <- nc_open(fnewname, write = TRUE)
z_data <- ncvar_put(nc_data, varid = "z",vals = z_new)
nc_close(nc_data)