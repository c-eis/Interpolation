fitmodel <- function(x,y,z,d){
  #png(paste(d,"corr.png"))
  #plot(data.frame(x,y,z))
  #dev.off()
  model_lin_xy <- lm(z~poly(x,degree=1,raw=TRUE)*poly(y,degree=1,raw=TRUE))
  r2_lin_xy <- summary(model_lin_xy)$adj.r.squared
  model_quad_xy <- lm(z~poly(x,degree=2,raw=TRUE)*poly(y,degree=2,raw=TRUE))
  r2_quad_xy <- summary(model_quad_xy)$adj.r.squared
  model_cub_xy <- lm(z~poly(x,degree=3,raw=TRUE)*poly(y,degree=3,raw=TRUE))
  r2_cub_xy <- summary(model_cub_xy)$adj.r.squared
  index <- which.max(c(r2_lin_xy,r2_quad_xy,r2_cub_xy))
  if (index == 1){
    return(model_lin_xy)
  }
  else if (index==2){
    return(model_quad_xy)
  }
  else if (index==3){
    return(model_cub_xy)
  }
  else {
    print("Error while computing the best fitting model")
  }
}