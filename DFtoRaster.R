DFtoRaster <- function(x,xcol=1,ycol=2,d=3) {
  x <- x[,c(xcol,ycol,d)]
  a <- colnames(x)
  names(x) <- c('x','y',a[d])
  coordinates(x) <- ~ x+y
  gridded(x) <- TRUE
  rasterDF <- raster(x)
  return(rasterDF)
}


