#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#' @import plsRglm
#'

# Decorrelation Stretching raster images in R
#  based on https://gist.github.com/fickse/82faf625242f6843249774f1545d7958
decorrelation_stretch <- function(pathr, outdir="~/Documents/GitHub/NeonTreeEvaluation/evaluationDecor"){
  #load raster
  tryCatch({
    r <- raster::brick(pathr)#[[c(13,55,134)]]

    #get plot name
    fname <- substring(pathr, nchar(pathr)-12+1, nchar(pathr)-4)
    # r must be a >= 3 band raster
    # determine eigenspace
    means_per_layer <- lapply(1:dim(r)[3], function(x) fill_gaps(r[[x]])) #
    r <- do.call(raster::brick, means_per_layer)
    pc <- princomp(r[])
    # get inverse rotation matrix
    R0 <- solve(pc$loadings)
    # 'stretch' values in pc space, then transform back to RGB space
    fun <- function(x){(x-min(x))/(max(x)-min(x))*255}
    scp  <- apply(predict(pc), 2, function(x) scale(ecdf(x)(x), scale = FALSE))
    scpt <- scp %*% R0
    r[] <- apply(scpt, 2, fun)
    raster::writeRaster(r, filename = paste(outdir, "/", fname, ".tif", sep=""),
                        datatype = 'INT2U',
                        overwrite=TRUE)
  })
}
# example
# b <- brick(system.file("external/rlogo.grd", package="raster"))
# dc <- decorrelation_stretch(b)
# plotRGB(dc)


fill_gaps <- function(r){
  mean_r <- raster::cellStats(r, 'mean',  na.rm=TRUE)
  raster::values(r)[is.na(raster::values(r))] = mean_r
  return(r)
}


decorrelation_stretch <- function(pathr, outdir="~/Documents/GitHub/NeonTreeEvaluation/evaluationDecor"){
  #load raster
  r <- raster::brick(pathr)#[[c(13,55,134)]]

  #get plot name
  fname <- substring(pathr, nchar(pathr)-12+1, nchar(pathr)-4)
  # r must be a >= 3 band raster
  # determine eigenspace
  means_per_layer <- lapply(1:dim(r)[3], function(x) fill_gaps(r[[x]])) #
  r <- do.call(raster::brick, means_per_layer)


  raster::writeRaster(r, filename = paste(outdir, "/", fname, ".tif", sep=""),
                      datatype = 'INT2U',
                      overwrite=TRUE)
}
