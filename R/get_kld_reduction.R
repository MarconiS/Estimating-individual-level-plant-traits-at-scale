#' r wrapper to use functions from hsi toolkit
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#' @import plsRglm
#'

get_kld_reduction <- function(path_hsi, bands = 15, outdir = NA, mat_args = NULL){
  fname <- substring(path_hsi, nchar(path_hsi)-12+1, nchar(path_hsi))
  hsi <- brick(path_hsi)
  source_python("../hsi_toolkit_py/dim_reduction/hdr.py")
  kld_array <- dimReduction(as.array(t(hsi)), bands)

  #make matrix from array
  v <- sapply(kld_array, function(x) as.vector(t(x)))

  x <- hsi[[1:bands]]
  # create a rasterBrick with the reduced dimensions
  kld_ras <- setValues(x, v)
  #kld_ras <- rescale_aop(kld_ras, 0, 10000)
  #plotRGB(hsi, r = 15, g = 58, b = 104, stretch = "lin")
  plotRGB(kld_ras, r = 1, g = 3, b = 7, stretch = "lin")

  if(!is.na(outdir)){
    raster::writeRaster(kld_ras, #
            filename = paste(outdir, fname, sep="/"),
            #datatype = 'INT2U',
            overwrite=TRUE)
  }
  return(fname)
}
