#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
decorrelation_bands <- function(raster_path, raster_outpath){
  list_rasters <- list.files(raster_path, recursive = F, full.names = T, pattern = ".tif")
  pblapply(list_rasters, decorrelation_stretch, outdir=raster_outpath)
}

normalize_bands <- function(raster_path, raster_outpath){
  list_rasters <- list.files(raster_path, recursive = F, full.names = T, pattern = ".tif")
  pblapply(list_rasters, normalize_stretch, outdir=raster_outpath)
}
