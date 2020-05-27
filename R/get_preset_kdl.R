#' r wrapper to use functions from hsi toolkit
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#' @import plsRglm
#'

get_preset_kdl <- function(path_hsi, bands = 15, outdir = "./outdir/kld", mat_args = NULL){
  fname <- substring(path_hsi, nchar(path_hsi)-12+1, nchar(path_hsi))
  hsi <- brick(path_hsi)
  kld_groups <- readr::read_csv("./indir/Misc/bands_10clusters_kld.csv") %>%
    data.frame
  bands = kld_groups %>% unique %>% dim %>% max
  kld_array <- as.array(t(hsi))
  dims_mat <- dim(hsi)
  #make matrix from array
  v <- matrix(kld_array, dims_mat[1]* dims_mat[2], dims_mat[3]) %>% data.frame
  v <- cbind.data.frame(kld_groups, t(v)) %>%
    group_by(KLD) %>% summarise_all(mean) %>%
    mutate_all(~ . /10000) %>% t
  v <- v[-1,]
  x <- hsi[[1:bands]]
  # create a rasterBrick with the reduced dimensions
  kld_ras <- setValues(x, v)
  #kld_ras <- rescale_aop(kld_ras, 0, 10000)
  #plotRGB(hsi, r = 15, g = 58, b = 104, stretch = "lin")
  #plotRGB(kld_ras, r = 1, g = 3, b = 7, stretch = "lin")

  raster::writeRaster(kld_ras, #
                      filename = paste(outdir, fname, sep="/"),
                      #datatype = 'INT2U',
                      overwrite=TRUE)

  return(fname)
}
