#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
build_model <- function(dat_pt = "./indir/Spectra/CrownBrdfSpectra.csv"
                        ,loop=1,  nrmlz = F, trait = c("N_pct", "C_pct", "P_pct", "LMA_g.m2")){
  library(tidyverse)
  library(plsRglm)
  source("./R/clean_spectra.R")
  source("./R/clean_permutations.R")
  source("./R/get_random_bags.R")
  source("./R/cut_set.R")
  source("./R/pls_glm.R")

  #clean pixels in the dataset using ndvi and nir fitler, and maybe  detecting outliers from pca
  spectra = read_csv(dat_pt)
  reduced_spectra = clean_spectra(spectra, ndvi = 0.5, nir = 0.2, outlier = F)
  spectra = cbind.data.frame(spectra[reduced_spectra$good_pix, 1:2], reduced_spectra$refl)

  #are random combinations of pixels already set up?
  combinations <- file.exists(paste('./indir/Permutations/onePix1Crown_', loop, ".csv", sep=""))
  if(!combinations){
    #extract n combinations of pixles by extracting one per bag
    get_random_bags(spectra, lp = loop)
  }
  # run the pls glm on training bags for each random extractions
  for(tr in trait){
    random_bag_pls <- pls_glm(trait = tr,
                              ll = loop, nrmlz = nrmlz)
  }
}
