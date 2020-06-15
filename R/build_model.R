#!/bin/bash
build_model <- function(loop=1, dat_pt = "./indir/Spectra/CrownBrdfSpectra.csv"
                        , nrmlz = F,
                        trait = c("Npercent", "LMA", "Ppercent", "Cpercent")){
  library(tidyverse)
  library(plsRglm)
  source("./R/clean_spectra.R")
  source("./R/clean_permutations.R")
  source("./R/get_random_bags.R")
  source("./R/cut_set.R")
  source("./R/pls_glm.R")
  wrangle = function(x){
    ifelse(is.character(x), x[1], mean(x))
  }
  #clean pixels in the dataset using ndvi and nir fitler, and maybe  detecting outliers from pca
  spectra = read_csv(dat_pt)
  reduced_spectra = clean_spectra(spectra, ndvi = 0.5, nir = 0.2, outlier = F)
  spectra = cbind.data.frame(spectra[reduced_spectra$good_pix, 1:2], reduced_spectra$refl)
  spectra_ave = spectra %>% group_by(individualID) %>% summarize_all(wrangle)
  #are random combinations of pixels already set up?
  combinations <- file.exists(paste('./indir/Permutations/onePix1Crown_', loop, ".csv", sep=""))
  #fprint(trait, loop, nrmlz)
  #if(!combinations){
  #extract n combinations of pixles by extracting one per bag
  get_random_bags(spectra, lp = loop)
  #}
  # run the pls glm on training bags for each random extractions
  for(trait in tr){
    random_bag_pls <- pls_glm(trait = trait, ll = loop, nrmlz = nrmlz)
  }
}

2
3

args <- commandArgs()
print(args)
build_model(loop=args[6])
