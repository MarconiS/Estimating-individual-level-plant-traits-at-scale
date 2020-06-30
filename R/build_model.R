#!/bin/bash
build_model <- function(loop=1, dat_pt = "./indir/Spectra/CrownBrdfSpectra.csv"
                        #, nrmlz = F,
                        ,tr = c("LMA", "Npercent", "Ppercent", "Cpercent")){
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
  reduced_spectra = clean_spectra(spectra, ndvi = 0.4, nir = 0.2, outlier = F)
  spectra = cbind.data.frame(spectra[reduced_spectra$good_pix, 1:2], reduced_spectra$refl)
  spectra_ave = spectra %>% group_by(individualID) %>% summarize_all(wrangle)
  #all_spectra = readr::write_csv(spectra, "./indir/Spectra/reflectance_all.csv")
  #extract n combinations of pixles by extracting one per bag
  get_random_bags(spectra, lp = loop)
  # run the pls glm on training bags for each random extractions
  for(trait in tr){
    random_bag_pls <- pls_glm(trait = trait, ll = loop, nrmlz = nrmlz)
  }
}


args <- commandArgs()
print(args)
build_model(loop=as.integer(args[6]))
