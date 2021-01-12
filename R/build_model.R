#!/bin/bash
build_model <- function(loop=1, dat_pt = "./indir/Spectra/field_delineated_spectra.csv"
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
  spectra = spectra[complete.cases(spectra),]
  reduced_spectra = clean_spectra(spectra,  ndvi = 0.7, nir = 0.3)
  spectra = cbind.data.frame(spectra[reduced_spectra$good_pix, 1:2], reduced_spectra$refl)
  spectra = spectra %>% group_by(individualID)%>% top_n(1, wt = flpt)
  spectra_ave = spectra %>% group_by(individualID) %>% summarize_all(wrangle)
  spectra = spectra %>% dplyr::select(-one_of("flpt"))
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
