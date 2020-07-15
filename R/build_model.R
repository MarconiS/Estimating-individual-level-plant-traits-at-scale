#!/bin/bash
build_model <- function(loop=1, dat_pt = "./indir/Spectra/silva_final_15.csv"
                        ,tr = c("LMA", "Npercent", "Ppercent", "Cpercent")){
  library(tidyverse)
  library(plsRglm)
  source("./R/clean_spectra.R")
  source("./R/clean_permutations.R")
  source("./R/get_random_bags.R")
  source("./R/cut_set.R")
  source("./R/pls_glm.R")
  #rand_range = readr::read_table("random_example.txt", col_names=T)
  #rand_range = colnames(rand_range)
  #loop = as.integer(rand_range[loop])
  wrangle = function(x){
    ifelse(is.character(x), x[1], mean(x))
  }
  #clean pixels in the dataset using ndvi and nir fitler, and maybe  detecting outliers from pca
  spectra = read_csv(dat_pt)
  spectra = spectra[complete.cases(spectra),]
  reduced_spectra = clean_spectra(spectra,  ndvi = 0.7, nir = 0.3)
  spectra = cbind.data.frame(spectra[reduced_spectra$good_pix, 1:2], reduced_spectra$refl)
  spectra = spectra %>% filter(!individualID %in% c(319, 211, 367, 347, 2111, 2234))
  spectra = spectra %>% group_by(individualID)%>% top_n(1, wt = flpt)
  spectra_ave = spectra %>% group_by(individualID) %>% summarize_all(wrangle)
  spectra = spectra %>% dplyr::select(-one_of("flpt"))
  #spectra = readr::read_csv("./indir/Spectra/plot_reflectance.csv")
  #readr::write_csv(spectra, "./indir/Spectra/reflectance_silva.csv")
  #readr::write_csv(spectra_ave, "./indir/Spectra/plot_save.csv")
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
