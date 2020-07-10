#!/bin/bash

args <- commandArgs()
print(args)
f = args[6]
siteID = args[7]
trait = args[8]
nbags = args[9]
get_transformation = F
# make prediction for a tile
library(tidyverse)
library(raster)
source("./R/fasRas_to_dt.R")
source("./R/clean_spectra.R")

#retrieve number of snaps
softmax <- function(x) {
  x <- x[!is.na(x)]

  exp(x) / sum(exp(x))
}


# get spectra and clean it
# f = "465000_3637000_145036.tif"
# siteID = "TALL"
# trait = "Ppercent"
# nbags = 100

#epsg = 32617
sites = c("OSBS", "TALL")
if(siteID == "OSBS"){
  pt = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/corrHSI/"
  outdir = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/traits/"
  tmp_dir = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/tmp/"
  epsg = 32617
}else{
  pt = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D08/TALL/L4/corrHSI/"
  outdir = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D08/TALL/L4/traits/"
  tmp_dir = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D08/TALL/L4/tmp/"
  epsg = 32616
}
dat = raster::brick(paste(pt,f, sep="/"))
rbbox = dim(dat)
dat = as.data.table.raster(dat)
colnames(dat) = paste("band", 1:369, sep="_")
if(get_transformation == T){
  dat
  reduced_spectra = clean_spectra(dat, ndvi = 0.7, nir = 0.3)
  saveRDS(reduced_spectra, paste(tmp_dir, f, sep ="/"))
}else{
  reduced_spectra = readRDS(paste(tmp_dir, f, sep ="/"))
  dat = reduced_spectra$refl
  reduced_spectra = reduced_spectra$good_pix
  if( nrow(dat) !=0){
    #dat = dat[complete.cases(dat),]
    #add site band
    #check = apply(dat, 1, min)
    bnd_site <- rep(siteID, nrow(dat)) %>% factor(levels=sites) %>% fastDummies::dummy_cols()
    colnames(bnd_site) <- stringr::str_replace(colnames(bnd_site), ".data_", "band_")
    dat = cbind.data.frame(bnd_site, dat)
    rm(bnd_site)
    model_stack <- readRDS(paste("/orange/ewhite/s.marconi/Chapter1/hiPyRneon/outdir/EPBMs/100_", trait, ".rds", sep=""))
    mod_col = model_stack[[1]]$mod$dataX %>% names
    dat = dat %>% dplyr::select(mod_col)
    # calculate, scale the dAIC to rank and weight each model using a softmax function
    mod.aic=rep(0,length(model_stack))
    for(bb in 1: length(model_stack)){
      mod.aic[bb] <- model_stack[[bb]]$mod$FinalModel$aic
    }
    #selected_mods = which(mod.aic %in% sort(mod.aic, decreasing = F)[1:nbags])
    #mod.aic = mod.r2
    mod.aic <- scale(mod.aic)#[selected_mods])
    delta.aic <- mod.aic - min(mod.aic)
    weights <- softmax(-0.5*delta.aic)

    output.daic = matrix(0,nrow(dat), 3)
    for(bb in 1:length(weights)){
      md = bb #selected_mods[bb]
      pls.mod.train <- model_stack[[md]]$mod
      optim.ncomps <- model_stack[[md]]$ncomp
      #make predictions using the ith model
      newdata = dat
      newdata <- sweep(sweep(newdata, 2, attr(pls.mod.train$ExpliX, "scaled:center")),
                       2, attr(pls.mod.train$ExpliX, "scaled:scale"), "/")
      newdata <- as.matrix(newdata)
      nrnd <- nrow(newdata)
      newdata = lapply(1:nrnd, function(x)c(newdata[x,] %*% pls.mod.train$wwetoile[, 1:optim.ncomps],
                                            rep(0, pls.mod.train$computed_nt - optim.ncomps)))
      newdata = do.call(rbind, newdata)
      colnames(newdata) <- NULL
      newdata <- data.frame(tt = newdata)
      pred_int = HH::interval(pls.mod.train$FinalModel, newdata=newdata, type="response")
      pred_int = pred_int[,c(1,4,5)]
      colnames(pred_int) = c("fit","lwr","upr")
      output.daic <- output.daic + pred_int * weights[bb]
    }
    #rm(model_stack)
    dat = matrix(NA, length(reduced_spectra), 3)
    dat[reduced_spectra,] = exp(output.daic)
    rm(output.daic)

    dim(dat) = c(rbbox[2:1],3)
    lyr = raster(paste(pt,f, sep="/"))
    dat = raster::brick(dat, xmn=lyr@extent[1], xmx=lyr@extent[2], #nl = 9,
                        ymn=lyr@extent[3], ymx=lyr@extent[4], crs=paste('+init=epsg:', epsg, sep=""), transpose=TRUE)
    names(dat) = paste(trait, c("hat", "lw","up"), sep="_")
    writeRaster(dat, paste(outdir, trait, f, sep ="/"), overwrite = T)
  }
}
