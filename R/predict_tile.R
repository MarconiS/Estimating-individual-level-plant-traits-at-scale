# make prediction for a tile
library(tidyverse)
library(raster)
source("./R/fasRas_to_dt.R")
source("./R/clean_spectra.R")
source("")
#retrieve number of snaps
softmax <- function(x) {
  x <- x[!is.na(x)]

  exp(x) / sum(exp(x))
}


# get spectra and clean it
f = "394000_3283000_152940.tif"
siteID = "OSBS"
trait = "LMA"
sites = c("OSBS", "TALL")
pt = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/corrHSI/"
outdir = "//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/traits/"
dat = raster::brick(paste(pt,f, sep="/"))
nbags = 100
rbbox = dim(dat)
dat = as.data.table.raster(dat)
colnames(dat) = paste("band", 1:369, sep="_")
reduced_spectra = clean_spectra(dat, ndvi = 0.5, nir = 0.2)
dat = reduced_spectra$refl
reduced_spectra = reduced_spectra$good_pix
#dat = dat[complete.cases(dat),]
#add site band
bnd_site <- rep(siteID, nrow(dat)) %>% factor(levels=sites) %>% fastDummies::dummy_cols()
colnames(bnd_site) <- stringr::str_replace(colnames(bnd_site), ".data_", "band_")
dat = cbind.data.frame(bnd_site, dat)
rm(bnd_site)
model_stack <- readRDS(paste("/orange/ewhite/s.marconi/Chapter1/hiPyRneon/outdir/EPBMs/1002", trait, ".rds", sep=""))
mod_col = model_stack[[1]]$mod$dataX %>% names
dat = dat %>% dplyr::select(mod_col)
# calculate, scale the dAIC to rank and weight each model using a softmax function
mod.aic=rep(0,length(model_stack))
for(bb in 1: length(model_stack)){
  mod.aic[bb] <- model_stack[[bb]]$mod$FinalModel$aic
}
selected_mods = which(mod.aic %in% sort(mod.aic, decreasing = F)[1:nbags])
#mod.aic = mod.r2
mod.aic <- scale(mod.aic[selected_mods])
delta.aic <- mod.aic - min(mod.aic)
weights <- softmax(-0.5*delta.aic)

output.daic = matrix(0,nrow(dat), 3)
for(bb in 1:length(weights)){
  md = selected_mods[bb]
  pls.mod.train <- model_stack[[md]]$mod
  optim.ncomps <- model_stack[[md]]$ncomp
  #make predictions using the ith model
  # ith_mod_prediction <- pls_glm_predict(pls.mod.train, newdata = dat,
  #                                       wt = rep(1, nrow(test.PLS)),
  #                                       ncomp=optim.ncomps,  type='response')
  ith_mod_prediction <- plsglm_fast_predict(pls.mod.train, newdata = dat,
                                        wt = rep(1, nrow(test.PLS)),
                                        ncomp=optim.ncomps,  type='response')
  output.daic <-ith_mod_prediction * weights[bb]
}

dat = matrix(NA, length(reduced_spectra), 3)
dat[reduced_spectra,] = output.daic
rm(output.daic)

dim(dat) = c(rbbox[1:2],3)
lyr = raster(f)
dat = raster::brick(dat, xmn=lyr@extent[1], xmx=lyr@extent[2], #nl = 9,
                    ymn=lyr@extent[3], ymx=lyr@extent[4], crs=lyr@crs, transpose=FALSE)
names(dat) = paste(trait, c("hat", "lw","up"), sep="_")
brick <- raster::stack(brick, lyr)
raster::writeRaster(brick, paste(outdir, trait, f, sep ="/"))
