# make prediction for a tile

# get spectra and clean it
f = "/Users/sergiomarconi/Documents/Data/bdrf/141_2.tif"
siteID = "TALL"
sites = c("OSBS", "TALL")
dat = raster::brick(f)
dat = as.data.table.raster(dat)
colnames(dat) = paste("band", 1:369, sep="_")
reduced_spectra = clean_spectra(dat, ndvi = 0.5, nir = 0.2)
dat = cbind.data.frame(dat[reduced_spectra$good_pix, 1:2], reduced_spectra$refl)
#add site band
bnd_site <- rep(siteID, nrow(dat)) %>% factor(levels=sites) %>% fastDummies::dummy_cols()
colnames(bnd_site) <- stringr::str_replace(colnames(bnd_site), ".data_", "band_")
dat = cbind.data.frame(bnd_site, dat)
model_stack <- readRDS(paste("./outdir/EPBMs/100", trait, ".rds", sep=""))

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

for(bb in 1:length(weights)){
  md = selected_mods[bb]
  pls.mod.train <- model_stack[[md]]$mod
  optim.ncomps <- model_stack[[md]]$ncomp
  #make predictions using the ith model
  ith_mod_prediction <- pls_glm_predict(pls.mod.train, newdata = test.PLS,
                                        wt = rep(1, nrow(test.PLS)),
                                        ncomp=optim.ncomps,  type='response')
  ith_mod_prediction=(ith_mod_prediction)
  output.daic$yhat <- output.daic$yhat + ith_mod_prediction[,1] * weights[bb]
  output.up.daic$yhat <- output.up.daic$yhat + ith_mod_prediction[,3] * weights[bb]
  output.lw.daic$yhat <- output.lw.daic$yhat + ith_mod_prediction[,2] * weights[bb]
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind.data.frame(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
  }else{
    output <- rbind.data.frame(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]),
                                                       as.vector(pred.val.data$fit))))
  }
}
