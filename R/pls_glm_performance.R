
trait = "Ppercent"
nbags = 100
pls_glm_predict <- function(object,newdata,
                            comps=object$computed_nt, type=c("link", "response", "terms", "scores", "class", "probs"),
                            se.fit=FALSE, wt = NULL, dispersion = NULL,methodNA="adaptative",verbose=TRUE,...)
{
  nrnd <- nrow(newdata)
  newdataNA <- !is.na(newdata)
  newdata <- sweep(sweep(newdata, 2, attr(object$ExpliX, "scaled:center")),
                   2, attr(object$ExpliX, "scaled:scale"), "/")
  newdata <- as.matrix(newdata)
  #newdata[!newdataNA] <- 0
  newdata = lapply(1:nrnd, function(x)c(newdata[x,] %*% object$wwetoile[, 1:comps],
                                        rep(0, object$computed_nt - comps)))
  newdata = do.call(rbind, newdata)
  colnames(newdata) <- NULL
  newdata <- data.frame(tt = newdata)
  pred_int = HH::interval(object$FinalModel, newdata=newdata, type="response")
  pred_int = pred_int[,c(1,4,5)]
  colnames(pred_int) = c("fit","lwr","upr")
  return(pred_int[])
}

#pls_glm_performance <- function(trait = "N_pct", nbags = 10, normz=F){
#retrieve number of snaps
softmax <- function(x) {
  x <- x[!is.na(x)]
  exp(x) / sum(exp(x))
}

rmse <- function(x_p, x_o){
  ssd = sum((x_p - x_o)^2)
  rmse = sqrt(ssd/length(x_o))
}
loops <- list.files("./outdir/PBMs", pattern = trait, full.names = T)

#check if Ensemble model has been created yet
model_stack <- file.exists(paste("./outdir/EPBMs/", trait, ".rds", sep=""))

#if(!model_stack || length(model_stack) != nbags){
#get delta_aic and rank best models
mod.aic <- rep(NA, length(loops))
mod.r2 <- rep(NA, length(loops))

model_stack <- NULL
for(bb in 1: length(mod.aic)){
  foo <- readRDS(loops[bb])
  mod.aic[bb] <- foo$mod$FinalModel$aic
  mod.r2[bb] <- - foo$pR2
  #create model stack with all random bags models
  #model_stack[[bb]] <-foo
}
#subset the best n models to be used for performance
mask <- cbind(which(mod.r2 %in% head(sort(mod.r2), nbags)),
              mod.r2[mod.r2 %in% head(sort(mod.r2), nbags)])

model_stack <- NULL
i=0
for(bb in mask[,1]){
  i=i+1
  foo <- readRDS(loops[bb])
  model_stack[[i]] <-foo
}
saveRDS(list(mod.aic, mod.r2), paste("./outdir/eval_silva_", trait, ".rds", sep=""))
#save the model ensemble in an R object
saveRDS(model_stack, file = paste("./outdir/EPBMs/100_s_", trait, ".rds", sep=""))

# }else{
model_stack <- readRDS(paste("./outdir/EPBMs/100_s_", trait, ".rds", sep=""))
# }

library(tidyverse)
# calculate, scale the dAIC to rank and weight each model using a softmax function
mod.aic=rep(0,length(model_stack))
mod.r2=rep(0,length(model_stack))
for(bb in 1: length(model_stack)){
  mod.aic[bb] <- model_stack[[bb]]$mod$FinalModel$aic
  mod.r2[bb] <- - model_stack[[bb]]$pR2

}
#selected_mods = which(mod.aic %in% sort(mod.aic, decreasing = F)[1:nbags])
#mod.aic = mod.r2
mod.aic <- scale(mod.aic)
delta.aic <- mod.aic - min(mod.aic)
weights <- softmax(-0.5*delta.aic)

#get response for training crowns
train.data.y <- read.csv("./indir/Traits/Chapter1_field_data.csv") %>%
  dplyr::select(c(individualID, trait))
#just_for_
train.data.y = train.data.y[complete.cases(train.data.y),]
#get OOB data crowns
test.data.y = readr::read_csv("./indir/Misc/oob_ids.csv")
#cleaning out of bag test Y and X
# test.data.y <- read.csv("./indir/Traits/Chapter1_field_data.csv") %>%
#   filter(individualID %in% oob$individualID) %>%
#   dplyr::select(c("individualID", trait, "SITE")) %>%
#   group_by(individualID)%>%
#   summarize_if(is.numeric, mean)
test.data.x <- read.csv("./indir/Spectra/reflec2ce_all.csv") %>%
  filter(individualID %in% test.data.y$individualID)

aug.mat = inner_join(test.data.x, test.data.y)
# aug.mat = aug.mat %>%
#   dplyr::group_by(individualID, band_site) %>%
#   #slice_min(n = 20)
#   top_n(20, wt = "band_124") #
tmp_features<- aug.mat[grepl("band", names(aug.mat))]
tmp_variables <- aug.mat[names(aug.mat) %in% trait]
tmp_variables <- (round(tmp_variables,2))
bnd_site <- aug.mat[["SITE"]] %>% factor %>% fastDummies::dummy_cols()
colnames(bnd_site) <- stringr::str_replace(colnames(bnd_site), ".data_", "band_")
#aug.spectra <- data.frame(aug.spectra["individualID"], aug.spectra["taxonID"], bnd_site[-1], tmp_features, tmp_variables)

test.data.x <- data.frame(bnd_site[-1], tmp_features, tmp_variables)
test.data.x <- test.data.x %>% dplyr::select(-one_of(trait))

crownID = aug.mat["individualID"]
test.data.x <- dplyr::select(test.data.x, colnames(model_stack[[1]]$mod$dataX))
#transform features to mirror train features structure
#test.data.x[test.data.x==0] = 0.0000001
#test.data.x=test.data.x[, colSums(is.na(test.data.x)) == 0]
# if(normz==T){
#   foot <-t(diff(t(log(test.data.x[,-c(1:nsites)])),differences=1, lag=3))
#   test.data.x <- cbind(test.data.x[,c(1:nsites)], foot)
# }

test.PLS = as.matrix(test.data.x)

#initialize variabiles
rm(output)
out <- list()
pred.val.data <- list()
output.daic =  output.up.daic = output.lw.daic =  crownID
output.daic$yhat = output.up.daic$yhat = output.lw.daic$yhat =0

for(bb in 1:length(weights)){
  md = bb
  pls.mod.train <- model_stack[[md]]$mod
  optim.ncomps <- model_stack[[md]]$ncomp
  #make predictions using the ith model
  ith_mod_prediction <- pls_glm_predict(pls.mod.train, newdata = test.PLS,
                                        ncomp=optim.ncomps,  type='response')
  ith_mod_prediction=exp(ith_mod_prediction)

  out$upper <-   ith_mod_prediction[,3]
  out$lower <-  ith_mod_prediction[,2]
  output.daic$yhat <- output.daic$yhat + ith_mod_prediction[,1] * weights[bb]
  output.up.daic$yhat <- output.up.daic$yhat + ith_mod_prediction[,3] * weights[bb]
  output.lw.daic$yhat <- output.lw.daic$yhat + ith_mod_prediction[,2] * weights[bb]
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind.data.frame(crownID, rep(bb, dim(test.data.x)[1]), as.vector(ith_mod_prediction[,1]))
  }else{
    output <- rbind.data.frame(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]),
                                                       as.vector(ith_mod_prediction[,1]))))
  }
}
colnames(output) <- c("individualID", "modelID", "yhat")
colnames(output.daic) <- c("individualID", "yhat")
#compare pixel predicions with crowns
output <- inner_join(output, test.data.y[c("individualID", trait)], by = "individualID")
output.daic <- inner_join(output.daic, test.data.y[c("individualID", trait)], by = "individualID")
output = output[complete.cases(output),]
output.daic = output.daic[complete.cases(output.daic),]

pbm_all <- 1 - sum((as.numeric(output[["yhat"]]) - (output[[trait]]))^2) / sum((output[[trait]] - mean(output[[trait]]))^2)
epbm_r2 <-  1 - sum((output.daic$yhat - (output.daic[[trait]]))^2) /
  sum((output.daic[[trait]] - mean(output.daic[[trait]]))^2)

#crown aggregation R2
crown.based.daic <- output.daic %>%
  dplyr::select(c("yhat", trait, "individualID")) %>%
  group_by(individualID) %>%
  summarise_if(is.numeric, mean)
ceam_r2 <-  1 - sum((crown.based.daic$yhat - (crown.based.daic[[trait]]))^2) /
  sum((crown.based.daic[[trait]] - mean(crown.based.daic[[trait]]))^2)

# #calculation of RMSE
rmse_crown <- rmse(crown.based.daic[[trait]], crown.based.daic$yhat)
rmse_pix_ensamble <- rmse(output.daic[[trait]], output.daic$yhat)
rmse_pix <- rmse( output[[trait]], output$yhat)

nrmse_crown = rmse_crown / abs(diff(range(crown.based.daic[[trait]])))
nrmse_ens = rmse_pix_ensamble / abs(diff(range(output.daic[[trait]])))
nrmse_pix = rmse_pix / abs(diff(range(output[[trait]])))
#calculate probability intervals
pix.up <- inner_join(output.up.daic, test.data.y[c("individualID", trait)], by = "individualID")
colnames(pix.up) <- c("individualID", "y_up", "y")
cr.up <- pix.up %>%
  group_by(individualID) %>%
  summarise(y95pi = mean(y_up))
cr.up <- inner_join(cr.up,  test.data.y[c("individualID", trait)], by = "individualID")
colnames(cr.up) <- c("individualID", "y_up", "y")


pix.lw <- inner_join(output.lw.daic,  test.data.y[c("individualID", trait)], by = "individualID")
colnames(pix.lw) <-  c("individualID", "y_lw", "y")
cr.lw <- pix.lw %>%
  group_by(individualID) %>%
  summarise(y5pi = mean(y_lw))
cr.lw <- inner_join(cr.lw,  test.data.y[c("individualID", trait)], by = "individualID")
colnames(cr.lw) <- c("individualID", "y_lw", "y")

test_results = list(y = list(pbm = output[,4], epbm = output.daic[,3], ceam = crown.based.daic[,3]),
                    y_hat = list(pbm = output$yhat, epbm = output.daic$yhat, ceam = crown.based.daic$yhat),
                    y_95 = list(epbm = pix.up, ceam = cr.up),
                    y_5 = list(epbm = pix.lw, ceam = cr.lw),
                    r2 = list(pbm = pbm_all, epbm = epbm_r2, ceam = ceam_r2)
                    ,rmse = list(pbm = rmse_pix, epbm = rmse_pix_ensamble, ceam = rmse_crown)
)
saveRDS(test_results, paste("./outdir/performance_f100_", trait, ".rds", sep=""))
#  return(test_results)
#}




rm(output)
for(bb in 1:length(weights)){
  md = bb
  pls.mod.train <- model_stack[[md]]$mod
  optim.ncomps <- model_stack[[md]]$ncomp
  #make predictions using the ith model
  ith_mod_prediction <- pls_glm_predict(pls.mod.train, newdata = test.PLS,
                                        ncomp=optim.ncomps,  type='response')
  ith_mod_prediction=exp(ith_mod_prediction)
  #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
  if(!exists("output")){
    output <- cbind.data.frame(crownID, rep(bb, dim(test.data.x)[1]), ith_mod_prediction, aug.mat[trait])
    colnames(output) = c("individualID", "mod", "yhat", "y5", "y95", "y")
  }else{
    foo <- cbind.data.frame(crownID, rep(bb, dim(test.data.x)[1]), ith_mod_prediction, aug.mat[trait])
    colnames(foo) = c("individualID", "mod", "yhat", "y5", "y95", "y")
    output <- rbind.data.frame(output, foo)
  }
}
lw = output$y > output$y5
up = output$y < output$y95
sum(lw*up)/nrow(output)

lw = pix.lw$y > pix.lw$y_lw
up = pix.up$y < pix.up$y_up
sum(lw*up)/nrow(pix.lw)

lw = cr.lw$y > cr.lw$y_lw
up = cr.up$y < cr.up$y_up
sum(lw*up)/nrow(cr.up)

ceam_r2
epbm_r2
pbm_all

rmse_crown
rmse_pix_ensamble
rmse_pix

nrmse_crown
nrmse_ens
nrmse_pix
