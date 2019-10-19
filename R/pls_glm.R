#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#'
#'
pls_glm <- function(trait = "N_pct",cores = 1, loops = 1, nrmlz=T){

  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())

  ll = as.character(1:loops)
  mod.out <- foreach(laps = ll, .verbose = T, .multicombine = TRUE) %dopar% {

    #get traits data
    cr.Traits <- read_csv(paste("./indir/Traits/CrownTraits.csv",sep=""))

    #get random flip spectra data
    aug.spectra <- read_csv(paste('./indir/Permutations/onePix1Crown_', laps, '.csv', sep = ''))
    aug.spectra <- merge(cr.Traits, aug.spectra, by = "pixel_crownID")
    tmp_features<- aug.spectra[grepl("band", names(aug.spectra))]
    tmp_features <- tmp_features[-(which(colnames(tmp_features) %in% c("band_1", "band_370")))]
    tmp_variables <- aug.spectra[names(aug.spectra) %in% trait]
    tmp_variables <- (round(tmp_variables,3))
    aug.X <- data.frame(aug.spectra$name, tmp_features, tmp_variables)

    # Subset data into cal/val by site
    eval.set <- cut.set(aug.X,out.dir, aug.spectra$pixel_crownID)
    train.data <- eval.set$train
    test.data <- eval.set$test
    colnames(train.data)[1] <- "name"
    colnames(test.data)[1] <- "name"
    #predictors matrix
    X <- as.matrix(train.data[grepl("band", names(train.data))])
    X[X==0] = 0.0000001
    #responses matrix
    Y <- as.vector(train.data[,names(train.data) %in% trait])
    #K=nrow(Y)
    #get sites and set up sites fixed effects
    nsites <- length(unique(aug.spectra$site))
    if(nrmlz==T){
      X.n <-t(diff(t(log(X[,-c(1:nsites)])),differences=1, lag=3))
      X <- cbind(X[,c(1:nsites)], X.n)
    }

    #perform a cross-valiadation on train-validation set
    train.PLS<- cv.plsRglm(dataY=(Y),dataX=X,
                           nt=15,NK=1, K=10,
                           modele="pls-glm-family",family=gaussian(), verbose = F)
    out <- list()
    #define number of components to chose by calculating PRESS statistics
    press <- unlist(kfolds2Press(train.PLS))
    out["ncomp"] = which(press == min(press))

    # retrain on chosen number of components
    mod <- plsRglm(dataY=(Y),dataX=X,as.integer(out["ncomp"]),modele="pls-glm-gaussian", verbose = F)
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    X.tst[X.tst==0] = 0.0000001
    if(nrmlz==T){
      X.ntst <-t(diff(t(log(X.tst[,-c(1:nsites)])),differences=1, lag=3))
      X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    }

    #store validation metrics
    Y.test <- as.vector(test.data[,names(test.data) %in% trait])
    out[["pred"]] <- predict(mod,newdata=X.tst,type="response",comps=as.integer(out["ncomp"]))
    out[["mod"]] <- mod
    #out[["score"]] <- gtrait(out$pred ~ (Y.test), family = gaussian())
    #out["aic"]<- AIC(out$score)
    out
  }
  stopCluster(cl)
  return(mod.out)
}
