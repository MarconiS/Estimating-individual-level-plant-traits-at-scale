#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#' @import plsRglm
pls_glm <- function(ll = NULL, trait = NULL, nrmlz=F){
    #get traits data
    cr.Traits <- readr::read_csv(paste("./indir/Traits/CrownTraits.csv",sep=""))

    #get random flip spectra data
    aug.spectra <- readr::read_csv(paste('./indir/Permutations/onePix1Crown_', ll, '.csv', sep = ''))
    aug.spectra <- merge(cr.Traits, aug.spectra, by = "individualID")
    tmp_features<- aug.spectra[grepl("band", names(aug.spectra))]
    tmp_variables <- aug.spectra[names(aug.spectra) %in% trait]
    tmp_variables <- (round(tmp_variables,3))
    bnd_site <- aug.spectra[["siteID"]] %>% factor %>% fastDummies::dummy_cols()
    colnames(bnd_site) <- stringr::str_replace(colnames(bnd_site), ".data_", "band_")
    aug.X <- data.frame(aug.spectra["taxonID"], bnd_site[-1], tmp_features, tmp_variables)

    # Subset data into cal/val by site
    eval.set <- cut_set(aug.X, c.id = aug.spectra[["individualID"]])
    train.data <- eval.set$train
    test.data <- eval.set$test
    colnames(train.data)[1] <- "taxonID"
    colnames(test.data)[1] <- "taxonID"
    #predictors matrix
    X <- as.matrix(train.data[grepl("band", names(train.data))])
    X[X==0] = 0.0000001
    #responses matrix
    Y <- as.vector(train.data[,names(train.data) %in% trait])
    #K=nrow(Y)
    #get sites and set up sites fixed effects
    nsites <- nrow(unique(aug.spectra["siteID"]))
    if(nrmlz==T){
      X.n <-t(diff(t(log(X[,-c(1:nsites)])),differences=1, lag=3))
      X <- cbind(X[,c(1:nsites)], X.n)
    }
    if(nsites==1){
      X = X[,-1]
    }
    #perform a cross-valiadation on train-validation set
    train.PLS<- plsRglm::cv.plsRglm(dataY=(Y),dataX=X,
                           nt=15,NK=1, K=10,
                           modele="pls-glm-family",family=gaussian(), verbose = F)
    out <- list()
    #define number of components to chose by calculating PRESS statistics
    press <- unlist(plsRglm::kfolds2Press(train.PLS))
    out["ncomp"] = which(press == min(press))

    # retrain on chosen number of components
    mod <- plsRglm::plsRglm(dataY=(Y),dataX=X,as.integer(out["ncomp"]),
                            modele="pls-glm-gaussian", verbose = F)
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    X.tst[X.tst==0] = 0.0000001
    if(nrmlz==T){
      X.ntst <-t(diff(t(log(X.tst[,-c(1:nsites)])),differences=1, lag=3))
      X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    }
    if(nsites==1){
      X.tst = X.tst[,-1]
    }

    #store validation metrics
    Y.test <- as.vector(test.data[,names(test.data) %in% trait])
    out[["pred"]] <- predict(mod, newdata=X.tst,
                             type="response",comps=as.integer(out["ncomp"]))
    out[["pR2"]] <- 1 - sum((out[["pred"]] - (Y.test))^2) /
      sum((Y.test - mean(Y.test))^2)
    out[["mod"]] <- mod
    saveRDS(out, paste("./outdir/PBMs/pls_glm_", trait, ll, ".rds", sep=""))
  return(out)
}
