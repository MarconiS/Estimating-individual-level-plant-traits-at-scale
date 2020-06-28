#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#' @import plsRglm
pls_glm <- function(ll = NULL, trait = NULL, nrmlz=F){

    wrangle = function(x){
      ifelse(is.character(x), x[1], mean(x))
    }
    #get traits data
    cr.Traits <- readr::read_csv(paste("./indir/Traits/Chapter1_field_data.csv",sep=""))
    cr.Traits = cr.Traits %>% group_by(individualID) %>% summarize_all(wrangle)
    train_ids = readr::read_csv("./indir/Misc/train_ids.csv")
    calib_ids = readr::read_csv("./indir/Misc/calibration_ids.csv")
    test_ids = readr::read_csv("./indir/Misc/oob_ids.csv")
      #cr.Traits %>% filter(individualID %in% spectra$individualID) %>% group_by(taxonID, SITE)%>%
     # sample_frac(0.2)
    #get random flip spectra data
    aug.spectra <- readr::read_csv(paste('./indir/Permutations/onePix1Crown_', ll, '.csv', sep = ''))
    #aug.spectra = spectra_ave
    aug.spectra <- inner_join(cr.Traits, aug.spectra, by = "individualID")
    aug.spectra = aug.spectra %>% filter(!individualID %in% test_ids)
    tmp_features<- aug.spectra[grepl("band", names(aug.spectra))]
    tmp_variables <- aug.spectra[names(aug.spectra) %in% trait]
    tmp_variables <- (round(tmp_variables,2))
    bnd_site <- aug.spectra[["SITE"]] %>% factor %>% fastDummies::dummy_cols()
    colnames(bnd_site) <- stringr::str_replace(colnames(bnd_site), ".data_", "band_")
    aug.spectra <- data.frame(aug.spectra["individualID"], aug.spectra["taxonID"], bnd_site[-1], tmp_features, tmp_variables)
    aug.spectra = aug.spectra[complete.cases(aug.spectra),]
    # Subset data into cal/val by site
    train.data = aug.spectra %>% filter(!individualID %in% calib_ids$calib_id)
    test.data = aug.spectra %>% filter(individualID %in% calib_ids$calib_id)
    #eval.set <- cut_set(aug.X, c.id = aug.spectra[["individualID"]])
    #train.data <- eval.set$train
    #test.data <- eval.set$test
    #colnames(train.data)[1] <- "taxonID"
    #colnames(test.data)[1] <- "taxonID"
    #predictors matrix
    train.data= train.data %>% select(-one_of("band_site"))
    test.data= test.data %>% select(-one_of("band_site"))

    X <- as.matrix(train.data[grepl("band", names(train.data))])
    #responses matrix
    Y <- as.vector(train.data[,names(train.data) %in% trait])
    #K=nrow(Y)
    #get sites and set up sites fixed effects
    nsites <- nrow(unique(aug.spectra["band_site"]))
    # if(nrmlz==T){
    #   X.n <- t(diff(t((X[,-c(1:nsites)])),differences=1, lag=3))
    #   X <- cbind(X[,c(1:nsites)], X.n)
    # }
    if(nsites==1){
      X = X[,-1]
    }
    #perform a cross-valiadation on train-validation set
    train.PLS<- plsRglm::cv.plsRglm(dataY=(Y),dataX=X, scaleY = T, verbose=F,
                                    nt=15,NK=1, K=5,
                                    modele="pls-glm-family",family=gaussian(link = "log"))
    out <- list()
    #define number of components to chose by calculating PRESS statistics
    press <- unlist(plsRglm::kfolds2Press(train.PLS))
    out["ncomp"] = which(press == min(press))

    # retrain on chosen number of components
    mod <- plsRglm::plsRglm(dataY=(Y),dataX=X,as.integer(out["ncomp"]), scaleY = T,
                            modele="pls-glm-family",family=gaussian(link = "log"))
    X.tst <- as.matrix(test.data[grepl("band", names(test.data))])
    # if(nrmlz==T){
    #   X.ntst <-t(diff(t((X.tst[,-c(1:nsites)])),differences=1, lag=3))
    #   X.tst <- cbind(X.tst[,c(1:nsites)], X.ntst)
    # }
    if(nsites==1){
      X.tst = X.tst[,-1]
    }

    #store validation metrics
    Y.test <- as.vector(test.data[,names(test.data) %in% trait])
    out[["pred"]] <- (predict(mod, newdata=X.tst,
                             type="response",comps=as.integer(out["ncomp"])))
    out[["pR2"]] <- 1 - sum((out[["pred"]] - (Y.test))^2) /
      sum((Y.test - mean(Y.test))^2)
    out[["mod"]] <- mod
    saveRDS(out, paste("./outdir/PBMs/pls_glm_", trait, ll, ".rds", sep=""))
  return(out)
}
