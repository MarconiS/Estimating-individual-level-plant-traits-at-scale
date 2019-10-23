#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
calculate_model_performance <- function(trait = "P_pct", snap_path = "", nbags = 1){
  #retrieve number of snaps
  loops <- list.files("./outdir/PBMs", pattern = ".rds", full.names = T)

  #check if Ensemble model has been created yet
  model_stack <- file.exists(paste("./outdir/EPBMs/", trait, ".rds", sep=""))

  if(!model_stack || length(model_stack) != nbags){
    #get delta_aic and rank best models
    mod.aic <- rep(NA, length(loops))
    model_stack <- NULL
    for(bb in 1: length(mod.aic)){
      foo <- readRDS(loops[bb])
      mod.aic[bb] <- foo$mod$AIC[foo$ncomp]
      #create model stack with all random bags models
      model_stack[[bb]] <-foo
    }
    #subset the best n models to be used for performance
    mask <- cbind(which(mod.aic %in% head(sort(mod.aic), nbags)),
                  mod.aic[mod.aic %in% head(sort(mod.aic), nbags)])
    #save the model ensemble in an R object
    saveRDS(model_stack[mask[,1]],
            file = paste("./outdir/EPBMs/", trait, ".rds", sep=""))
  }else{
    model_stack <- readRDS(paste("./outdir/EPBMs/", trait, ".rds", sep="")))
  }

}
