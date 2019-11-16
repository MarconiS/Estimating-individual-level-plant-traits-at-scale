#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#' @import Metrics
pls_glm_performance <- function(trait = "N_pct", nbags = 20, normz=F){
  #retrieve number of snaps
  loops <- list.files("./outdir/PBMs", pattern = trait, full.names = T)

  #check if Ensemble model has been created yet
  model_stack <- file.exists(paste("./outdir/EPBMs/", trait, ".rds", sep=""))

  if(!model_stack || length(model_stack) != nbags){
    #get delta_aic and rank best models
    mod.aic <- rep(NA, length(loops))
    model_stack <- NULL
    for(bb in 1: length(mod.aic)){
      foo <- readRDS(loops[bb])
      mod.aic[bb] <- foo$mod$FinalModel$aic
      #create model stack with all random bags models
      #model_stack[[bb]] <-foo
    }
    #subset the best n models to be used for performance
    mask <- cbind(which(mod.aic %in% head(sort(mod.aic), nbags)),
                  mod.aic[mod.aic %in% head(sort(mod.aic), nbags)])

    model_stack <- NULL
    i=0
    for(bb in mask[,1]){
      i=i+1
      foo <- readRDS(loops[bb])
      model_stack[[i]] <-foo
    }
    #save the model ensemble in an R object
    #model_stack <- model_stack[mask[,1]]
    saveRDS(model_stack,
            file = paste("./outdir/EPBMs/", trait, ".rds", sep=""))
  }else{
    model_stack <- readRDS(paste("./outdir/EPBMs/", trait, ".rds", sep=""))
  }
  # calculate, scale the dAIC to rank and weight each model using a softmax function
  mod.aic=rep(0,length(model_stack))
  for(bb in 1: length(model_stack)){
    mod.aic[bb] <- model_stack[[bb]]$mod$FinalModel$aic
  }
  mod.aic <- scale(mod.aic)
  delta.aic <- mod.aic - min(mod.aic)
  weights <- tsensembler::softmax(-0.5*delta.aic)

  #get response for training crowns
  train.data.y <- read.csv("./indir/Traits/CrownTraits.csv") %>%
    dplyr::select(c("individualID", trait))

  #cleaning out of bag test Y and X
  test.data.y <- read.csv("./indir/Tests/CrownTraits_outBag.csv") %>%
    dplyr::select(c("individualID", trait, "siteID"))

  nsites <- length(unique(test.data.y["siteID"]))
  test.data.x <- read.csv("./indir/Tests/CrownPix_outBag_new.csv")
  crownID = test.data.x["individualID"]
  test.data.x <- dplyr::select(test.data.x, colnames(foo$mod$dataX))

  #transform features to mirror train features structure
  test.data.x=test.data.x[, colSums(is.na(test.data.x)) == 0]
  if(normz==T){
    foot <-t(diff(t(log(test.data.x[,-c(1:nsites)])),differences=1, lag=3))
    test.data.x <- cbind(test.data.x[,c(1:nsites)], foot)
  }

  test.PLS = as.matrix(test.data.x)

  #initialize variabiles
  rm(output)
  out <- list()
  pred.val.data <- list()
  output.daic =  output.up.daic = output.lw.daic =  crownID
  output.daic$yhat = output.up.daic$yhat = output.lw.daic$yhat =0

  for(bb in 1:length(weights)){
    pls.mod.train <- model_stack[[bb]]$mod
    optim.ncomps <- model_stack[[bb]]$ncomp
    #make predictions using the ith model
    ith_mod_prediction <- pls_glm_predict(pls.mod.train, newdata = test.PLS,
                                          wt = rep(1, nrow(test.PLS)),
                                          ncomp=optim.ncomps,  type='response')
    ith_mod_prediction=ith_mod_prediction[]
    pred.val.data$fit <- ith_mod_prediction[,1]
    pred.val.data$upper <- ith_mod_prediction[,3]
    pred.val.data$lower <- ith_mod_prediction[,2]
    out$output.daic <-  pred.val.data$fit
    out$upper <-  pred.val.data$upper
    out$lower <-  pred.val.data$lower
    output.daic$yhat <- output.daic$yhat + pred.val.data$fit * weights[bb]
    output.up.daic$yhat <- output.up.daic$yhat + pred.val.data$upper * weights[bb]
    output.lw.daic$yhat <- output.lw.daic$yhat + pred.val.data$lower * weights[bb]
    #you have then a vector of predicions whose legnth is sum(crID_i * pixels_i)
    if(!exists("output")){
      output <- cbind.data.frame(crownID, rep(bb, dim(test.data.x)[1]), as.vector(pred.val.data$fit))
    }else{
      output <- rbind.data.frame(output, as.matrix(cbind(crownID, rep(bb, dim(test.data.x)[1]),
                                              as.vector(pred.val.data$fit))))
    }
  }
  colnames(output) <- c("individualID", "modelID", "yhat")
  colnames(output.daic) <- c("individualID", "yhat")

  #compare pixel predicions with crowns
  output <- inner_join(output, test.data.y, by = "individualID")
  output.daic <- inner_join(output.daic, test.data.y, by = "individualID")
  pbm_all <- 1 - sum((as.numeric(output[["yhat"]]) - (output[[trait]]))^2) / sum((output[[trait]] - mean(output[[trait]]))^2)
  epbm_r2 <-  1 - sum((output.daic$yhat - (output.daic[[trait]]))^2) /
    sum((output.daic[[trait]] - mean(output.daic[[trait]]))^2)

  #crown aggregation R2
  crown.based.daic <- output.daic %>%
    group_by(individualID) %>%
    summarise(yhat = median(yhat))
  crown.based.daic <- inner_join(crown.based.daic, test.data.y, by = "individualID")
  ceam_r2 <-  1 - sum((crown.based.daic$yhat - (crown.based.daic[[trait]]))^2) /
    sum((crown.based.daic[[trait]] - mean(crown.based.daic[[trait]]))^2)

  # #calculation of RMSE
  # rmse_crown <- rmse(crown.based.daic[[trait]], crown.based.daic$yhat)
  # rmse_pix_ensamble <- rmse(output.daic[[trait]], output.daic$yhat)
  # rmse_pix <- rmse( output[[trait]], output$yhat)

  #calculate probability intervals
  pix.up <- inner_join(output.up.daic, test.data.y, by = "individualID")
  colnames(pix.up)[2] <- "y_up"
  cr.up <- pix.up %>%
    group_by(individualID) %>%
    summarise(y95pi = median(y_up))
  cr.up <- inner_join(cr.up, test.data.y, by = "individualID")
  colnames(cr.up) <- c("individualID", "y_up", "y")


  pix.lw <- inner_join(output.lw.daic, test.data.y, by = "individualID")
  colnames(pix.lw)[2] <- "y_lw"
  cr.lw <- pix.lw %>%
    group_by(individualID) %>%
    summarise(y5pi = median(y_lw))
  cr.lw <- inner_join(cr.lw, test.data.y, by = "individualID")
  colnames(cr.lw) <- c("individualID", "y_lw", "y")

  test_results = list(y_hat = list(pbm = output$yhat, epbm = output.daic$yhat, ceam = crown.based.daic$yhat),
                      y_95 = list(epbm = pix.up, ceam = cr.up),
                      y_5 = list(epbm = pix.lw, ceam = cr.lw),
                      r2 = list(pbm = pbm_all, epbm = epbm_r2, ceam = ceam_r2),
                     #rmse = list(pbm = rmse_pix, epbm = rmse_pix_ensamble, ceam = rmse_crown)
                      )
  return(test_results)
}
