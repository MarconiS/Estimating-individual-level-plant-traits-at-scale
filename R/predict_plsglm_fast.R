plsglm_fast_predict <- function(object,newdata,
                            comps=object$computed_nt, type=c("link", "response", "terms", "scores", "class", "probs"),
                            se.fit=FALSE, wt = NULL, dispersion = NULL,methodNA="adaptative",verbose=TRUE,...)
{
  type <- match.arg(type)
  nrnd <- nrow(newdata)
  newdataNA <- !is.na(newdata)
  newdata <- sweep(sweep(newdata, 2, attr(object$ExpliX, "scaled:center")),
                   2, attr(object$ExpliX, "scaled:scale"), "/")
  newdata <- as.matrix(newdata)
  #newdata[!newdataNA] <- 0
  ttpredY = lapply(1:nrnd, function(x)c(newdata[x,] %*% object$wwetoile[, 1:comps],
                                        rep(0, object$computed_nt - comps)))
  ttpredY = do.call(rbind, ttpredY)
  colnames(ttpredY) <- NULL
  ttpred <- data.frame(tt = ttpredY)
  pred_int = HH::interval(object$FinalModel, newdata=ttpred, type="response")
  pred_int = pred_int[,c(1,4,5)]
  colnames(pred_int) = c("fit","lwr","upr")
  return(pred_int[])
}


