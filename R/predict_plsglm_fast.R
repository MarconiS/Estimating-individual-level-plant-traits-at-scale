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


