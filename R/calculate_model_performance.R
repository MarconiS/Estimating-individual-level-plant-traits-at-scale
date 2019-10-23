#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
calculate_model_performance <- function(trait = "P_pct", snap_path = "", get_best = F){
  #retrieve number of snaps
  loops <- list.files("./outdir/snaps/") %>% length

  #get delta_aic and rank best models
  mod.aic <- rep(NA, loops)
  for(bb in 1: length(mod.aic)){
    foo <- readRDS()
    mod.aic[bb] <- foo[[bb]]$score$aic
   }
  mask <- cbind(which(mod.aic %in% head(sort(mod.aic), 100)), mod.aic[mod.aic %in% head(sort(mod.aic), 100)])

}
