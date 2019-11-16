#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
build_model <- function(loop=1,  nrmlz = F, trait = "N_pct"){
  #are random combinations of pixels already set up?
  combinations <- file.exists(paste('./indir/Permutations/onePix1Crown_', loop, ".csv", sep=""))
  if(!combinations){
    #extract n combinations of pixles by extracting one per bag
    get_random_bags(lp = loop)
  }
  # run the pls glm on training bags for each random extractions
  for(tr in trait){
    print(tr)
    random_bag_pls <- pls_glm(trait = tr,
                              ll = loop, nrmlz = nrmlz)
  }
}
