#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#'
#'maybe you want ot batchtool from here: get teh seed, get the split, build the model
build_model <- function(loop=1,  nrmlz = T, trait = "P_pct"){
  #are random combinations of pixels already set up?
  combinations <- file.exists(paste('./indir/Permutations/onePix1Crown_', loop, ".csv", sep=""))
  if(!combinations){
    #extract n combinations of pixles by extracting one per bag
    get_random_bags(lp= loop)
  }
  # run the pls glm on training bags for each random extractions
  random_bag_pls <- pls_glm(trait = trait,
                 ll = loop, nrmlz = nrmlz)
}
