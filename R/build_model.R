#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#'
#'

build_model <- function(loops=1000, cores = 1, nrmlz = T, trait = "P_pct"){
  #are random combinations of pixels already set up?
  combinations <- list.files('./indir/Permutations', pattern=".csv")
  if(length(combinations)!=loops){
    #extract n combinations of pixles by extracting one per bag
    get_random_bags(loops= loops)
  }

  # run the pls glm on training bags for each random extractions
  random_bag_pls <- pls_glm(trait = trait,
                 loops = loops,
                 cores = cores,
                 nrmlz = nrmlz)

  #save trained model object
  save(random_bag_pls, file = paste(out_dir,  "/pls_glm_nodiff", nm,  sep = ""))
}
