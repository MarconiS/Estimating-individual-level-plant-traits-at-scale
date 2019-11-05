#' clean the permutation folder from all bag extractions
#'
#'
#' @return
#' @export
#' @examples
clean_permutations <- function(){
  do.call(file.remove, list(list.files("./indir/Permutations/", full.names = TRUE)))
}
