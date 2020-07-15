#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples get_random_bags(1000)
get_random_bags<- function(allData, lp = 1, seed = 1987){

  unqCrown = unique(allData["individualID"])
  #get number of individual crowns
  #bootDat <- allData[1:nrow(unqCrown), ]
  tk = 1
  bootDat = list()
  for(i in unlist(unqCrown)){
    tmp_id = allData %>% filter(individualID ==i)
    if(nrow(tmp_id)==1){
      #only one pixel in the bag. Use that one for the permutation
      bootDat[[tk]] <- allData[which(allData["individualID"]==i),]
    }else{
      # set reproducible random seed different for each of the n extractions
      set.seed(seed)
      rraand <- lp *runif(1, 1, 10^3)
      set.seed(as.integer(rraand)) # todo: change laps * odd number?
      #extract one pixel per bag and create the ith dataset
      bootDat[[tk]] <- tmp_id %>% sample_n(1)
    }
    tk=tk+1
  }
  bootDat = do.call(rbind.data.frame, bootDat)
  write_csv(bootDat, paste('./indir/Permutations/onePix1Crown_', lp, '.csv', sep = ''))
  return(bootDat)
}
