#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples get_random_bags(1000)
get_random_bags<- function(allData, lp = 1, seed = 1987){

  unqCrown = unique(allData["individualID"])
  #get number of individual crowns
  bootDat <- allData[1:nrow(unqCrown), ]
  bootDat[] <- NA
  tk=1
  for(i in unlist(unqCrown)){
    if(sum(allData["individualID"]==i)==1){
      #only one pixel in the bag. Use that one for the permutation
      bootDat[tk,] <- allData[which(allData["individualID"]==i),]
    }else{
      # set reproducible random seed different for each of the n extractions
      rraand <- lp #*runif(1, 1, 10^6)
      set.seed(rraand) # todo: change laps * odd number?
      #extract one pixel per bag and create the ith dataset
      bootDat[tk,] <- allData[sample(which(allData["individualID"]==i), 1),]
    }
    tk=tk+1
    write_csv(bootDat, paste('./indir/Permutations/onePix1Crown_', lp, '.csv', sep = ''))
  }
  return(bootDat)
}
