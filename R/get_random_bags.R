#' retrieve vegetation structure data  from NEON
#'
#'
#' @return
#' @export
#' @examples get_random_bags(1000)
#'
#'
get_random_bags<- function(loops = 1, seed = 1987){
  #get dataset paths
  data_path = list.files("./indir/Spectra/", pattern = ".csv", full.names = T)
  # Read data
  allData=read_csv(data_path)

  unqCrown = unique(allData["pixel_crownID"])
  #get number of individual crowns
  bootDat <- allData[1:nrow(unqCrown), ]
  bootDat[] <- NA
  for (laps in 1:loops){
    tk = 1
    for(i in unlist(unqCrown)){
      if(sum(allData["pixel_crownID"]==i)==1){
        #only one pixel in the bag. Use that one for the permutation
        bootDat[tk,] <- allData[which(allData["pixel_crownID"]==i),]
      }else{
        rraand <- laps *runif(1, 1, 10^6)
        set.seed(rraand) # todo: change laps * odd number?
        #extract one pixel per bag and create the ith dataset
        bootDat[tk,] <- allData[sample(which(allData["pixel_crownID"]==i), 1),]
      }
      tk = tk +1
    }
    print(laps)
    write_csv(bootDat, paste('./indir/Permutations/onePix1Crown_', laps, '.csv', sep = ''))
  }
  #return(bootDat)
}
