#' split data into train-validation set
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#'
#'
cut.set<-function(aug.X,out.dir, c.id, prop = 0.7){
  species <- unique(aug.X$aug.spectra.name)
  train.data <- 0
  test.data <- 0
  j <- 1
  cr.id <- NULL
  #loop through species to create a stratified dataset split
  for (i in as.character(species)){
    set.seed(1)
    #subset by species
    temp.data <- aug.X[which(aug.X$aug.spectra.name==i),]
    rows <- sample(1:nrow(temp.data),ceiling(prop*nrow(temp.data)))
    foo <- c.id[which(aug.X$aug.spectra.name==i)]
    cr.id <- c(cr.id, foo[rows])
    cal.data = droplevels(temp.data[rows,])
    val.data = droplevels(temp.data[-rows,])

    if(j==1){
      train.data <- cal.data
      test.data <- val.data
    } else {
      train.data <- rbind(train.data,cal.data)
      test.data <- rbind(test.data,val.data)
    }

    j <- j+1
  }
  return(list(train=train.data, test=test.data, cr.id=cr.id))
}
