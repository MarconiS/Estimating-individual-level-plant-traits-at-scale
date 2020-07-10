clean_spectra <- function(brick, ndvi = 0.3, nir = 0.1, outlier = F){
  # filter for no data
  brick = brick %>% ungroup %>% dplyr::select(contains("band"))
  brick = brick %>% dplyr::select(-one_of("band_site"))
  mask1 = apply(brick[,1:369], 1, function(x)all(x>0.001))
  mask2 = apply(brick[,1:369], 1, function(x)all(x<1))
  brick[!as.logical(mask1 *mask2), ] = NA

  #filter for greennes and shadows
  ndvi <- (brick[,"band_90"]- brick[,"band_58"])/(brick[,"band_58"] + brick[,"band_90"]) <ndvi
  nir860 <- (brick[,"band_96"] + brick[,"band_97"])/2 < nir
  mask = as.logical(ndvi * nir860)
  mask[is.na(mask)] = T
  brick[mask,] = NA
  rm(mask, ndvi, nir860)
  brick = brick %>% dplyr::select(one_of(paste("band", 10:360, sep="_")))
  normMat <- sqrt(apply(brick^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat <- matrix(data=rep(normMat,ncol(brick)),ncol=ncol(brick))
  brick=brick/normMat
  rm(normMat)

  #filter for known artifacts
  cnd = (brick[,298] > 0.03)
    idx <- (apply(data.frame(cnd), 1, any))
  if(length(idx) !=0){
    idx[is.na(idx)] = T
    brick[idx,] = NA
  }

  cnd = (brick[,24:45] > 0.03)
  idx <- (apply(cnd, 1, any))
  if(length(idx) !=0){
    idx[is.na(idx)] = T
    brick[idx,] = NA
  }
  cnd = (brick[,195:200] > 0.043)
  idx <- (apply(cnd, 1, any))
  if(length(idx) !=0){
    idx[is.na(idx)] = T
    brick[idx,] = NA
  }
  rm(cnd,idx)

  # # save pixel positions
  good_pix = !is.na(brick)
  good_pix = (apply(good_pix, 1, all))
  #
  brick = brick[complete.cases(brick),]
  return(list(refl = brick, good_pix = good_pix))
}
