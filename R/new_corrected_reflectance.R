#loop through polygons
library(tidyverse)
library(raster)
data = sf::read_sf("../Chapter1/corrected_res/TALL_crown_polygons/")

#extract_pixels = function(ii, data){
new_corrected_reflectance = list()
for(ii in 1:nrow(data)){
  itc = data[ii,]
  coords = (sf::st_coordinates(itc)[,1:2]/1000)
  coords =as.integer(unlist(coords)) %>% unique
  coords = paste(coords[1] *1000, coords[2]*1000,  sep="_")

  list_files = list.files("//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D08/TALL/L4/corrHSI/",
                          pattern = coords,full.names = T)
  itc_refl = list()
  for(pt in list_files){
    #make sure CRS is set
    f = raster::stack(pt)
    proj4string(f) <- CRS(raster::crs(itc))
    itc_dat = raster::extract(f, itc, df=T)
    if(nrow(itc_dat) !=0){
      colnames(itc_dat) = paste("band", 1:ncol(itc_dat), sep="_")
      rownames(itc_dat) = 1:nrow(itc_dat)
      itc_dat = cbind.data.frame(itc[["tree"]], itc_dat)
      itc_refl[[pt]] = itc_dat
    }
  }
  foo = do.call(rbind.data.frame, itc_refl)
  foo = foo[complete.cases(foo),]
  #return(foo)
  new_corrected_reflectance[[ii]] = foo
#}
}


cl <- makeCluster(mc <- getOption("cl.cores", 22))
clusterExport(cl, c('data'))
results <- parLapply(cl, 1:nrow(data), fun=extract_pixels, data = data )




itc = data[1,]
coords = (sf::st_coordinates(itc)[,1:2]/1000)
coords =as.integer(unlist(coords)) %>% unique
coords = paste(coords[1] *1000, coords[2] *1000,itc$path, sep="_")

list_files = list.files("//orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/corrHSI/",
                        pattern = coords,full.names = T)

f = raster::stack(list_files)
itc_dat = raster::extract(f, itc, df=T)
