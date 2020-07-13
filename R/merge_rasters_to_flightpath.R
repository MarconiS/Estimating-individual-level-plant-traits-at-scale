#merge rasters to flightpaths
fld = "/Users/sergiomarconi/Dropbox (UFL)/Traits_Maps/OSBS/Cpercent"
ls.fl = list.files(fld, pattern = ".tif")
flpths = do.call(rbind.data.frame, strsplit(ls.fl, split = "_"))
colnames(flpths) = c("X", "Y", "pt")
pt = unique(flpths$pt)

library(gdalUtils)
library(rgdal)
library(raster)
library(tidyverse)
for (i in pt){
  lf = flpths %>% dplyr::filter(pt == i)
  ras = list()
  for(jj in 1:nrow(lf)){
    ras[[jj]] = raster::brick(paste(fld, paste(lf[jj,], collapse = "_"), sep="/"))
  }
  m <- do.call(raster::merge, ras)
  m[m$layer.1 >56,]=NA
  writeRaster(m, paste("/Users/sergiomarconi/Dropbox (UFL)/Traits_Maps/OSBS/Cpercent_", i, sep=""), overwrite=TRUE)
}
