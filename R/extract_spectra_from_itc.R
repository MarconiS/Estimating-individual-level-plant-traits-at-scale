itcs = readr::read_csv("../neonVegWrangleR/outdir/field_data/vst_latest_2m.csv")
itcs = readr::read_csv("../Chapter2/indir/final_field_dataset.csv")
itcs = itcs %>% dplyr::filter(canopyPosition %in% c("Full sun", "Partially shaded", "Open grown"))
itcs = sf::st_as_sf(itcs, coords = c("longitude", "latitude"), crs = 4326)
plot(itcs[1])
raster_path = "/Users/smarconi/Documents/Data/bdrf"
extract_spectra_from_itcs <- function(itcs, raster_path){
  data = NULL
  for(plt in unique(itcs$ID_x)){
    st_ic = itcs %>% dplyr::filter(ID_x == plt)
    tryCatch({
      tiles = list.files(raster_path, pattern = as.character(plt), full.names = T)
      for(pt in tiles){
        hsi = raster::brick(pt)
        raster::crs(hsi) = raster::crs(st_ic)
        #extract hsi with a buffer if point, or exact polygon if itc
        #if(class(st_geometry(st_ic))[1] == "sfc_POINT"){
        foo = raster::extract(hsi, st_ic)#, buffer = 1)
        # add individual id to extracte pixels
        individualID = st_ic$ID_x
        foo = lapply(1:length(individualID), function(x) cbind(individualID[[x]], foo[[x]]))
        foo = do.call(rbind.data.frame, foo)

        colnames(foo) = c("individualID", paste("band", 1:369, sep = "_"))
        data[[plt]] = foo
      }
    },error=function(e){message(paste(plt,"tile is missing"))})
  }
  data = do.call(rbind.data.frame, data)
  readr::write_csv(data, "./indir/spectra/silva2px_traits_reflectance.csv")
}

extract_spectra_from_itcs(itcs, raster_path)
