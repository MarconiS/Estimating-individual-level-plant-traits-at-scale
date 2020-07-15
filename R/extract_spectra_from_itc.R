# itcs = readr::read_csv("../neonVegWrangleR/outdir/field_data/vst_latest_2m.csv")
# itcs = readr::read_csv("../Chapter2/indir/final_field_dataset.csv")
# itcs = itcs %>% dplyr::filter(canopyPosition %in% c("Full sun", "Partially shaded", "Open grown"))
#itcs = sf::st_as_sf(itcs, coords = c("longitude", "latitude"), crs = 4326)
raster_path = "../../Data/ch1_plots"
#itcs = sf::read_sf("../../Data/Dimensions/OSBS_crown_polygons/OSBS_sample_polygons_edits_Feb2018.shp")
itcs = sf::read_sf("../../Data/Dimensions/silva_tall.shp")
itcs$ID_x = itcs$treeID
itcs = itcs %>%filter(!is.na(ID_x))
extract_spectra_from_itcs <- function(itcs, raster_path){
  data = pix = NULL
  for(plt in unique(itcs$ID_x)){
    st_ic = itcs %>% dplyr::filter(ID_x == plt)
      tiles = list.files(raster_path, pattern = as.character(plt), full.names = T)
      for(pt in tiles){
        tryCatch({
        hsi = raster::brick(pt)
        raster::crs(hsi) = raster::crs(st_ic)
        #extract hsi with a buffer if point, or exact polygon if itc
        #if(class(st_geometry(st_ic))[1] == "sfc_POINT"){
        st_ic = sf::st_centroid(st_ic)
        foo = raster::extract(hsi, st_ic, buffer = 2)
        # add individual id to extracte pixels
        individualID = st_ic$ID_x
        foo = lapply(1:length(individualID), function(x) cbind(individualID[[x]], foo[[x]]))
        foo = do.call(rbind.data.frame, foo)
        colnames(foo) = c("individualID", paste("band", 1:369, sep = "_"))
        data[[pt]] = foo
        pix[[pt]] = nrow(foo)
      },error=function(e){message(paste(pt,"tile is missing"))})
    }
  }
  dat = do.call(rbind.data.frame, data)
  flpt = do.call(rbind.data.frame, strsplit(rownames(dat), split = "/"))[5]
  flpt = do.call(rbind.data.frame, strsplit(unlist(flpt), split = "_"))[1]
  colnames(flpt) = "flpt"
  dat = cbind.data.frame(flpt, dat)
  dat["band_site"] = "TALL"
  rownames(dat)=NULL
  readr::write_csv(final, "./indir/spectra/silva_save_15.csv")
}

extract_spectra_from_itcs(itcs, raster_path)
