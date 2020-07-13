clean_data = list()
plot_spectra<-function(spectra){
  #plot reflectances
  id = spectra$individualID %>% unique
  id = c(146, 192, 238, 242, 260, 326, 328, 331, 339, 343, 349, 350, 357, 365, 368, 2108, 2118, 2119, 2128, 2130, 2131, 2132, 2137, 2151, 2159, 2165, 2174, 2178, 2180, 2200, 2202, 2205, 2242)
  ii = ii+1
  plot_data <- spectra %>%
    dplyr::select(-one_of(c("site_ID", "species_ID",  "band_site","band_species", "flightpath"))) %>%
    dplyr::filter(individualID == id[ii])
  outlrs = remove_outliers(plot_data)
  if(length(outlrs$outliers) > 0){
    plot_data = plot_data[-outlrs$outliers,]
  }
  #clean_data[[ii]] = plot_data
  plot_data <- plot_data[-1] %>%
    t %>%
    data.frame
  colnames(plot_data) = id[ii] #unlist(spectra[1]) # the first row will be the header
  plot_data <- data.frame(bnd = 1:dim(plot_data)[1], plot_data)
  ggdat <- tidyr::gather(plot_data, treeID, Reflectance,-bnd)
  ggp <- ggplot(ggdat, aes(x = bnd, y = Reflectance)) +
    geom_line(aes(color = factor(treeID), alpha= 1), size = 0.2) +
    theme_bw()+
    theme(legend.position="none")
  ggp
  return(ggp)

}


library(ForestTools)

pt = "////orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L3/LiDAR/CHM/"
get_tiles = list.files(pt, pattern = ".tif")
chm = raster::raster("~/Documents/2014_OSBS_1_412000_3280000_CHM.tif")
lin <- function(x){x * 0.15 + 0.6}

chm = raster::raster(paste(pt,r, sep="/"))
ttops <- vwf(CHM = chm, winFun = lin, minHeight = 2)
# tops = sf::st_as_sf(ttops)
# ttops= sf::write_sf(tops, "~/Documents/2014_OSBS.shp")
# tops = sf::read_sf("////orange/ewhite/s.marconi/Chapter1/hiPyRneon/indir/itcs/TALL_tops.shp")
preds = lidR::silva2016(chm, ttops)()
writeRaster(preds, "~/Documents/testITC.tif")


# check which flightpath is missing::
used_pts = list.files("~/Dropbox (UFL)/Chapter1_results/OSBS/Npercent/", ".tif")
missing_pts  = list.files("~/Documents/Data/ITCS/OSBS/ITCs/", ".shp")
itcs_pt = do.call(rbind.data.frame, strsplit(missing_pts, split = "_"))[3]
itcs_pt = unlist(itcs_pt)  %>% unique
used_pts = do.call(rbind.data.frame, strsplit(used_pts, split = "_"))[3]
used_pts = unlist(used_pts) %>% unique



library(sf)
osbs = data.table::fread("/Volumes/Athena/BackUps/May2020/Data/FULL_TOH.csv")
tall = osbs %>% filter(Site == "TALL") %>% select(ID, CHM, Slope, Albedo, Aspect, CA, geometry)
osbs2 = osbs %>% filter(Site == "OSBS") %>% select(ID, CHM, Slope, Albedo, Aspect, CA, geometry)
write_csv(osbs2, "~/Documents/Data/OSBS.csv")
write_csv(tall, "~/Documents/Data/TALL.csv")
