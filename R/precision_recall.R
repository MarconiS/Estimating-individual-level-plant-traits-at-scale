library(sf)
#filter by the crowns in the datset
test.data.x <- read.csv("./indir/Spectra/reflectance_all.csv")
train.data.y <- read.csv("~/Dropbox (UFL)/Chapter1_results/Chapter1_field_data.csv")
train.data.y = train.data.y %>% dplyr::filter(CRLIGHT != "shade")
ids = test.data.x %>% dplyr::filter(test.data.x$individualID %in% train.data.y$individualID)
ids = ids %>% dplyr::select(individualID)%>% unique

tall_itcs = read_sf("~/Documents/Data/Dimensions/TALL_crown_polygons/TALL_sample_crowns_edits_May2018.shp")
tall_itcs = tall_itcs %>% dplyr::filter(tall_itcs$tree %in% ids$individualID)
write_sf(tall_itcs, "~/Documents/Data/Dimensions/field_tall.shp")


tall_tops = data.table::fread("~/Documents/Data/TALL.csv")
tall_tops$DN = 1:nrow(tall_tops)

tall_tops = sf::st_as_sf(tall_tops, wkt = "geometry")
st_crs(tall_tops) = st_crs(tall_itcs)

tall_trees = sf::st_join(tall_itcs, tall_tops, left=T)
ids = tall_trees$DN %>% unique
tall_field = tall_tops %>% filter(DN %in% ids)
write_sf(tall_field, "~/Documents/Data/Dimensions/silva_tall.shp")




id_trees = sf::st_join(tall_tops, tall_itcs, left=T, join = st_within)
id_trees = id_trees %>% dplyr::filter(!is.na(tree))

missing = tall_itcs %>% dplyr::filter(!tall_itcs$tree %in% id_trees$tree)
undetected = nrow(missing)/nrow(tall_itcs)
write_sf(id_trees, "~/Documents/Data/TALL_ids_dt.shp")
write_sf(tall_tops, "~/Documents/Data/tall_tops.shp")


osbs_itcs = read_sf("~/Documents/Data/Dimensions/OSBS_crown_polygons/OSBS_sample_polygons_edits_Feb2018.shp")
osbs_itcs = osbs_itcs %>% dplyr::filter(osbs_itcs$ID %in% ids$individualID)
osbs_tops = sf::st_as_sf(osbs, wkt = "geometry")
st_crs(osbs_tops) = st_crs(osbs_itcs)

osbsid_trees = sf::st_join(osbs_itcs, osbs_tops, left=T)
ids = osbsid_trees$DN %>% unique
osbs_field = osbs_tops %>% filter(DN %in% ids)
write_sf(osbs_field, "~/Documents/Data/Dimensions/silva_osbs.shp")
osbsid_trees = osbsid_trees %>% dplyr::filter(!is.na(ID))

osbsmissing = osbs_itcs %>% dplyr::filter(!osbs_itcs$ID %in% osbsid_trees$ID)
osbsundetected = nrow(osbsmissing)/nrow(osbs_itcs)
write_sf(osbsid_trees, "~/Documents/Data/OSBS_ids.shp")
write_sf(osbs_tops, "~/Documents/Data/OSBS_tops.shp")
