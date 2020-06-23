library(sf)
toh = sf::read_sf("/Users/smarconi/Documents/Data/TOH/latlon_toh.shp")
osbs = toh %>% filter(Site == "OSBS")
field = sf::read_sf("/Users/smarconi/Documents/Data/Dimensions/OSBS.shp")
field = field %>% select(ID_x)
osbs = st_transform(osbs, crs = st_crs(field))
interosbs = st_join(osbs, field)
interosbs = interosbs %>% filter(!is.na(ID_x))
write_sf(interosbs, "~/Documents/Data/shp/silva_osbs.shp")

tall = toh %>% filter(Site == "TALL")
field = sf::read_sf("/Users/smarconi/Documents/Data/Dimensions/TALL_crown_polygons/TALL_sample_crowns_edits_May2018.shp")
field = field %>% select(tree)
tall = st_transform(tall, crs = st_crs(field))
interstall = st_join(tall, field)
interstall = interstall %>% filter(!is.na(tree))
write_sf(interstall, "~/Documents/Data/shp/silva_tall.shp")

#
interstall = st_transform(interstall, crs = st_crs(toh))
interosbs = st_transform(interosbs, crs = st_crs(field))

rbind(interstall, interosbs)
