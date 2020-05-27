library(sf)

#link polygons to field_Data
field_data = readr::read_csv("~/Documents/Data/Dimensions/trait database/plant_info_2018_0628_dimensions_OSBS_TALL_MLBS.csv")
trait_data= readxl::read_excel("~/Documents/Data/Dimensions/trait database/leaf_traits_2018_0628_dimensions_OSBS_TALL_MLBS.xlsx")

#load_OSBS
OSBS = read_sf("~/Documents/Data/Dimensions/OSBS_crown_polygons/OSBS_sample_polygons_edits_Feb2018.shp")
TALL = read_sf("~/Documents/Data/Dimensions/TALL_crown_polygons/TALL_sample_crowns_edits_May2018.shp")
MLBS = read_sf("~/Documents/Data/Dimensions/MLBS_crown_polygons/MLBS_sample_polygons.shp")

field_data = field_data %>% select(ID, SITE, CRLIGHT, GENUS, SPECIES, LAT, LON, dbh1)
field_data = field_data %>% filter(!is.na(LAT))
field_data = st_as_sf(field_data, coords = c("LON", "LAT"), crs = 4326)
traitDB = inner_join(field_data, trait_data) %>% unique
foo = traitDB %>% filter(SITE == "OSBS")
foo = st_transform(foo, crs = st_crs(OSBS))
OSBS_crwn = sf::st_join(OSBS, foo, join = st_contains)
st_write(OSBS_crwn, "~/Documents/Data/Dimensions/OSBS.shp", delete_layer=TRUE)

foo = traitDB %>% filter(SITE == "TALL")
foo = st_transform(foo, crs = st_crs(TALL))
TALL_crwn = sf::st_join(foo, TALL)
st_write(TALL_crwn, "~/Documents/Data/Dimensions/TALL.shp", delete_layer=TRUE)

foo = traitDB %>% filter(SITE == "MLBS")
foo = st_transform(foo, crs = st_crs(MLBS))
MLBS_crwn = sf::st_join(foo, MLBS)
st_write(MLBS_crwn, "~/Documents/Data/Dimensions/MLBS.shp", delete_layer=TRUE)

OSBS_crwn = OSBS_crwn %>% select(SITE, CRLIGHT, GENUS, SPECIES, DATE, geometry,
                                 LMA, SLA, Ppercent, Parea, Npercent, Narea,
                                 Cpercent, Carea, d15Npermil,d13Cpermil, ID.x)

TALL_crwn = TALL_crwn %>% select(SITE, CRLIGHT, GENUS, SPECIES, DATE, geometry,
                                 LMA, SLA, Ppercent, Parea, Npercent, Narea,
                                 Cpercent, Carea, d15Npermil,d13Cpermil, tree)
colnames(OSBS_crwn)[17] = "individualID"
colnames(TALL_crwn)[17] = "individualID"
final = rbind.data.frame(OSBS_crwn, TALL_crwn)
final$taxonID = paste(final$GENUS, final$SPECIES, sep = "_")
write_csv(final, "~/Documents/GitHub/hiPyRneon/indir/Traits/Chapter1_field_data.csv")
