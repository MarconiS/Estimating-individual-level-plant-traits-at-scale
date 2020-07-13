#train calibration and test split
Y <- read.csv("./indir/Traits/Chapter1_field_data.csv")
X <- read.csv("./indir/Spectra/reflectance_all.csv")
Y_dat = Y %>% filter(individualID %in% X$individualID)
Y_dat = Y_dat %>% filter(CRLIGHT != "shade")
Y_dat = Y_dat %>% group_by(individualID, SITE, taxonID)%>%
  summarize_if(is.numeric, mean)

Y_train = Y_dat %>% group_by(taxonID, SITE) %>% sample_frac(0.7)
Y_test = Y_dat %>% filter(!individualID %in% Y_train$individualID)
Y_oob = Y_test %>% group_by(taxonID) %>% sample_frac(0.5)
Y_test = Y_test %>% filter(!individualID %in% Y_oob$individualID)

readr::write_csv(Y_train,"./indir/Misc/train_ids.csv")
readr::write_csv(Y_oob, "./indir/Misc/calibration_ids.csv")
readr::write_csv(Y_test, "./indir/Misc/oob_ids.csv")
