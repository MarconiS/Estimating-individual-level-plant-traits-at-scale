#apply silva 2016
library("lidR")
library(ForestTools)
library(raster)
colnames(test_result$y_95$ceam)[2:14] = c("y_up","siteID", "taxonID", "LMA", "SLA","Ppercent", "Parea",
                                          "Npercent","Narea","Cpercent","Carea","d15Npermil", "d13Cpermil"  )
colnames(test_result$y_5$ceam)[2:14] = c("y_lw","siteID", "taxonID", "LMA", "SLA","Ppercent", "Parea",
                                          "Npercent","Narea","Cpercent","Carea","d15Npermil", "d13Cpermil"  )
up = test_result$y_95$ceam$y_up > test_result$y_95$ceam$Npercent
lw = test_result$y_5$ceam$y_lw < test_result$y_5$ceam$Npercent
sum(up*lw)/length(up)

up = test_result$y_95$ceam$y_up < test_result$y_95$ceam$Npercent
lw = test_result$y_5$ceam$y_up > test_result$y_95$ceam$Npercent
sum(up*lw)/nrows(up)


colnames(test_result$y_95$epbm)[2:14] = c("y_up","siteID", "taxonID", "LMA", "SLA","Ppercent", "Parea",
                                          "Npercent","Narea","Cpercent","Carea","d15Npermil", "d13Cpermil"  )
colnames(test_result$y_5$epbm)[2:14] = c("y_lw","siteID", "taxonID", "LMA", "SLA","Ppercent", "Parea",
                                         "Npercent","Narea","Cpercent","Carea","d15Npermil", "d13Cpermil"  )
up = test_result$y_95$epbm$y_up > test_result$y_95$epbm$Npercent
lw = test_result$y_5$epbm$y_lw < test_result$y_5$epbm$Npercent
sum(up*lw)/length(up)

up = test_result$y_95$epbm$y_up < test_result$y_95$epbm$Npercent
lw = test_result$y_5$epbm$y_up > test_result$y_95$epbm$Npercent
sum(up*lw)/length(up)

library(sf)
chm = raster::raster("/Users/sergiomarconi/OneDrive - University of Florida/Site_rasters/OSBS_CHM.tif")
#get treetops
lin <- function(x){x * 0.05 + 0.6}
ttops <- vwf(CHM = chm, winFun = lin, minHeight = 2)
saveRDS(ttops, "./osbs_ttops.rds")
osbs = silva2016(chm, ttops, max_cr_factor = 0.6, exclusion = 0.3, ID = "treeID")()
