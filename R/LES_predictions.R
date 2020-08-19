osbs = sf::read_sf("/Users/sergiomarconi/Documents/traits_on_heaven/all_but_ele/OSBS_full_centr.shp")
globnet <- readr::read_csv("../toh/globnet.csv")[2:4]
osbs_elevation = sf::read_sf("/Users/sergiomarconi/Documents/traits_on_heaven/final_osbs.shp", as_tibble = T)
osbs_elevation = data.frame(osbs_elevation)
tall_elevation = sf::read_sf("/Users/sergiomarconi/Documents/traits_on_heaven/finaltall.shp", as_tibble = T)
tall_elevation = data.frame(tall_elevation)
osbs$geometry = as.character(osbs$geometry)
osbs_elevation$geometry = as.character(osbs_elevation$geometry)
foo = left_join(osbs, osbs_elevation)
write_csv(foo, "~/Documents/Data/fOSBS.csv")

tall$geometry = as.character(tall$geometry)
tall_elevation$geometry = as.character(tall_elevation$geometry)
foo2 = left_join(tall, tall_elevation)
write_csv(foo2, "~/Documents/Data/fTALL.csv")

foo = foo %>% dplyr::select(DBH, CHM, CA, LMA, C, N, P, Slope, Aspect, Albedo, OSBS_eleva)
foo2 = foo2 %>% dplyr::select(DBH, CHM, CA, LMA, C, N, P, Slope, Aspect, Albedo, TALL_eleva)
colnames(foo)[11] = "Elev"
colnames(foo2)[11] = "Elev"

foo$Slope = (1-cos(foo$Slope)) + (1-sin(foo$Slope))
foo$Aspect = (1-cos(foo$Aspect)) + (1-sin(foo$Aspect))
foo2$Aspect = (1-cos(foo2$Aspect)) + (1-sin(foo2$Aspect))
foo2$Slope = (1-cos(foo2$Slope)) + (1-sin(foo2$Slope))


corr <- round(cor(fcorr[complete.cases(fcorr),]), 2)
ggcorrplot(corr, hc.order = F, type = "lower", lab = TRUE)

library(propagate)

tall = osbs
dat = osbs
plot_les_comparison <- function(dat){
  library(ggpubr)
  globnet <- readr::read_csv("/Users/sergiomarconi/Documents/GitHub/toh/glopnet.csv")[2:4]
  dat <- dat %>%
    dplyr::select("LMA", "N", "P")
  dat = data.frame(dat) %>% dplyr::select(-one_of("geometry"))
  dat[c("LMA", "N", "P")] =   log10(dat[c("LMA", "N", "P")])
  dat$source <- "predictions"
  globnet= log10(globnet)
  globnet$source <- "globnet"
  # quantile(globnet[3],probs=c(.025,.975), na.rm=T)
  # range(dat$P, na.rm=T)
  # range(globnet[3], na.rm=T)
  colnames(globnet) <- c("LMA", "N", "P", "source")
  globnet = globnet %>%filter(P !=0)
  out_data <- dat %>%
    rbind(globnet)

  lma_N_plot <- ggscatter(out_data, x = "LMA", y = "N", add = "reg.line",
                          conf.int = TRUE, color = "source", palette = "jco", alpha = 0.7,
                          shape = "source") + ylim(-1,1) + xlim(0.2, 4)#stat_cor(aes(color = source))
  P_N_plot <- ggscatter(out_data, x = "P", y = "N", add = "reg.line",
                        conf.int = TRUE, color = "source", palette = "jco", alpha = 0.7,
                        shape = "source") + stat_cor(aes(color = source))
  lma_P_plot <- ggscatter(out_data, x = "LMA", y = "P", add = "reg.line",
                          conf.int = TRUE, color = "source", palette = "jco", alpha = 0.7,
                          shape = "source") + xlim(0.2, 4)#stat_cor(aes(color = source), label.x = 3)
  ggarrange(lma_N_plot, P_N_plot, lma_P_plot + rremove("x.text"),
            labels = c("A", "B", "C"),
            ncol = 3, nrow = 1)
}

plot(globnet$`log LMA`, globnet$`log Nmass`)
plot(table$LMA, table$`%N`)
plot(globnet$`log LMA`, globnet$`log Pmass`)
plot(globnet$`log Pmass`, globnet$`log Nmass`)

library(tidyverse)
fit1 <- lm(globnet$`log LMA` ~  globnet$`log Nmass`)
fit2 <- lm(log(table$LMA) ~ log(table$`%N`))

ggplot(globnet, aes(x = `log LMA`, y = `log Nmass`)) +
  stat_smooth(method = "lm", col = "red")

ggplot(table, aes(x = log(LMA), y = log(`%N`))) +
  stat_smooth(method = "lm", col = "red")



