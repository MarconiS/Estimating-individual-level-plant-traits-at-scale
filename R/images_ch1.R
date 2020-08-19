FULL <- readr::read_csv("~/Downloads/FULL_TOH.csv")
head(FULL)
FULL[FULL==-9999] = NA
library(tidyverse)
covariance_plots <- function(dat, title){
  #dat <- cleantableCrown
  library(corrplot)
  dat <- dat[complete.cases(dat), ]
  M <- cor(dat)
  diag(M) = NA
  corrplot(M, type = "lower", title = title, na.label = ".", tl.srt = 45)
}


a <- FULL %>% select(LMAg_m2, Npct, Ppct, CA, CHM, LAI, Albedo, DTM, Aspect, Slope)
colnames(a)[c(1,2,3,4,5,8)] <- c("LMA", "%N", "%P", "Cr Area", "Height", "Elevation")
a%>% covariance_plots(title = "")

a <- FULL %>% filter(Site == "OSBS") %>%
  select(LMAg_m2, Npct, Ppct, CA, CHM, LAI, Albedo, DTM, Aspect, Slope)
colnames(a)[c(1,2,3,4,5,8)] <- c("LMA", "%N", "%P", "Cr Area", "Height", "Elevation")
a%>% covariance_plots(title = "OSBS")

a <- FULL %>% filter(Site == "TALL") %>%
  select(LMAg_m2, Npct, Ppct, CA, CHM, LAI, Albedo, DTM, Aspect, Slope)
colnames(a)[c(1,2,3,4,5,8)] <- c("LMA", "%N", "%P", "Cr Area", "Height", "Elevation")
a%>% covariance_plots(title = "TALL")


covariance_plots <- function(dat, title){
  #dat <- cleantableCrown
  library(corrplot)
  dat <- dat[complete.cases(dat), ]
  M <- cor(dat)
  diag(M) = NA
  corrplot(M, type = "lower", title = title, na.label = ".", tl.srt = 45)
}


a <- FULL %>% select(LMAg_m2, Npct, Ppct, CA, CHM, LAI, Albedo, DTM, Aspect, Slope)
colnames(a)[c(1,2,3,4,5,8)] <- c("LMA", "%N", "%P", "Cr Area", "Height", "Elevation")
a%>% covariance_plots(title = "")

a <- FULL %>% filter(Site == "OSBS") %>%
  select(LMAg_m2, Npct, Ppct, CA, CHM, LAI, Albedo, DTM, Aspect, Slope)
colnames(a)[c(1,2,3,4,5,8)] <- c("LMA", "%N", "%P", "Cr Area", "Height", "Elevation")
a%>% covariance_plots(title = "OSBS")

a <- FULL %>% filter(Site == "TALL") %>%
  select(LMAg_m2, Npct, Ppct, CA, CHM, LAI, Albedo, DTM, Aspect, Slope)
colnames(a)[c(1,2,3,4,5,8)] <- c("LMA", "%N", "%P", "Cr Area", "Height", "Elevation")
a%>% covariance_plots(title = "TALL")

plot_les_comparison <- function(dat, ssize){
  library(ggpubr)
  globnet <- readr::read_csv("~/Downloads/globnet.csv")[2:4]
  colnames(globnet) <- c("LMA", "N", "P")

  dat <- dat %>%
    dplyr::select("LMA", "N", "P") %>%
    mutate_all(log10) %>%
    filter(LMA > (min(globnet$LMA, na.rm = T)-0.5)) %>%
    filter(N > (min(globnet$N, na.rm = T)-0.5)) %>%
    filter(P > (min(globnet$P, na.rm = T)-0.5)) %>%
    sample_n(ssize)

  dat$source <- "ITCs"
  globnet$source <- "LES"
  colnames(globnet) <- c("LMA", "N", "P", "source")
  out_data <- dat %>%
    rbind(globnet)

  lma_N_plot <- ggscatter(out_data, x = "LMA", y = "N", add = "reg.line",
                          conf.int = TRUE, color = "source", palette = "jco", alpha = 0.3,
                          shape = "source") + stat_cor(aes(color = source), label.x = 2.2) +theme_bw()#+ geom_density2d(colour=c("black", "gray"))
  P_N_plot <- ggscatter(out_data, x = "P", y = "N", add = "reg.line",
                        conf.int = TRUE, color = "source", palette = "jco", alpha = 0.3,
                        shape = "source") + stat_cor(aes(color = source)) +theme_bw()#+ geom_density2d(colour=c("black", "gray"))
  lma_P_plot <- ggscatter(out_data, x = "LMA", y = "P", add = "reg.line",
                          conf.int = TRUE, color = "source", palette = "jco", alpha = 0.3,
                          shape = "source") + stat_cor(aes(color = source)) +theme_bw()#+ geom_density2d(colour= c("black", "gray"))
  ggarrange(lma_N_plot, P_N_plot, lma_P_plot + rremove("x.text"),
            labels = c("A", "B", "C"),
            ncol = 3, nrow = 1)
}

dat <- FULL %>% select(LMAg_m2, Npct, Ppct)
colnames(dat) <- c("LMA", "N", "P")
plot_les_comparison(dat, 100000)

dat <- FULL %>% filter(Site == "OSBS") %>% select(LMAg_m2, Npct, Ppct)
colnames(dat) <- c("LMA", "N", "P")
plot_les_comparison(dat, 100000)

dat <- FULL %>% filter(Site == "TALL") %>%  select(LMAg_m2, Npct, Ppct)
colnames(dat) <- c("LMA", "N", "P")
plot_les_comparison(dat, 100000)

unmixed_distributions <- function(dat, tr_nm, modes = 2, ssize = 100){
  library(mixtools)
  dat <- dat %>% data.frame %>% sample_n(ssize) %>% as.matrix
  # mixt<-normalmixEM(dat,k=modes,fast=TRUE, maxit = 10000)
  # dat2gauss <- data.frame(value = rnorm(length(dat)/(modes*5),
  #                     mean = mixt$mu[1], sd = mixt$sigma[1]),
  #                      tl = paste("UnMix", "1", sep=""))
  # for(mx in 2:modes){
  #   tmp <- data.frame(value = rnorm(length(dat)/(modes*5),
  #                            mean = mixt$mu[mx], sd = mixt$sigma[mx]),
  #                     tl = paste("UnMix", mx, sep=""))
  #   dat2gauss <- rbind(dat2gauss, tmp)
  # }
  # #pi<-mixt$lambda[1]
  # # mu1<-mixt$mu[1]
  # # mu2<-mixt$mu[2]
  # # sigma1<-mixt$sigma[1]
  # # sigma2<-mixt$sigma[2]
  # #
  # # gauss1 <- data.frame( value = rnorm(10000, mean = mu1, sd = sigma1), tl = "UnMix1")
  # # gauss2 <- data.frame( value = rnorm(10000, mean = mu2, sd = sigma2), tl = "UnMix2")
  # # dat2gauss = rbind(gauss1, gauss2) %>% data.frame
  # dat2gauss$tl <- factor(dat2gauss$tl)
  # colnames(dat2gauss) <- c(tr_nm, "mixture")
  dat <- data.frame(dat, "MixD")
  colnames(dat) <- c(tr_nm, "Site")
  #dat2gauss <- rbind(dat, dat2gauss)
  dat2gauss <- dat
  ggdensity(dat2gauss, x = tr_nm, y = "..count..",
            add = "mean", rug = TRUE,
            color = "mixture", palette = "jco", fill = "mixture", alpha = 0.5) +
    theme_bw() +theme(legend.position="none") + labs(x = tr_nm, y = "n ITCs")
}

tr_istributions <- function(dat, ssize = 300000, full=F){
  library(ggpubr)
  library(ggsci)

  dat <- dat %>% sample_n(ssize) %>% filter(CA > 4)
  dat = dat %>% filter(Height> 0)
  dat = dat[complete.cases(dat), ]
  if(full ==T){
    dt2 <- dat
    dt2$Site <- "All"
    dat = rbind(dt2, dat)
  }
  lma <- ggdensity(dat, x = "LMA", y = "..count..",
                   add = "mean", rug = TRUE,
                   color = "Site", palette = "jco", fill = "Site", alpha = 0.5) +
    theme_bw()  + labs(x = "LMA", y = "n ITCs") +theme(legend.position="none") +
    xlim(-50,500)

  lma_h <- dat %>%  sample_n(10000) %>% ggplot(aes(x = LMA, y = Height)) +
    stat_ellipse(geom = "polygon",  alpha = 0.3, aes(color = Site, fill = Site))  + xlim(-50,500) +
    scale_fill_manual(values = c('#99CCFF', '#FCE516')) +
    geom_smooth(method = "lm", aes(color = Site)) +
    theme(legend.position="none") + theme_bw() + scale_color_jco()  +theme(legend.position="none")

  ggarrange(lma, lma_h,
            ncol = 1, nrow = 2)

  n <- ggdensity(dat, x = "N", y = "..count..",
                 add = "mean", rug = TRUE,
                 color = "Site", palette = "jco", fill = "Site", alpha = 0.5) + xlim(0,4) +
    theme_bw() + labs(x = "%N", y = "n ITCs") + theme(legend.position="none")

  n_h <- dat %>%  sample_n(10000) %>% ggplot(aes(x = N, y = Height)) +
    stat_ellipse(geom = "polygon",  alpha = 0.3, aes(color = Site, fill = Site))  + xlim(0,4) +
    scale_fill_manual(values = c('#99CCFF', '#FCE516')) +
    geom_smooth(method = "lm", aes(color = Site)) +
    theme(legend.position="none") + theme_bw() + scale_color_jco()  +theme(legend.position="none")

  ggarrange(n, n_h,
            ncol = 1, nrow = 2)

  p <- ggdensity(dat, x = "P", y = "..count..",
                 add = "mean", rug = TRUE,
                 color = "Site", palette = "jco", fill = "Site", alpha = 0.5) + xlim(0,0.2) +
    theme_bw()  + labs(x = "%P", y = "n ITCs")  + theme(legend.position="none")

  p_h <- dat %>%  sample_n(10000) %>% ggplot(aes(x = P, y = Height)) +
    stat_ellipse(geom = "polygon",  alpha = 0.3, aes(color = Site, fill = Site))  + xlim(0,0.2) +
    scale_fill_manual(values = c('#99CCFF', '#FCE516')) +
    geom_smooth(method = "lm", aes(color = Site)) +
    theme(legend.position="none") + theme_bw() + scale_color_jco()  +theme(legend.position="none")
  ggarrange(p, p_h,
            ncol = 1, nrow = 2)

  c <- ggdensity(dat, x = "C", y = "..count..",
                 add = "mean", rug = TRUE,
                 color = "Site", palette = "jco", fill = "Site", alpha = 0.5) + xlim(47,54) +
    theme_bw()  + labs(x = "%C", y = "n ITCs")  + theme(legend.position="none")

  c_h <- dat %>%  sample_n(10000) %>% ggplot(aes(x = C, y = Height)) +
    stat_ellipse(geom = "polygon",  alpha = 0.3, aes(color = Site, fill = Site))  +xlim(47,54) +
    scale_fill_manual(values = c('#99CCFF', '#FCE516')) +
    geom_smooth(method = "lm", aes(color = Site)) +
    theme(legend.position="none") + theme_bw() + scale_color_jco()  +theme(legend.position="none")
  ggarrange(p, p_h,
            ncol = 1, nrow = 2)

  ggarrange(lma, n, p, c,
            labels = c("A", "B", "C", "D"),
            ncol = 4, nrow = 1)


  #split by "BLF" "NDL" using an hard cut (best a mixed model, but ok)
  dat$PFT = "BDL"
  dat$PFT[dat$C > 51] = "NDL"

  ggplot(dat, aes(x = LMA, y = Height, group = as.factor(PFT))) +
    #geom_bin2d() +
    stat_density_2d(geom = "polygon",aes(alpha = ..level.., fill = as.factor(PFT)), bins = 6) +
    facet_grid(. ~ Site) +   theme(legend.position="none") + theme_bw() #geom_smooth(method = "lm", aes(colour = PFT)) +

  ggplot(dat, aes(x = C, y = Height, group = as.factor(PFT))) +
    #geom_bin2d() +
    stat_density_2d(geom = "polygon",aes(alpha = ..level.., fill = as.factor(PFT)), bins = 6) +
    facet_grid(. ~ Site) +   theme(legend.position="none") + theme_bw() #geom_smooth(method = "lm", aes(colour = PFT)) +

  ggplot(dat, aes(x = P, y = Height, group = as.factor(PFT))) +
    #geom_bin2d() +
    stat_density_2d(geom = "polygon",aes(alpha = ..level.., fill = as.factor(PFT)), bins = 6) +
    facet_grid(. ~ Site) +   theme(legend.position="none") + theme_bw() #geom_smooth(method = "lm", aes(colour = PFT)) +

  ggplot(dat, aes(x = N, y = Height, group = as.factor(PFT))) +
    #geom_bin2d() +
    stat_density_2d(geom = "polygon",aes(alpha = ..level.., fill = as.factor(PFT)), bins = 6) +
    facet_grid(. ~ Site) +   theme(legend.position="none") + theme_bw() #geom_smooth(method = "lm", aes(colour = PFT)) +

  dat %>% sample_n(1000) %>% ggplot(aes(x = LMA, y = Height)) +
    geom_smooth(method = "lm",aes(color = "Site", palette = "jco")) + theme(legend.position="none") + theme_bw()

  library("ggsci")
  dat %>%  sample_n(10000) %>% ggplot(aes(x = LMA, y = Height)) +
    stat_ellipse(geom = "polygon",  alpha = 0.3, aes(color = PFT, fill = PFT))  +
    scale_fill_manual(values = c("blue", "yellow")) +
    geom_smooth(aes(color = Site)) + theme(legend.position="none") + theme_bw() +
    theme(legend.position="none") + theme_bw() + scale_color_jco()

  dat %>%  sample_n(10000) %>% ggplot(aes(x = P, y = Height)) +
    stat_ellipse(geom = "polygon",  alpha = 0.3, aes(color = Site, fill = Site))  + ylim +
    scale_fill_manual(values = c("blue", "yellow")) +
    geom_smooth(method = "lm", aes(color = Site)) +
    theme(legend.position="none") + theme_bw() + scale_color_jco()
  # + scale_fill_manual(values = c("yellow", "blue"))

}
library(ggpubr)

dat <- fcorr %>% dplyr::select(LMA, N, P, C, CHM, CA, SITE)
colnames(dat) <- c("LMA", "N", "P", "C", "Height", "CA", "Site")
tr_istributions(dat = dat, ssize = 150000)
tr_istributions(dat = dat, ssize = 300000, full = T)
ggsave("./tr_distributions.jpg", width = 24, height = 9, units = "cm")

dat <- FULL %>% filter(Site == "OSBS") %>% select(LMAg_m2, Npct, Ppct)
colnames(dat) <- c("LMA", "N", "P")
lma <- unmixed_distributions(dat = dat$LMA, modes = 2, tr_nm = "LMA",ssize = 100000)
n <- unmixed_distributions(dat = dat$N, modes = 2, tr_nm = "N",ssize = 100000)
p <- unmixed_distributions(dat = dat$P, modes = 2, tr_nm = "P",ssize = 100000)

ggarrange(lma, n, p + rremove("x.text"),
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

dat <- FULL %>% filter(Site == "TALL") %>% select(LMAg_m2, Npct, Ppct)
colnames(dat) <- c("LMA", "N", "P")
lma <- unmixed_distributions(dat = dat$LMA, modes = 2, tr_nm = "LMA",ssize = 100000)
n <- unmixed_distributions(dat = dat$N, modes = 2, tr_nm = "N",ssize = 100000)
p <- unmixed_distributions(dat = dat$P, modes = 2, tr_nm = "P",ssize = 100000)

ggarrange(lma, n, p + rremove("x.text"),
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


sample_from_unc <- function(dat){
  sampled <- matrix(nrow = dim(dat)[1], ncol = 3)

  for(bb in 1:dim(dat)[1]){
    lma_ave = dat$LMAg_m2[bb]
    lma_sdM = dat$LMA_95[bb]
    lma_sdm = dat$LMA_5[bb]
    uLMA = rnorm(1, mean = lma_ave, sd = (lma_sdM - lma_sdm)/4)
    #lma_unc2 = (rnorm(1, mean = lma_ave, sd = (lma_sdM - lma_ave)/2) + rnorm(1, mean = lma_ave, sd = (lma_ave - lma_sdm)/2))/2

    P_ave = dat$Ppct[bb]
    P_sdM = dat$P_95[bb]
    P_sdm = dat$P_5[bb]
    uP = rnorm(1, P_ave, sd = (P_sdM - P_sdm)/4)
    #P_unc2 = (rnorm(1, P_ave, sd = (P_sdM - P_ave)/2) + rnorm(1, P_ave, sd = (P_ave = P_sdm)/2))/2

    N_ave = dat$Npct[bb]
    N_sdM = dat$N_95[bb]
    N_sdm = dat$N_5[bb]
    uN = rnorm(1, N_ave, sd = (N_sdM - N_sdm)/4)
    #N_unc2 = (rnorm(1, N_ave, sd = (N_sdM - N_ave)/2) + rnorm(1, N_ave, sd = (N_ave = N_sdm)/2))/2

    sampled[bb,] <- c(uLMA, uN, uP)

  }
  return(as.data.frame(sampled))
}

dat <- FULL %>% select(LMAg_m2, LMA_95, LMA_5, Npct, N_95, N_5, Ppct, P_95, P_5)
#udat <- sample_from_unc(dat)


unc_propagation <- function(dat, ssize = 100){
  library(ggpubr)
  dat <- dat %>% data.frame %>% sample_n(ssize)
  udat <- sample_from_unc(dat)
  colnames(udat) <- c("LMAg_m2", "Npct", "Ppct")
  dat$source <- "point estimate"
  udat$source <- "propagation"

  tmp_dat <- dat %>% select(LMAg_m2, source)
  tmp_udat <- udat %>% select(LMAg_m2, source)
  dat2gauss <- rbind(tmp_dat, tmp_udat)
  lma <- ggdensity(dat2gauss, x = "LMAg_m2", y = "..count..",
                   add = "mean", rug = TRUE,
                   color = "source", palette = "jco", fill = "source", alpha = 0.5) +
    theme_bw() +theme(axis.title=element_blank(), legend.position="none")


  tmp_dat <- dat %>% select(Npct, source)
  tmp_udat <- udat %>% select(Npct, source)
  dat2gauss <- rbind(tmp_dat, tmp_udat)
  n <- ggdensity(dat2gauss, x = "Npct", y = "..count..",
                 add = "mean", rug = TRUE,
                 color = "source", palette = "jco", fill = "source", alpha = 0.5)+
    theme_bw() +theme(axis.title=element_blank(), legend.position="none")


  tmp_dat <- dat %>% select(Ppct, source)
  tmp_udat <- udat %>% select(Ppct, source)
  dat2gauss <- rbind(tmp_dat, tmp_udat)
  p <- ggdensity(dat2gauss, x = "Ppct", y = "..count..",
                 add = "mean", rug = TRUE,
                 color = "source", palette = "jco", fill = "source", alpha = 0.5) + xlim(0, 0.2) +
    theme_bw() +theme(axis.title=element_blank(), legend.position="none")

  ggarrange(lma, n, p, legend = "bottom",
            ncol = 3, nrow = 1)
}

dat <- FULL %>% select(LMAg_m2, LMA_95, LMA_5, Npct, N_95, N_5, Ppct, P_95, P_5) %>%
  unc_propagation(150000)
dat
crop <- FULL %>% dplyr::select(LMAg_m2, LMA_95, LMA_5, Npct, N_95, N_5, Ppct, P_95, P_5, lat, lon, geometry) %>%
  dplyr::filter(lat < -81.98997 & lat > -82.01288) %>%
  dplyr::filter(lon < 29.71288 & lon > 29.69686)

readr::write_csv(crop, "~/Desktop/crop_osbs.csv")
#tall
dat <- FULL %>% filter(Site == "TALL") %>% dplyr::select(LMAg_m2, LMA_95, LMA_5, Npct, N_95, N_5, Ppct, P_95, P_5) %>%
  unc_propagation(100000)
dat
#osbs
dat <- FULL %>% filter(Site == "OSBS") %>% dplyr::select(LMAg_m2, LMA_95, LMA_5, Npct, N_95, N_5, Ppct, P_95, P_5) %>%
  unc_propagation(150000)
dat

library(tidyverse)
library("ggsci")

r2_comparison <- readr::read_csv("~/Documents/GitHub/TraitsOnHeaven/misc/modelPerformance.csv")
deltar2 <- readr::read_csv("~/Documents/GitHub/TraitsOnHeaven/misc/deltaR2.csv")
r2_comparison <-  gather(r2_comparison, Trait, R2, `%P`:`%C`, factor_key=TRUE)
deltar2 <-  gather(deltar2, Trait, R2, `%P`:`%C`, factor_key=TRUE)

ggplot(r2_comparison, aes(x = Trait, y = R2, fill = Model)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + theme(legend.position="none") + theme_bw() +
  scale_fill_manual(values=c('#FCE516','#99CCFF', '#999999')) +
  labs(x=parse(text='Trait'), y=parse(text='R^2')) +   geom_hline(yintercept = 1, linetype = "dashed")

ggplot(deltar2, aes(x = Trait, y = R2, fill = Model)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + theme(legend.position="none") + theme_bw() +
  scale_fill_manual(values=c('#FCE516','#99CCFF', '#999999')) +
  labs(x=parse(text='Trait'), y=parse(text='dR^2')) +
  geom_hline(yintercept = 0,  linetype = "dashed")
