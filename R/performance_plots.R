library(tidyverse)
crown_ls <- list.files(path = out_dir, pattern = "crown_res.csv")
traits = c("LMA_g.m2", "N_pct", "P_pct", "C_pct")

for(ii in traits){
  crown_dat <- read.csv(paste(out_dir,
                              list.files(path = out_dir, pattern = paste(ii, "crown_res.csv", sep="_")), sep = ""))
  crown_dat$X = "ITC"
  pix_dat <- read.csv(paste(out_dir,
                            list.files(path = out_dir, pattern = paste(ii, "pix_res.csv", sep="_")), sep = ""))
  pix_dat$X = "pix"
  max_tr <-  max(max(crown_dat$y), max(crown_dat$yhat))
  min_tr <-  min(min(crown_dat$y), min(crown_dat$yhat))

  regr_dat <- rbind(pix_dat, crown_dat)
  p <- ggplot(data = regr_dat, aes(x= yhat, y = y, color=X)) + geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    #geom_smooth(method=lm,  linetype="dashed", se=FALSE) +
    xlim(min_tr, max_tr) + ylim(min_tr, max_tr) + theme_bw() +
    guides(fill=FALSE, color=FALSE) +
    labs(y = "Observed", x = "Predicted")
  ggsave(paste(out_dir, ii, "_oneOne.svg", sep=""), p, device = "svg")
}

#get coverage plots
for(ii in traits){
  crown_dat <- read.csv(paste(out_dir,ii, "_Coverage.csv", sep=""))

  cr_lw <- crown_dat %>% group_by(pixel_crownID, sim) %>% mutate(PI5, min(PI5))
  cr_up <- crown_dat %>% group_by(pixel_crownID, sim) %>% mutate(PI95, max(PI95))

  crown_dat$PI5 <- cr_lw$`min(PI5)`
  crown_dat$PI95 <- cr_up$`max(PI95)`
  p <- ggplot(data = crown_dat, aes(x= as.numeric(factor(pixel_crownID)), y = y)) + geom_point() +
    geom_ribbon(aes(ymin = PI5, ymax = PI95, color=sim),alpha=0.3) +
    theme_bw() +
    guides(fill=FALSE, color=FALSE) +
    labs(y = "Observed", x = "Held out Crown")
  p
  ggsave(paste(out_dir, ii, "_Coverage.svg", sep=""), p, device = "svg")
}



for(ii in traits){
  crown_dat <- read.csv(paste(out_dir,ii, "_Coverage.csv", sep=""))

  cr_lw <- crown_dat %>% group_by(sim) %>% mutate(PI5, min(PI5))
  cr_up <- crown_dat %>% group_by(sim) %>% mutate(PI95, max(PI95))

  crown_dat$PI5 <- cr_lw$`min(PI5)`
  crown_dat$PI95 <- cr_up$`max(PI95)`
  ggplot(data = crown_dat, aes(x= factor(sim), y = y)) +
    geom_errorbar(aes(ymin = PI5, ymax = PI95, color=sim)) + geom_point() +
    theme_minimal() +
    guides(fill=FALSE, color=FALSE) +
    labs(y = "Observed", x = "Crown ID")
  ggsave(paste(out_dir, ii, "_Coverage.jpg", sep=""), p, device = "svg")
}
#plot coverage
cover = list()
for(ii in traits){
  crown_dat <- read.csv(paste(out_dir,ii, "_Coverage.csv", sep=""))
  crown_dat$in_bnd <- (crown_dat$y < crown_dat$PI95) * (crown_dat$y > crown_dat$PI5)
  crown_dat <- crown_dat %>% group_by(sim) %>% mutate(cover = sum(in_bnd), n =n())
  cc <- crown_dat %>% select(sim, cover, n) %>% unique
  cover[[ii]] = cc$cover / cc$n
  #cover[[ii]] <- crown_dat$cover %>% unique
}

coverage <- do.call(cbind.data.frame, cover)
coverage$sim <- c("pixel", "crown")
coverage <- gather(coverage, trait, coverage, LMA_g.m2:C_pct, factor_key=TRUE)
coverage <- coverage %>% filter(Model != "EPBM")
ggplot(coverage, aes(x = factor(Trait), y = Coverage, color = Model)) + geom_point() +
  geom_hline(yintercept=0.95) + theme_linedraw()


