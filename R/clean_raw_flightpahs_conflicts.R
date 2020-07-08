remove_outliers <- function(spectra, method = "lof"){
  library(bigutilsr)
  pca <- prcomp(spectra[-1], scale. = TRUE)
  U <- pca$x
  # library(ggplot2)
  # theme_set(bigstatsr::theme_bigstatsr(0.8))
  # qplot(U[, 1], U[, 2]) + coord_equal()
  #
  if(method == "pca_dist"){
    ids <- apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 )) %>%
      Reduce(union, .)
  }
  if(method == "lof"){
    llof <- LOF(U)
    div_lim = hist_out(llof)$lim[2]
    ids = which(llof > div_lim)
  }

  eigval <- pca$sdev^2
  npc = pca_nspike(eigval)
  return(list(outliers = ids, npc = npc))
}

clean_data = list()

  #plot reflectances
  id = spectra$individualID %>% unique
  ii = ii+1
  plot_data <- spectra %>%
    dplyr::select(-one_of(c("site_ID", "species_ID",  "band_site","band_species", "flightpath"))) %>%
    dplyr::filter(individualID == id[ii])
  outlrs = remove_outliers(plot_data[2:370])
  if(length(outlrs$outliers) > 0){
    plot_data = plot_data[-outlrs$outliers,]
  }
  plot_spectra(plot_data)
  }
  #from large to long
plot_spectra<-function(plot_data){
  plot_data <- plot_data %>%
    dplyr::select(-one_of(c( "band_site","band_species", "flightpath")))
  ggdat = reshape2::melt(data.frame(plot_data), id.vars = c("individualID","flpt" ))
  ggp <- ggplot(ggdat, aes(x = variable, y = value)) +
    geom_line(aes(color = factor(flpt), alpha= 1), size = 0.2) +
    theme_bw()+
    theme(legend.position="none")
  ggp
  return(ggp)

}
