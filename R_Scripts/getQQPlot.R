getQQPlot <- function(gwas,x_lim,y_lim,
                      color = "#2166AC",
                      dot_size = 0.1,
                      reduce_dataset = 0.01
){
  
  # QQ plot ----------------------------------------------------------------------
  n <- nrow(gwas)
  ci <- .95
  
  dat <- data.frame(
    observed = sort(gwas$LOG10P),
    expected = sort(-log10(ppoints(n))),
    clower   = sort(-log10(qbeta(p = (1 - ci) / 2, shape1 = seq(n), shape2 = rev(seq(n))))),
    cupper   = sort(-log10(qbeta(p = (1 + ci) / 2, shape1 = seq(n), shape2 = rev(seq(n)))))
  )
  
  # Reduce dataset size
  data_top <- dat %>% filter(observed > 1)
  data_low <- dat %>% filter(observed <= 1) %>% sample_frac(reduce_dataset)
  dat <- data_top %>% full_join(data_low)
  
  # Customize qqplot
  plot2 <- ggplot(dat, aes(x = expected, y = observed)) +
    scale_x_continuous(limits = c(0,x_lim), expand = c(0,0), breaks = seq(0,7,2)) + 
    scale_y_continuous(limits = c(0,y_lim), expand = c(0,0), breaks = seq(0,y_lim,4)) +
    geom_segment(data = . %>% filter(expected == max(expected)),
                 aes(x = 0, xend = x_lim, y = 0, yend = x_lim),
                 size = 0.25, color = "grey30", lineend = "round",alpha = 0.7) +
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    geom_point(color = color,  size = dot_size) +
    labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
         y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_rect(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title = element_text(size = 5),
      axis.text  = element_text(size = 5)
    ) 
  
  return(plot2)
}
