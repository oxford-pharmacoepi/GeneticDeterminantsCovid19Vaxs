getManhattanPlot <- function(gwas, 
                             y_lim,
                             dot_size = 0.1,
                             colors = c("#92C5DE","#4393C3","#2166AC"),
                             chr_len = 22,
                             reduce_dataset = 0.01,
                             h_line_color = "red"
){
  # Reduce dataset
  gwas_top <- gwas %>% filter(LOG10P > 1)
  gwas_low <- gwas %>% filter(LOG10P <= 1) %>% sample_frac(reduce_dataset)
  
  gwas1 <- gwas_top %>% full_join(gwas_low) 
  
  # Manhattan plot -------
  don <- gwas1 %>%
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwas1, ., by = c("CHR" ="CHR")) %>%
    # Add the cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum = BP+tot)
  
  # x axis
  axisdf = don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 )
  
  plot1 <- ggplot(don, aes(x=BPcum, y=LOG10P)) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), size = dot_size) +
    scale_color_manual(values = rep(colors, chr_len)) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center, expand = c(0.015,0.015)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,y_lim,2)) +     # remove space between plot area and x axis
    coord_cartesian(ylim = c(0,y_lim)) +
    # # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.text.x = element_text(margin = margin(t = 1), size = 5),
      axis.text.y  = element_text(size = 6),
      axis.title = element_text(size = 7),
    ) +
    #Plot a red horizontal line at 5e-8
    geom_hline(yintercept = -log10(5e-8), color = h_line_color, linewidth = 0.3) +
    labs(x = 'Chromosome', y = expression(-log[10](P)))
  
  return(plot1)
}

