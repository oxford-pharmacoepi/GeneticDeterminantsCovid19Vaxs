# Manhattan plot
dir <- c('oneDose',
         'twoDose',
         'breakthroughSusceptibility',
         'breakthroughSeverity')

tit <- c('A) Seroconversion - One dose', 'B) Seroconversion - Two doses',
         'C) Breakthrough susceptibility', 'D) Breakthrough severity')

ylim_MH <- c(30,22,70,19)
ylim_QQ <- c(14,18,18,14)
gws_y    <- c(-log10(1.5e-8),-log10(2.5e-8),-log10(3e-8),-log10(1e-8))
y_QQ_min <- c(15,11,37,10)
y_QQ_max <- c(32,25,80,21.75)
alph     <- 0.7
#colors <- c('#79A7A5','#A3C6BC','#CEE5D4')

for (ii in 1:length(dir)){
  gwas <-  as_tibble(read.table(paste0(dir_results,"GWAS/",dir[ii],"_assoc.regenie.merged.txt"), header = TRUE)) |>
    mutate(P = 10^(-LOG10P))  %>%
    select(CHR = CHROM, BP = GENPOS, SNP = ID, LOG10P, P) %>%
    filter(!is.na(LOG10P))
  
  source(here('R_Scripts',"getManhattanPlot.R"))
  plot1 <- getManhattanPlot(gwas,y_lim = ylim_MH[ii])
  
  source(here('R_Scripts',"getQQPlot.R"))
  plot2 <- getQQPlot(gwas, x_lim = 6.5, y_lim = ylim_QQ[ii])
  
  gwas_top <- gwas %>% filter(LOG10P > 1)
  gwas_low <- gwas %>% filter(LOG10P <= 1) %>% sample_frac(0)
  gwas1 <- gwas_top %>% full_join(gwas_low) 
  
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
  
  axisdf <- don %>%
    group_by(CHR) %>%
    summarize(center=(max(BPcum) + min(BPcum) ) / 2)
  
  if (ii == 1){
    plot1 <- plot1 +
      annotate('text',
               don %>% filter(CHR == 6) |> filter(LOG10P==max(LOG10P, na.rm = TRUE)) %>% select(BPcum) %>% as.numeric(),
               don %>% summarise(LOG10P=max(LOG10P)) %>% as.numeric()+1,
               label = 'rs9275109',
               size = 1.8, color = 'black', hjust=0.5, vjust = 0) +
      annotate('text',
               don |> filter(CHR == 2) |> filter(LOG10P==max(LOG10P, na.rm = TRUE)) %>% select(BPcum) %>% as.numeric(),
               don |> filter(CHR == 2) |> summarise(LOG10P=max(LOG10P)) %>% as.numeric()+1,
               label = 'rs79510369',
               size = 1.8, color = 'black', hjust=0.5, vjust = 0.5) +
      ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
  }else if(ii == 2){
    plot1 <- plot1 +
      annotate('text',
               don %>% filter(CHR == 6, BP == 32634226) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 6, BP == 32634226) |> summarise(LOG10P=max(LOG10P)) %>% as.numeric()+.25,
               label = 'rs68033958',
               size = 1.8, color = 'black', hjust=0.5, vjust = 0) +
      annotate('rect',
               xmin = don %>% filter(CHR == 6, BP == 32634226) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don %>% filter(CHR == 6, BP == 32634226) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = 9.5, ymax = 10.4,
               alpha = alph, fill = "white") +
      annotate('text',
               don %>% filter(CHR == 6, BP == 30332146) %>% select(BPcum) %>% as.numeric(),
               9.7,
               label = 'rs3094055',
               size = 1.8, color = 'black', hjust=0.5, vjust = 0) +
      ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
  }else if(ii == 3){
    plot1 <- plot1 +
      annotate('text',
               don %>% filter(CHR == 1, BP == 155123837) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 1, BP == 155123837) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs6676150',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',
               don %>% filter(CHR == 3, BP == 45835417) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 3, BP == 45835417) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+2.3,
               label = 'rs73062389',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('rect',
               xmin = don %>% filter(CHR == 3, BP == 101547733) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don %>% filter(CHR == 3, BP == 101547733) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = -log10(1.3e-10)-0.3+1.2, ymax = -log10(1.3e-10)+1.2+3,
               alpha = alph, fill = "white") +
      annotate('text',
               don %>% filter(CHR == 3, BP == 101547733) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 3, BP == 101547733) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.2,
               label = 'rs17347644',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',
               don %>% filter(CHR == 3, BP == 186696364) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 3, BP == 186696364) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.2,
               label = 'rs16861415',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('rect',
               xmin = don %>% filter(CHR == 3, BP == 195500549) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don %>% filter(CHR == 3, BP == 195500549) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = -log10(1.8e-12)-0.3+2.3, ymax = -log10(1.8e-12)+3+2.3,
               alpha = alph, fill = "white") +
      annotate('text',
               don %>% filter(CHR == 3, BP == 195500549) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 3, BP == 195500549) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+2.3,
               label = 'rs2550250',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',
               don %>% filter(CHR == 10, BP == 111975041) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 10, BP == 111975041) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs1977829',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',
               don %>% filter(CHR == 19, BP == 9007748) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 19, BP == 9007748) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs11673136',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('rect',
               xmin = don %>% filter(CHR == 19, BP == 5831724) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don %>% filter(CHR == 19, BP == 5831724) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = -log10(1.11e-22)-0.2, ymax = -log10(1.11e-22)+2,
               alpha = alph, fill = "white") +
      annotate('text',
               don %>% filter(CHR == 19, BP == 5831724) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 19, BP == 5831724) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric(),
               label = 'rs112313064',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('rect',
               xmin = don %>% filter(CHR == 19, BP == 45418790) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don %>% filter(CHR == 19, BP == 45418790) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = -log10(4.7e-9)-0.2, ymax = -log10(4.7e-9)+2,
               alpha = alph, fill = "white") +
      annotate('text',
               don %>% filter(CHR == 19, BP == 45418790) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 19, BP == 45418790) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric(),
               label = 'rs5117',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('rect',
               xmin = don %>% filter(CHR == 19, BP == 49206462) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don %>% filter(CHR == 19, BP == 49206462) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = -log10(1.4e-17)-0.2+1.25, ymax = -log10(1.4e-17)+3+1.25,
               alpha = alph, fill = "white") +
      annotate('text',
               don %>% filter(CHR == 19, BP == 49206462) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 19, BP == 49206462) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs681343',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0,ylim_MH[ii],5)) +
      ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
  }else if(ii == 4){
    plot1 <- plot1 +
      annotate('text',
               don %>% filter(CHR == 19, BP == 45411941) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 19, BP == 45411941) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1,
               label = 'rs429358',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
    
  }
  
  # Merge both plots -----------------------------------------------------------
  p <- plot1 + annotation_custom(ggplotGrob(plot2),
                                 ymin = y_QQ_min[ii], ymax = y_QQ_max[ii],
                                 xmin = axisdf$center[axisdf$CHR == 12],
                                 xmax = don %>% filter(BPcum == max(BPcum)) %>% reframe(BPcum) %>% as.numeric()
                                 ) 

  ggsave(paste0(dir_results,'Figures/',dir[ii],'.png'), plot = p, width = 80, height = 54, dpi = 300, units = 'mm')
}




