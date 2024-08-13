# ============================================================================ #
#                           CREATEFIGURE_VALIDATION                            #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #

ch  <- c('oneDose', 
         'twoDose',
         'breakthroughSusceptibility',
         'breakthroughSeverity')
nam <- c('Seroconversion - One dose',
         'Seroconversion - Two dose',
         'Breakthrough susceptibility',
         'Breakthrough severity')

getOrder <- function(genomicRL){
  order <- genomicRL %>% 
    select(rsID,LeadSNPs,p) %>%
    arrange(p) %>%
    mutate(num = str_count(LeadSNPs,";")) %>%
    separate_rows(LeadSNPs, sep = ";") %>%
    mutate(rsID = if_else(num > 0 & rsID == LeadSNPs, 'NA',rsID)) %>%
    filter(rsID != 'NA') %>%
    group_by(rsID) %>%
    mutate(LeadSNPs = paste0(LeadSNPs, collapse = ";")) %>% 
    distinct() %>%
    mutate(a = 1) %>%
    group_by(a) %>%
    summarise(SNP = paste0(rsID, ";", LeadSNPs, collapse = ";")) %>%
    select(-a) %>%
    separate_rows(SNP, sep = ";") %>% 
    distinct()
  return(order)
}


for(i in 1:4){
  
  leadSnps    <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/leadSNPs.txt'))
  annovar     <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/annov.txt'))
  genomicRL   <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/GenomicRiskLoci.txt'))
  
  genomicRL <- getOrder(genomicRL = genomicRL) |>
    inner_join(genomicRL |> select(SNP = rsID), by = "SNP") %>%
    left_join(leadSnps |>
                select('uniqID',
                       'CHR' = 'chr',
                       'BP'  = 'pos',
                       'SNP' = 'rsID'),
              by = 'SNP') |>
    mutate('Phenotype' = nam[i]) %>%
    left_join(annovar |>
                group_by(uniqID) |>
                summarise('symbol' = paste0(symbol, collapse=" - ")) |>
                ungroup(),
              by = "uniqID")  |>
    select(rsID = SNP, symbol)
  
  
  if(i == 1){
    order <- genomicRL
  }else{
    order <- order %>%
      rbind(genomicRL)
    }
}

order <- data.frame('SNP' = rev(order$rsID), "symbol" = rev(order$symbol))
gwas <- order %>% rename(ID = SNP) %>%
  left_join(
    read.table(paste0(dir_results,'Validation/Validation.txt')),
    by = 'ID')

dodge <- 0.5
w <- 0.3
ww <- 0.6
shape_main <- 15
shape_val  <- 16
c_val <- c('#72ACB9','#D57100','#EFBA00')

plots <- list()
gwas <- gwas %>% 
  mutate(Vali = if_else(sign(Validation_OR-1) == sign(Main.analysis_OR-1),1,0)) %>%
  mutate(Vali = if_else(sign(Validation_OR-1) == sign(Main.analysis_OR-1) & Validation_P.Value <= 0.05,2,Vali))

for(i in c(1:length(nam))){
  tab <- gwas %>% filter(Phenotype == nam[i]) %>%
    select(c('CHROM','ID','Phenotype','Vali', "symbol") | contains('Validation')) %>%
    mutate(Type = 'Validation') %>%
    rename('N' = 'Validation_N',
           'OR' = 'Validation_OR',
           'P.Value' = 'Validation_P.Value',
           'Lower' = 'Validation_Lower',
           'Upper' = 'Validation_Upper')
  tab <- tab %>%
    mutate(order = 1:nrow(tab)) %>%
    full_join(gwas %>% filter(Phenotype == nam[i]) %>%
                select(c('CHROM','ID','Phenotype','Vali',"symbol") | contains('Main.analysis')) %>%
                mutate(Type = 'Main') %>%
                rename('N' = 'Main.analysis_N',
                       'OR' = 'Main.analysis_OR',
                       'P.Value' = 'Main.analysis_P.Value',
                       'Lower' = 'Main.analysis_Lower',
                       'Upper' = 'Main.analysis_Upper') %>%
                mutate(order = 1:nrow(tab))) %>%
    mutate(Type = factor(Type, levels = c('Validation','Main'))) %>%
    mutate(Vali = factor(Vali, levels = c('0','1','2'))) %>%
    ungroup() %>%
    select("Phenotype", "ID", "symbol", "order", "Type", "Lower", "Upper", "Vali", "OR")


  plots[[i]] <- ggplot(tab, aes(x = order,y = OR)) +
    geom_errorbar(
      aes(ymin = Lower, ymax = Upper, color = Vali, shape = Type),
      position = position_dodge(width = dodge), width = w, cex = ww
    ) +
    geom_point(aes(y = OR, color = Vali, shape = Type, size = Type),
               position = position_dodge(width = dodge)) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
    scale_shape_manual(values = c(shape_main,shape_val)) +
    scale_color_manual(values = c("2" = '#72ACB9',"1"='#D57100',"0"='#EFBA00')) +
    scale_size_manual(values = c(2.5,2.5)) +
    theme_gray()
  
  if(i == 4){
  plots[[i]] <- plots[[i]] +
    facet_wrap(~Phenotype, strip.position = 'top') +
    theme(strip.text.x = element_text(size = 11, face = "bold"),
          axis.title.x = element_text(),
          axis.text.x  = element_text(),
          axis.text.y  = element_text(size = 10),
          axis.title.y = element_blank(),
          legend.position = "none",
          text = element_text(family = 'Calibri')) +
    scale_x_continuous(
      breaks = tab$order,
      labels = tab$ID,
      expand = c(0,0),
      limits = c(0.5,nrow(tab)/2+.5),
      sec.axis = sec_axis(~.,
                          breaks = tab$order,
                          labels = tab$symbol)
    ) +
    scale_y_log10(breaks = c(0.5,0.75,1,1.25,1.5, 1.75), limits = c(0.45, 1.85), expand = c(0,0)) +
    labs(y = 'Odds Ratio')+
    coord_flip()
  }else{
    plots[[i]] <- plots[[i]] +
      facet_wrap(~Phenotype, strip.position = 'top') +
      theme(strip.text.x = element_text(size = 11, face = "bold"),
            axis.title.x = element_blank(),
            axis.text.x  = element_text(),
            axis.text.y  = element_text(size = 10),
            axis.title.y = element_blank(),
            legend.position = "none",
            text = element_text(family = 'Calibri')) +
      scale_x_continuous(
        breaks = tab$order,
        labels = tab$ID,
        expand = c(0,0),
        limits = c(0.5,nrow(tab)/2+.5),
        sec.axis = sec_axis(~.,
          breaks = tab$order,
          labels = tab$symbol)
      ) +
      scale_y_log10(breaks = c(0.5,0.75,1,1.25,1.75), limits = c(0.45, 1.85), expand = c(0,0)) +
      coord_flip()
  }

}


leg1 <- data.frame(
  x = c(1,2,3),
  y = c(2,3,4),
  type =  factor(c('Fully validated','Partially validated','Not validated'), levels = c('Fully validated','Partially validated','Not validated'))
)
leg2 <- data.frame(
  or = c(0.1,0.1,0.1), lower = c(0,0,0), upper = c(0.2, 0.2,.2),
  y = factor(c(3,2,1)),
  type = factor(c('Main','Validation','Main'), levels = c('Main','Validation')),
  vali = factor(c('Fully validated','Partially validated','Not validated'), levels = c('Fully validated','Partially validated','Not validated'))
)


legq1 <- ggplot(leg1, aes(x,y, fill = type)) +
  geom_bar(stat="identity") +
  theme_bw()+
  scale_fill_manual(values = c('Fully validated' = '#72ACB9',
                               'Partially validated' = '#D57100',
                               'Not validated' = '#EFBA00'))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        text = element_text(family = 'Calibri')) +
  scale_y_continuous(expand = c(0,0)) 

legq2 <- ggplot(leg2, aes(or, y)) +
  geom_errorbar(aes(xmin = lower, xmax = upper, shape = type),
                width = w, cex = ww) +
  geom_point(aes(shape = type, size = type)) +
  theme_bw()+
  scale_shape_manual(values = c('Main' = 16,'Validation' = 15)) +
  scale_size_manual(values = c('Main' = 5, 'Validation' = 5)) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        text = element_text(family = 'Calibri'))


aux <- ggplot() + theme_void() + ylim(-9,45) + xlim(-1.2,1) +
  annotation_custom(ggplotGrob(legq1), # Legend
                    ymin = 40, ymax = 48,
                    xmin = -.7, xmax = 1) +
  annotation_custom(ggplotGrob(legq2),
                    ymin = 40.5, ymax = 46.75,
                    xmin = -.7, xmax = 1) +
  annotation_custom(ggplotGrob(plots[[1]]), # One dose
                    ymin = 40, ymax = 45.52,
                    xmin = -.735,xmax = 1.0775) +
  annotation_custom(ggplotGrob(plots[[2]]), # Two doses
                    ymin = 35.2, ymax = 40.7,
                    xmin = -.735,xmax = 0.93) +
  annotation_custom(ggplotGrob(plots[[3]]), # Breakthrough susceptibility
                    ymin = 15.1, ymax = 35.8,
                    xmin = -.75, xmax = 1.085) +
  annotation_custom(ggplotGrob(plots[[4]]), # Breakthrough severity
                    ymin = 12.2, ymax = 15.7,
                    xmin = -.7075, xmax = 0.865)

ggsave(paste0(dir_results,'Figures/Validation.png'), plot = aux, width = 30, height = 55, dpi = 600, units = 'cm')

