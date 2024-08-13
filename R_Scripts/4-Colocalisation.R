# ==============================================================================
#                            Colocalisation analysis 
#                          Marta Alcalde-Herraiz, 2023
# ==============================================================================
ch  <- c('oneDose', 
         'twoDose',
         'breakthroughSusceptibility',
         'breakthroughSeverity')
nam <- c('Immune response - One dose',
         'Immune response - Two dose',
         'Breakthrough susceptibility',
         'Breakthrough severity')

genRisk <- tibble("Phenotype"    = as.character(),
                  "GenomicLocus" = as.integer(),
                  "uniqID"       = as.character(),
                  "rsID"         = as.character(),
                  "chr"          = as.integer(),
                  "pos"          = as.integer(),
                  "p"            = as.double(),
                  "log10p"       = as.double())

gwasList <- list()
for(i in 1:4){
  genRisk <- genRisk |>
    union_all(
      read.table(paste0(dir_results,"GWAS/FUMA_",ch[i],"/GenomicRiskLoci.txt"), 
                 header = TRUE, sep = "\t") |> as_tibble() |>
        select("GenomicLocus", "uniqID", "rsID", "chr", "pos", "p") |>
        arrange(p) |>
        mutate(log10p = -log10(p)) |>
        mutate(Phenotype = ch[i])
    ) |>
    distinct()
  
  gwasList[[i]] <- read.delim(paste0(dir_results,'GWAS/',ch[i],'_assoc.regenie.merged.txt'),sep = "\t") |> as_tibble()
}

w <- 250e3
n <- matrix(0,nrow(genRisk),4)
c <- matrix(0,nrow(genRisk),4)

for(i in 1:nrow(genRisk)){
  d1 <- getColocFormat(gwasList[[1]], genRisk, i, w)
  d2 <- getColocFormat(gwasList[[2]], genRisk, i, w)
  d3 <- getColocFormat(gwasList[[3]], genRisk, i, w)
  d4 <- getColocFormat(gwasList[[4]], genRisk, i, w)
  
  if(genRisk[i,]$Phenotype == "oneDose"){
    d <- d1
  }else if(genRisk[i,]$Phenotype == "twoDose"){
    d <- d2
  }else if(genRisk[i,]$Phenotype == "breakthroughSusceptibility"){
    d <- d3
  }else if(genRisk[i,]$Phenotype == "breakthroughSeverity"){
    d <- d4
  }

  c1 <- coloc.abf(dataset1 = d, dataset2 = d1, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c2 <- coloc.abf(dataset1 = d, dataset2 = d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c3 <- coloc.abf(dataset1 = d, dataset2 = d3, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c4 <- coloc.abf(dataset1 = d, dataset2 = d4, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

  c[i,1] <- round(c1$summary[6]*100, digits = 2)
  c[i,2] <- round(c2$summary[6]*100, digits = 2)
  c[i,3] <- round(c3$summary[6]*100, digits = 2)
  c[i,4] <- round(c4$summary[6]*100, digits = 2)

  n[i,1] <- c1$summary[1]
  n[i,2] <- c2$summary[1]
  n[i,3] <- c3$summary[1]
  n[i,4] <- c4$summary[1]
}

coloc <- genRisk |>
  select("Phenotype", "SNP" = "rsID", "CHR" = "chr", "POS" = "pos") |>
  mutate(Phenotype = case_when(
    Phenotype == "oneDose" ~ "Immune response - One dose",
    Phenotype == "twoDose" ~ "Immune response - Two dose",
    Phenotype == "breakthroughSusceptibility" ~ "Breakthrough susceptibility",
    Phenotype == "breakthroughSeverity" ~ "Breakthrough severity"
  )) |>
  mutate("Immune response_One dose_N" = n[,1],
         'Immune response_One dose_P(H4) (%)' = c[,1],
         'Immune response_Two dose_N' = n[,2],
         'Immune response_Two dose_P(H4) (%)' = c[,2],
         'Breakthrough_Susceptibility_N' = n[,3],
         'Breakthrough_Susceptibility_P(H4) (%)' = c[,3],
         'Breakthrough_Severity_N' = n[,4],
         'Breakthrough_Severity_P(H4) (%)' = c[,4]) |>
  flextable() |>
  span_header(sep = "_") |>
  bold(bold = TRUE, part="header") |>
  align(align = 'center', part = "all") |>
  align(i = 1, align = 'center', part = "header") |>
  align(i = 2, align = 'center', part = "header") |>
  align(i = 3, align = 'center', part = "header") |>
  align(align = "center",part = "all") |>
  width(j = 'Phenotype', width = 3.45, unit = 'cm') |>
  width(j = 'SNP', width = 3.25, unit = 'cm') |>
  width(j = 'CHR', width = 1.6, unit = 'cm') |>
  width(j = 'POS', width = 2.35, unit = 'cm') |>
  width(j = 'Immune response_One dose_N', width = 1.45, unit = 'cm') |>
  width(j = 'Immune response_Two dose_N', width = 1.45, unit = 'cm') |>
  width(j = 'Breakthrough_Susceptibility_N', width = 1.45, unit = 'cm') |>
  width(j = 'Breakthrough_Severity_N', width = 1.45, unit = 'cm') |>
  width(j = 'Immune response_One dose_P(H4) (%)', width = 2, unit = 'cm') |>
  width(j = 'Immune response_Two dose_P(H4) (%)', width = 2, unit = 'cm') |>
  width(j = 'Breakthrough_Susceptibility_P(H4) (%)', width = 2, unit = 'cm') |>
  width(j = 'Breakthrough_Severity_P(H4) (%)', width = 2, unit = 'cm') |>
  bg(j = 'Immune response_One dose_N', bg = "#EFEFEF", part = "all") |> 
  bg(j = 'Immune response_Two dose_N', bg = "#EFEFEF", part = "body") |>
  bg(i = 3, j = c(7,9,11), bg = "#EFEFEF", part = "all") |>
  bg(j = 'Breakthrough_Susceptibility_N', bg = "#EFEFEF", part = "body") |>
  bg(i = 2, j = 9, bg = "#EFEFEF", part = "header")  |>
  bg(j = 'Breakthrough_Severity_N', bg = "#EFEFEF", part = "body") |>
  fontsize(size = 11, part = "all") |>
  font(fontname='Calibri', part = "all") |>
  hline(i = c(2,4,14), part = "body") |>
  vline(j = c(4,6,8,10))
save_as_docx(coloc, path = paste0(dir_results,'Tables/Colocalisation.docx'))

  
# Study of those variants that colocalise
# genRisk <- genRisk %>% filter(rsID == "rs681343")
# gwas1 <- gwasList[[1]] %>%
#   filter(GENPOS >= genRisk$pos[1]-w-1 & GENPOS <= genRisk$pos[1]+w+1 & CHROM == genRisk$chr[1]) %>%
#   mutate(f = if_else(ID == "rs681343",1,0)) %>%
#   arrange(f)
# gwas3 <- gwasList[[3]] %>%
#   filter(GENPOS >= genRisk$pos[1]-w-1 & GENPOS <= genRisk$pos[1]+w+1 & CHROM == genRisk$chr[1]) %>%
#   mutate(f = if_else(ID == "rs681343",1,0)) %>%
#   arrange(f)
# 
# library(extrafont)
# font_import()
# windowsFonts("Calibri" = windowsFont("Calibri"))
# 
# 
# gplot1 <- ggplot(gwas1, aes(x = GENPOS, y = LOG10P)) +
#   # Show all points
#   geom_point(size = 1.25, aes(colour = as.factor(f)), shape = 16) +
#   scale_color_manual(values = c("#4393C3",  "red")) +
#   # custom X axis:
#   scale_y_continuous(breaks = seq(0,6,1), limits = c(0,6.2), expand = c(0,0)) +     # remove space between plot area and x axis
#   # # Custom the theme:
#   theme_bw() +
#   theme(
#     legend.position="none",
#     # panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.line.x = element_line(color = "black", linewidth = 0.5),
#     axis.line.y = element_line(color = "black", linewidth = 0.5),
#     axis.text.x = element_text(margin = margin(t = 1), size = 11),
#     axis.text.y  = element_text(size = 11),
#     axis.title = element_text(size = 11),
#     plot.title = element_text(face = "bold"),
#     text = element_text(family = "Calibri")
#   ) +
#   #Plot a red horizontal line at 5e-8
#   geom_hline(yintercept = -log10(5e-8), linewidth = 0.3) +
#   labs(x = 'Chromosome 19 [BP]', y = expression(-log[10](P))) +
#   ggtitle('(A) Immune response - one dose')
# 
# gplot3 <- ggplot(gwas3, aes(x = GENPOS, y = LOG10P)) +
#   # Show all points
#   geom_point(size = 1.25, aes(colour = as.factor(f)), shape = 16) +
#   scale_color_manual(values = c("#4393C3",  "red")) +
#   # custom X axis:
#   scale_y_continuous(breaks = seq(0,50,2), limits = c(0,50.2), expand = c(0,0)) +     # remove space between plot area and x axis
#   # # Custom the theme:
#   theme_bw() +
#   theme(
#     legend.position="none",
#     # panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.line.x = element_line(color = "black", linewidth = 0.5),
#     axis.line.y = element_line(color = "black", linewidth = 0.5),
#     axis.text.x = element_text(margin = margin(t = 1), size = 11),
#     axis.text.y  = element_text(size = 11),
#     axis.title = element_text(size = 11),
#     plot.title = element_text(face = "bold"),
#     text = element_text(family = "Calibri")
#   ) +
#   #Plot a red horizontal line at 5e-8
#   labs(x = 'Chromosome 19 [BP]', y = expression(-log[10](P))) +
#   ggtitle('(B) Breakthrough susceptibility')
# 
# ggsave('mh_1.png', plot = gplot1, path = paste0(dir_results,'Figures/mh_1.png'), height = 6.7, width = 22.4, units = 'cm',dpi = 300)
# ggsave('mh_3.png', plot = gplot3, path = paste0(dir_results,'Figures/mh_3.png'), height = 6.7, width = 22.4, units = 'cm',dpi = 300)

