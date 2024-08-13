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

ch  <- c('oneDose', 'twoDose', 'breakthroughSusceptibility','breakthroughSeverity')
ch1 <- c('one_dose', "two_dose", "breakthrough_susceptibility", "breakthrough_severity")
nam <- c('Seroconversion - One dose', 'Seroconversion - Two dose', 
         'Breakthrough susceptibility', 'Breakthrough severity')
nam1 <- c('Seroconversion_One dose', 'Seroconversion_Two doses',
         'Breakthrough_Susceptibility','Breakthrough_Severity')

# Table 1 ======================================================================
outcome_name <- c("immuneResponse","immuneResponse","bt_infection","severity_index")
ukb_tableOne <- loadUKBTableOne() 

tableOne <- list()

for(i in 1:4){
  cohort <- as_tibble(read.table(paste0(dir_results,"Cohorts/", ch1[i], ".txt"), header = TRUE)) |>
    select("eid" = "FID", "Sex", "Age", "outcome" = outcome_name[i])  |>
    inner_join(ukb_tableOne, by = "eid") |>
    mutate(Sex = if_else(Sex == 0, "Female", "Male")) |>
    select(-c("sex", "year_of_birth")) |>
    mutate(outcome = as.factor(outcome)) |>
    rename("Body mass index" = "body_mass_index",
           "Index of multiple deprivation" = "index_of_multiple_deprivation",
           "Ethnic background" = "ethnic_background")
  
  tableOne[[ch[i]]] <- CreateTableOne(data = cohort,
                              vars = colnames(cohort[2:length(colnames(cohort))])) |>
    print(showAllLevels = TRUE) 
}

tab <- tibble(data.frame(tableOne[[ch[1]]])) |>
  rename("Immune response_One dose" = 'Overall') |>
  mutate(Variables = rownames(tableOne[[ch[1]]])) |>
  relocate('Variables') |>
  left_join(
    tibble(data.frame(tableOne[[ch[2]]])) |>
      rename("Immune response_Two doses" = 'Overall') |>
      mutate(Variables = rownames(tableOne[[ch[2]]])) |>
      relocate('Variables'),
    by = c('Variables','level')) |>
  left_join(
    tibble(data.frame(tableOne[[ch[3]]])) |>
      rename("Breakthrough_Susceptibility" = 'Overall') |>
      mutate(Variables = rownames(tableOne[[ch[3]]])) |>
      relocate('Variables'),
    by = c('Variables','level')) %>%
  left_join(
    tibble(data.frame(tableOne[[ch[4]]])) |>
      rename("Breakthrough_Severity" = 'Overall') |>
      mutate(Variables = rownames(tableOne[[ch[4]]])) |>
      relocate('Variables'),
      by = c('Variables','level')) |>
  rename(' ' = 'level') |>
  flextable() |>
  span_header(sep = "_") |>
  align(j = c(3:5), align = 'center', part = 'all')

save_as_docx(tab, path = paste0(dir_results,'Tables/TableOne.docx'))


# Table 1 - VALIDATION =========================================================
outcome_name <- c("immuneResponse","immuneResponse","bt_infection","severity_index")
ukb_tableOne <- loadUKBTableOne() 

tableOne <- list()

for(i in 1:4){
  cohort <- as_tibble(read.table(paste0(dir_results,"Cohorts/", ch1[i], "_validation.txt"), header = TRUE)) |>
    select("eid" = "FID", "Sex", "Age", "outcome" = outcome_name[i])  |>
    inner_join(ukb_tableOne, by = "eid") |>
    mutate(Sex = if_else(Sex == 0, "Female", "Male")) |>
    select(-c("sex", "year_of_birth")) |>
    mutate(outcome = as.factor(outcome)) |>
    rename("Body mass index" = "body_mass_index",
           "Index of multiple deprivation" = "index_of_multiple_deprivation",
           "Ethnic background" = "ethnic_background")
  
  tableOne[[ch[i]]] <- CreateTableOne(data = cohort,
                                      vars = colnames(cohort[2:length(colnames(cohort))])) |>
    print(showAllLevels = TRUE) 
}


tab <- tibble(data.frame(tableOne[[ch[1]]])) |>
  rename("Immune response_One dose" = 'Overall') |>
  mutate(Variables = rownames(tableOne[[ch[1]]])) |>
  relocate('Variables') |>
  left_join(
    tibble(data.frame(tableOne[[ch[2]]])) |>
      rename("Immune response_Two doses" = 'Overall') |>
      mutate(Variables = rownames(tableOne[[ch[2]]])) |>
      relocate('Variables'),
    by = c('Variables','level')) |>
  left_join(
    tibble(data.frame(tableOne[[ch[3]]])) |>
      rename("Breakthrough_Susceptibility" = 'Overall') |>
      mutate(Variables = rownames(tableOne[[ch[3]]])) |>
      relocate('Variables'),
    by = c('Variables','level')) %>%
  left_join(
    tibble(data.frame(tableOne[[ch[4]]])) |>
      rename("Breakthrough_Severity" = 'Overall') |>
      mutate(Variables = rownames(tableOne[[ch[4]]])) |>
      relocate('Variables'),
    by = c('Variables','level')) |>
  rename(' ' = 'level') |>
  flextable() |>
  span_header(sep = "_") |>
  align(j = c(3:5), align = 'center', part = 'all')

save_as_docx(tab, path = paste0(dir_results,'Tables/TableOne_validation.docx'))

# Table SNPs ===================================================================
tableList_SNPs <- list()

for(i in 1:4){
  gwas <- as_tibble(read.table(paste0(dir_results,"GWAS/",ch[i],".txt"), header = TRUE)) |>
    mutate(OR = exp(BETA), PVAL = 10^(-LOG10P))
  
  leadSnps    <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/leadSNPs.txt'))
  annovar     <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/annov.txt'))
  gwascatalog <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/gwascatalog.txt'))
  genomicRL   <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/GenomicRiskLoci.txt'))
  
  order <- getOrder(genomicRL = genomicRL)
  
  t <- order |>
    left_join(leadSnps |>
                select('uniqID',
                       'CHR' = 'chr',
                       'BP'  = 'pos',
                       'SNP' = 'rsID'),
              by = 'SNP') |>
    mutate('Phenotype' = nam[i]) %>%
    left_join(gwas |>
                select('CHR' = 'CHROM',
                       'BP'  = 'GENPOS',
                       'SNP' = 'ID',
                       'EAF' = 'A1FREQ',
                       'EA'  = 'ALLELE1',
                       'OA'  = 'ALLELE0',
                       'OR'  = 'OR',
                       'SE'  = 'SE',
                       'PVAL' = 'PVAL'),
              by = c('CHR','BP','SNP')) |>
    left_join(annovar |>
                group_by(uniqID, 'Function' = annot) |>
                summarise('Gene' = paste0(symbol, collapse=" - ")) |>
                ungroup(),
              by = "uniqID") 
  
  tableList_SNPs[[i]] <- t |>
    select(-uniqID) |>
    mutate('EAF' = round(EAF, digits=2),
           'OR'  = round(OR,  digits=2),
           'SE'  = round(SE,  digits=2),
           'PVAL' = formatC(PVAL, format = "e", digits = 1)) |> 
    select("Phenotype", "SNP", "CHR", "BP",  "EA", "OA","EAF","OR", "SE","PVAL",
           "Function", "Gene/Nearest genes" = "Gene") |>
    flextable() |>
    bg(i = 1, bg = "#EFEFEF", part = "head") %>%
    bold(bold = TRUE, part="header") %>%
    align(align = 'center', part = "all") %>%
    width(j = c('CHR','EAF','EA','OA','OR','SE'), width = 1.2, unit = 'cm') %>%
    width(j = c('SNP','BP','PVAL'), width = 1.7, unit = 'cm') %>%
    width(j = c('Function','Gene/Nearest genes'), width = 2, unit = 'cm') %>%
    hline(part = "body")
}

# Table SNPs (extra info) ------
tableList_SNPs_extra <- list()

for(i in 1:4){
  gwas <- as_tibble(read.table(paste0(dir_results,"GWAS/",ch[i],".txt"), header = TRUE)) |>
    mutate(OR = exp(BETA), PVAL = 10^(-LOG10P))
  
  leadSnps    <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/leadSNPs.txt'))
  annovar     <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/annov.txt'))
  gwascatalog <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/gwascatalog.txt'))
  genomicRL   <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/GenomicRiskLoci.txt'))

  order <- getOrder(genomicRL = genomicRL)
  
  t <- order |>
    left_join(leadSnps |>
                select('uniqID',
                       'CHR' = 'chr',
                       'BP'  = 'pos',
                       'SNP' = 'rsID'),
              by = 'SNP') |>
    #mutate('Trait' = nam[i]) %>%
    left_join(gwas |>
                select('CHR' = 'CHROM',
                       'BP'  = 'GENPOS',
                       'SNP' = 'ID',
                       'EAF' = 'A1FREQ',
                       'EA'  = 'ALLELE1',
                       'OA'  = 'ALLELE0',
                       'OR'  = 'OR',
                       'SE'  = 'SE',
                       'PVAL' = 'PVAL'),
              by = c('CHR','BP','SNP')) |>
    left_join(annovar |>
                group_by(uniqID, 'Function' = annot) |>
                summarise('Gene' = paste0(symbol, collapse=" - "), 
                          'Distance (b)' = paste0(dist, collapse=" - "))|>
                ungroup(),
              by = "uniqID") |>
    left_join(
      gwascatalog |> rename("SNP" = "IndSigSNP") |>
        select('SNP', 'Trait') |>
        distinct() |>
        group_by(SNP) |>
        summarise('Previous reported trait' = paste0(Trait, collapse = ", ")),
      by = 'SNP') |>
    mutate(`Previous reported trait` = if_else(is.na(`Previous reported trait`), '-',`Previous reported trait`))
  
  tableList_SNPs_extra[[i]] <- t |>
    select(-uniqID) |>
    mutate('EAF' = round(EAF, digits=2),
           'OR'  = round(OR,  digits=2),
           'SE'  = round(SE,  digits=2),
           'PVAL' = formatC(PVAL, format = "e", digits = 1)) |> 
    flextable() |>
    bg(i = 1, bg = "#EFEFEF", part = "head") %>%
    bold(bold = TRUE, part="header") %>%
    align(align = 'center', part = "all") %>%
    align(j = 13, align = 'left', part = 'body') %>%
    width(j = c('CHR','EAF','EA','OA','OR','SE'), width = 1.2, unit = 'cm') %>%
    width(j = c('SNP','BP','PVAL'), width = 1.7, unit = 'cm') %>%
    width(j = c('Function','Gene','Distance (b)'), width = 2, unit = 'cm') %>%
    width(j = c('Previous reported trait'), width = 5, unit = 'cm') %>%
    hline(part = "body")
}

# Table Validation -------------------------------------------------------------
tableListValidation <- list()
for(i in 1:4){
  genomicRL <- read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/GenomicRiskLoci.txt'))
  order <- getOrder(genomicRL)
  gwas_validation <- read_delim(paste0(dir_results,'Validation/',ch[i],'_validation.txt'))
  
  t <- order %>% rename('ID' = 'SNP') %>%
    left_join(gwas_validation,
              by = 'ID') %>%  
    mutate(OR = round(exp(BETA), digits = 2)) %>%
    mutate('P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1),
           'Lower'   = exp(BETA-1.96*SE),
           'Upper'   = exp(BETA+1.96*SE),
           "EAF"     = round(A1FREQ,2)) %>%
    select(CHROM, ID, 'Validation_N' = 'N', "EA" = "ALLELE1", 
           "Validation_OR" = 'OR', "Validation_EAF" = "EAF", 
           'Validation_P Value' = 'P Value',
           'Validation_Lower' = 'Lower', 'Validation_Upper' = 'Upper') %>%
    mutate(Phenotype = nam[i]) %>%
    left_join(
      read_delim(paste0(dir_results,'GWAS/',ch[i],'.txt')) %>%
        mutate('Main analysis_N' = N) %>%
        mutate("EAF" = round(A1FREQ,2)) %>%
        mutate('Main analysis_OR' = round(exp(BETA), digits = 2)) %>%
        mutate('Main analysis_P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
        mutate('Main analysis_Lower' = exp(BETA-1.96*SE)) %>%
        mutate('Main analysis_Upper' = exp(BETA+1.96*SE)) %>%
        select(ID, 'Main analysis_N',  "EA" = "ALLELE1",  'Main analysis_OR', "Main analysis_EAF" = "EAF", 
               'Main analysis_P Value', 'Main analysis_Lower', 'Main analysis_Upper'),
      by = c('ID',"EA")
    )
  
  tableListValidation[[i]] <- t
}

gwas <- tableListValidation[[1]] %>%
  add_row(tableListValidation[[2]]) %>%
  add_row(tableListValidation[[3]]) %>%
  add_row(tableListValidation[[4]])
write.table(gwas,paste0(dir_results,'/Validation/Validation.txt'))

tableList_Validation <- gwas %>%
  select(-CHROM,-`Main analysis_Lower`,-`Main analysis_Upper`,-`Validation_Upper`,-`Validation_Lower`) %>%
  select("Phenotype", "SNP" = "ID", "EA", 
         "Main analysis_N", "Main analysis_EAF", "Main analysis_OR", "Main analysis_P Value",
         "Validation_N", "Validation_EAF", "Validation_OR", "Validation_P Value") %>%
  flextable() %>%
  span_header(sep = "_") %>%
  align(align = 'center', part = 'all') %>%
  bg(j = c('Main analysis_N',
           'Main analysis_OR',
           'Validation_N', "Validation_OR"), bg = "#EFEFEF", part = "all")  %>%
  bg(i = 1, j = c(8:11), bg = "white", part = "header") %>%
  vline(j = c(3,7)) %>%
  color(i = ~ as.numeric(`Validation_P Value`) < 0.05, j = c('SNP', "EA",'Main analysis_N', 'Main analysis_EAF','Main analysis_OR','Main analysis_P Value',
                                                             'Validation_N','Validation_EAF','Validation_OR','Validation_P Value'), color = 'red') %>%
  hline(i=c(13,20,38)) %>%
  bold(i = c(1,2), part = "header")


# Table comparator ==============================================================
gwasList  <- list()
genomicRL <- data.frame(Phenotype = NA, CHR = NA, SNP = NA, uniqID = NA, Function = NA, Gene = NA)
annovarList <- list()
for (i in 1:4){
  aux <- getOrder(genomicRL = read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/GenomicRiskLoci.txt'))) %>% 
    left_join(read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/leadSNPs.txt')) %>% select(SNP = rsID, uniqID, CHR = chr), by = c('SNP')) %>%
    left_join(read_delim(paste0(dir_results,'GWAS/FUMA_',ch[i],'/annov.txt'))    %>% 
                group_by(uniqID, 'Function' = annot) %>%
                summarise('Gene' = paste0(symbol, collapse = " - ")) %>% 
                ungroup(),
              by = 'uniqID') %>%
    mutate(Phenotype = nam[i])
  
  genomicRL <- genomicRL %>% add_row(aux)
  
  gwasList[[i]] <- read_delim(paste0(dir_results,'GWAS/',ch[i],'_assoc.regenie.merged.txt')) %>% 
    mutate(OR = exp(BETA))
}

t <- genomicRL %>% 
  filter(!is.na(SNP)) %>% 
  rename('Lead SNP' = 'SNP') %>%
  left_join(gwasList[[1]] %>%
              mutate(OR = round(OR, digits = 2),
                     SE = round(SE, digits = 2),
                     P  = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
              select('Lead SNP' = 'ID',
                     'Immune response_One dose_OR' = OR,
                     'Immune response_One dose_SE' = SE,
                     'Immune response_One dose_P Value' = P),
            by = 'Lead SNP') %>%
  left_join(gwasList[[2]] %>%
              mutate(OR = round(OR, digits = 2),
                     SE = round(SE, digits = 2),
                     P  = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
              select('Lead SNP' = 'ID',
                     'Immune response_Two doses_OR' = OR,
                     'Immune response_Two doses_SE' = SE,
                     'Immune response_Two doses_P Value' = P),
            by = 'Lead SNP') %>%
  left_join(gwasList[[3]] %>%
              mutate(OR = round(OR, digits = 2),
                     SE = round(SE, digits = 2),
                     P  = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
              select('Lead SNP' = 'ID',
                     'Breakthrough_Susceptibility_OR' = OR,
                     'Breakthrough_Susceptibility_SE' = SE,
                     'Breakthrough_Susceptibility_P Value' = P),
            by = 'Lead SNP') %>%
  left_join(gwasList[[4]] %>%
              mutate(OR = round(OR, digits = 2),
                     SE = round(SE, digits = 2),
                     P  = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
              select('Lead SNP' = 'ID',
                     'Breakthrough_Severity_OR' = OR,
                     'Breakthrough_Severity_SE' = SE,
                     'Breakthrough_Severity_P Value' = P),
            by = 'Lead SNP')

tableList_Comparator <- t %>% 
  select(-"uniqID", -"Function") %>%
  flextable() %>%
  span_header(sep = "_") %>%
  bold(bold = TRUE, part="header") %>%
  align(align = 'center', part = "all") %>%
  align(i = 1, align = 'center', part = "header") %>%
  align(i = 2, align = 'center', part = "header") %>%
  align(i = 3, align = 'center', part = "header") %>%
  align(align = "center",part = "all") %>%
  # vline(i = 1, j = c(4:6), part = "header") %>%
  # vline(i = 2, j = c(4:9), part = "header") %>%
  # vline(i = 3, j = seq(4,10,2), part = "header") %>%
  # vline(j = seq(4,10,2), part = "body") %>%
  # bg(j = "Immune response_One dose_OR", bg = "#EFEFEF", part = "all") %>%
  # bg(j = "Immune response_One dose_P Value", bg = "#EFEFEF", part = "all") %>%
  # bg(j = "Immune response_Two dose_SE",    bg = "#EFEFEF", part = "body") %>%
  # bg(j = "Breakthrough_Susceptibility_OR",bg = "#EFEFEF", part = "body") %>%
  # bg(j = "Breakthrough_Susceptibility_P Value",bg = "#EFEFEF", part = "body") %>%
  # bg(j = "Breakthrough_Severity_SE",bg = "#EFEFEF", part = "body") %>%
  # bg(i = 3, j = c(5,7,9,11), bg = "#EFEFEF", part = "header") %>%
  # bg(i = 2, j = c(5,9), bg = "#EFEFEF", part = "header") %>%
  # bg(i = 1, j = c(5), bg = "#EFEFEF", part = "header") %>%
  width(j = 1, width = 1, unit = "in") %>%
  bold(bold = TRUE, part="header") %>%
  width(j = "Immune response_One dose_OR", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_Two doses_OR", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Susceptibility_OR", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Severity_OR", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_One dose_SE", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_Two doses_SE", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Susceptibility_SE", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Severity_SE", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_One dose_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Immune response_Two doses_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Breakthrough_Susceptibility_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Breakthrough_Severity_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Lead SNP", width = 2.7, unit = 'cm') %>%
  width(j = "CHR", width = 2.65, unit = 'cm') %>%
  width(j = "Gene", width = 2.3, unit = 'cm') %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname='Calibri',part = "all")



tList <- list(tableList_SNPs[[1]],
              tableList_SNPs[[2]],
              tableList_SNPs[[3]],
              tableList_SNPs[[4]],
              tableList_SNPs_extra[[1]],
              tableList_SNPs_extra[[2]],
              tableList_SNPs_extra[[3]],
              tableList_SNPs_extra[[4]],
              tableList_Validation, 
              tableList_Comparator)
save_as_docx(values = tList, path = paste0(dir_results,'Tables/AllTables.docx'))




