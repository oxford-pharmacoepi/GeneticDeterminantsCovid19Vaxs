# ==============================================================================
#                              Validation analysis 
#                          Marta Alcalde-Herraiz, 2023
# ==============================================================================

# Extract lead independent variants --------------------------------------------
ch  <- c("oneDose", "twoDose", "breakthroughSusceptibility", "breakthroughSeverity")
ch1 <- c("one_dose", "two_dose", "breakthrough_susceptibility", "breakthrough_severity")

nam <- c('Immune response - One dose',
         'Immune response - Two dose',
         'Breakthrough susceptibility',
         'Breakthrough severity')

outcome_name <- c("immuneResponse","immuneResponse","bt_infection","severity_index")

leadSNPs <- tibble(
  "chr" = as.integer(),
  "rsID" = as.character()
)

for(i in 1:4){
  leadSNPs <- leadSNPs |>
    union_all(
      read.table(paste0(dir_results,"GWAS/FUMA_",ch[i],"/leadSNPs.txt"), 
                 header = TRUE, sep = "\t") |> as_tibble() |>
        select("chr", "rsID")
    )
} 

for(i in unique(leadSNPs$chr)){
  write_delim(leadSNPs |>
                filter(chr == i) |>
                select(rsID),
              paste0(dir_results,"/Validation/validation_", i,".txt"),
              col_names = FALSE,
              progress = TRUE,
              quote = "none")
}


# Study of the cohort ----------------------------------------------------------
ukb_tableone <- loadUKBTableOne()

for(i in 1:4){
  cohort <- as_tibble(read.table(paste0(dir_results,"Cohorts/",ch1[i],"_validation.txt"), header = TRUE)) 
  
  cohort <- cohort |>
    rename("eid" = "FID") |>
    select(-"IID") |>
    inner_join(
      ukb_tableone,
      by = "eid"
    ) |>
    select("eid", outcome = outcome_name[i], "Sex", "Age", "body_mass_index", 
           "index_of_multiple_deprivation", "ethnic_background")
  
  write_delim(cohort,
              paste0(dir_results,"/Validation/tableOne_",ch[i],".txt"),
              col_names = FALSE,
              progress = TRUE,
              quote = "none")
}


