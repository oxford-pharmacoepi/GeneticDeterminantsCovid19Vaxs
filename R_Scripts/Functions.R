recordAttrition <- function(table, reason){
  n <- attr(table, "cohort_attrition") |> 
    filter(row_number() == max(row_number())) |> 
    pull(N)
  
  attr(table, "cohort_attrition") <- attr(table, "cohort_attrition") |> 
    union_all(
      tibble(
        "N"      = table |> nrow(),
        "Reason" = reason,
        "Excluded" = n -(table |> nrow())
      )
    ) |> 
    distinct()
  return(table)
}

loadVaccinationData <- function(){
  bd <- read.table(paste0(dir_data,"UKBiobank/vaccination_data.tab"), header=TRUE, sep="\t")
  bd$f.32041.0.0 <- as.Date(bd$f.32041.0.0)
  bd$f.32041.1.0 <- as.Date(bd$f.32041.1.0)
  bd$f.32041.2.0 <- as.Date(bd$f.32041.2.0)
  bd$f.32041.3.0 <- as.Date(bd$f.32041.3.0)
  bd$f.32041.4.0 <- as.Date(bd$f.32041.4.0)
  bd$f.32041.5.0 <- as.Date(bd$f.32041.5.0)
  bd$f.32041.6.0 <- as.Date(bd$f.32041.6.0)
  bd$f.32041.7.0 <- as.Date(bd$f.32041.7.0)
  bd$f.32041.8.0 <- as.Date(bd$f.32041.8.0)
  bd$f.32041.9.0 <- as.Date(bd$f.32041.9.0)
  
  return(bd)
}

loadImmuneResponseData <- function(){
  bd <- read.table(paste0(dir_data,"UKBiobank/immuneResponse_data.tab"), header=TRUE, sep="\t")
  bd$f.27982.0.0 <- as.Date(bd$f.27982.0.0)
  bd$f.27982.1.0 <- as.Date(bd$f.27982.1.0)
  bd$f.27984.0.0 <- as.Date(bd$f.27984.0.0)
  bd$f.27984.1.0 <- as.Date(bd$f.27984.1.0)
  bd$f.27986.0.0 <- as.Date(bd$f.27986.0.0)
  bd$f.27986.1.0 <- as.Date(bd$f.27986.1.0)
  
  return(bd)
}

loadCovariates <- function(){
  ukb_covar <- as_tibble(read.delim(paste0(dir_data,'UKBiobank/ukb_covariates.tab'), sep = "\t")) |>
    select("eid" = "f.eid",
           "Sex" = "f.31.0.0",
           "PC1" = "f.22009.0.1", "PC2" =	"f.22009.0.2","PC3" = "f.22009.0.3", 
           "PC4" = "f.22009.0.4", "PC5" =	"f.22009.0.5","PC6" = "f.22009.0.6", 
           "PC7" = "f.22009.0.7", "PC8" =	"f.22009.0.8","PC9" = "f.22009.0.9", 
           "PC10" = "f.22009.0.10", 
           "Year_of_birth" = "f.34.0.0",	  "Genetic_batch" = "f.22000.0.0",
           "Caucasian"	   = "f.22006.0.0", "sex_chromosome_aneuploidy" = "f.22019.0.0",	
           "kinship_to_other_participants" = "f.22021.0.0",	"genetic_sex" = "f.22001.0.0")
  
  return(ukb_covar)
}

loadCovidEnglandVaccination <- function(){
  read.table(paste0(dir_data,"UKBiobank/covid19_result_england.txt"), header = TRUE, sep = "\t") |>
    as_tibble() |>
    mutate(specdate = lubridate::dmy(specdate)) 
  
}

getColocFormat <- function(gwas,genRisk,i,w){
  gwas1 <- gwas %>%
    filter(CHROM == genRisk$chr[i]) %>%
    filter(GENPOS >= genRisk$pos[i]-w-1 & GENPOS <= genRisk$pos[i]+w+1)
  
  d <- list(
    beta = gwas1$BETA,
    varbeta = gwas1$SE^2,
    type = 'cc',
    snp  = gwas1$ID,
    position = gwas1$GENPOS
  )
  
  return(d)
}

loadUKBTableOne <- function(){
  
  bd <- as_tibble(read.table(paste0(dir_data,"UKBiobank/tableOne.tab"), header=TRUE, sep="\t")) |>
    rename("eid" = "f.eid",
           "sex" = "f.31.0.0",
           "year_of_birth" = "f.34.0.0",
           "ethnic_background" = "f.21000.0.0",
           "body_mass_index" = "f.21001.0.0") |>
    mutate("index_of_multiple_deprivation" = if_else(is.na(f.26410.0.0), f.26427.0.0, f.26410.0.0)) |>
    mutate("index_of_multiple_deprivation" = if_else(is.na(index_of_multiple_deprivation), f.26426.0.0, index_of_multiple_deprivation)) |>
    inner_join(
      as_tibble(read.table(paste0(dir_data,"UKBiobank/ethnic_background_coding.tsv"), header = TRUE, sep = "\t")) |>
      select("ethnic_background" = "coding", "meaning"),
      by = "ethnic_background"
    ) |>
    select(-c("ethnic_background")) |>
    rename("ethnic_background" = "meaning") |>
    select(-starts_with("f."))
    
  return(bd)
    
}