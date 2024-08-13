

for(phenotype in  c("oneDose","twoDose","breakthroughSusceptibility","breakthroughSeverity")){
  
  x <- read.table(paste0(dir_results,"/GWAS/", phenotype, "_assoc.regenie.merged.txt"),
                  header = TRUE)
  
  x <- x |> dplyr::as_tibble() |> dplyr::mutate(PVAL = 10^(-LOG10P))
  
  readr::write_delim(x,paste0(dir_results,"/GWAS/", phenotype,".txt"),
                     col_names = TRUE,
                     progress = TRUE,
                     quote = "none")
  
}