# ============================================================================ #
#                               IMMUNE RESPONSE                                #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
# Cohort -----------------------------------------------------------------------
print("load databases")

ukb_covid <- as_tibble(loadImmuneResponseData()) |>
  select("eid" = "f.eid",
         "Antibody_test_result"         = "f.27981.0.0",
         "Date_antibody_test_performed" = "f.27982.0.0",
         "Received_first_COVID19_vaccination"  = "f.27983.0.0",
         "Date_first_COVID19_vaccination"      = "f.27984.0.0",
         "Received_second_COVID19_vaccination" = "f.27985.0.0",
         "Date_second_COVID19_vaccination" = "f.27986.0.0",
         "Covid_antibody_test_results"     = "f.27990.0.0") |>
  mutate(Antibody_test_result = if_else(Antibody_test_result == 2,1,Antibody_test_result)) %>%
  mutate(Date_first_COVID19_vaccination  = as.Date(Date_first_COVID19_vaccination,"%Y-%m-%d"),
         Year_first_COVID19_vaccination  = year(Date_first_COVID19_vaccination),
         Date_second_COVID19_vaccination = as.Date(Date_second_COVID19_vaccination,"%Y-%m-%d"),
         Year_second_COVID19_vaccination = year(Date_second_COVID19_vaccination),
         Date_antibody_test_performed    = as.Date(Date_antibody_test_performed,"%Y-%m-%d"),
         Year_antibody_test_performed    = year(Date_antibody_test_performed))

ukb_covar <- as_tibble(loadCovariates())

# Attrition --------------------------------------------------------------------
print("load attrition")
attr(ukb_covid, "cohort_attrition") <- tibble("N" = ukb_covid |> nrow(), 
                                              "Reason"   = "UK Biobank participants",
                                              "Excluded" = 0)
# Main analysis ----------------------------------------------------------------
print("load main analysis")
vaccine_cohort <- ukb_covid |>
  filter(!is.na(Antibody_test_result)) |>
  recordAttrition(reason = "participants in the SARS-CoV-2 coronavirus antibody seroprevalence study") |>
  left_join(ukb_covar, by = "eid") |>
  filter(Caucasian == 1) |>
  recordAttrition(reason = "participants with European ancestry") |>
  filter(Sex == genetic_sex) |>
  recordAttrition(reason = "participants with the same sex and genetic sex registered") |>
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  recordAttrition(reason = "participants with no sex chromosome aneuploidy") |>
  # filter(kinship_to_other_participants == 0) |>
  # recordAttrition(reason = "no kinship to other participants") |>
  filter(Received_first_COVID19_vaccination == 1) %>%
  recordAttrition(reason = "participants that have received at least one dose of a COVID-19 vaccine") |>
  filter(Year_first_COVID19_vaccination >= 2021) |>
  recordAttrition(reason = "participants were vaccinated after 2021") |>
  filter(Covid_antibody_test_results == 0 | is.na(Covid_antibody_test_results)) |>
  recordAttrition(reason = "participants tested negative or did not test for COVID-19 infection") |>
  filter(Year_antibody_test_performed >= 2021) |>
  recordAttrition(reason = "participants COVID-19 test date is valid")

# One dose vaccine cohort ------------------------------------------------------
print("load one dose cohort")
one_vaccine_cohort <- vaccine_cohort |>
  filter(Received_second_COVID19_vaccination == 0) |>
  recordAttrition(reason = "participants have not received a second covid-19 vaccine dose") |> 
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_first_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine >= 8,
         Days_since_vaccine <= 56) |>
  # filter(Days_since_vaccine > 0, Days_since_vaccine <= 84) |>
  recordAttrition(reason = "participants have received the first COVID-19 vaccine between 7 and 56 days before the antibody test") |> 
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) |>
  mutate('FID' = eid, 'IID' = eid) |>
  relocate('FID', 'IID')
write.table(one_vaccine_cohort, paste0(dir_results, "Cohorts/one_dose_crude.txt"), row.names = FALSE)

one_vaccine_cohort <- one_vaccine_cohort |>
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         starts_with("PC"))

attr(one_vaccine_cohort,"cohort_attrition") <- attr(one_vaccine_cohort,"cohort_attrition") |>
  union_all(tibble("N" = c(one_vaccine_cohort |> filter(immuneResponse == 1) |> tally() |> pull(),
                           one_vaccine_cohort |> filter(immuneResponse == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(one_vaccine_cohort, paste0(dir_results,'Cohorts/one_dose.txt'), row.names = FALSE)
write.table(attr(one_vaccine_cohort,"cohort_attrition"), paste0(dir_results,'Cohorts/one_dose_attrition.txt'), row.names = FALSE)

# Two dose vaccine cohort ------------------------------------------------------
print("load two dose cohort")
two_vaccine_cohort <- vaccine_cohort |>
  filter(Received_second_COVID19_vaccination == 1) |>
  recordAttrition(reason = "participants have received a second COVID-19 vaccine dose") |> 
  filter(Year_second_COVID19_vaccination >= 2021) %>%
  recordAttrition(reason = "participants have received after 2021") |> 
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_second_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine <= 56 & Days_since_vaccine >= 8) %>%
  recordAttrition(reason = "participants have received the second vaccine dose between 8 and 56 days before antibody test") |> 
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID', 'IID')
write.table(two_vaccine_cohort, paste0(dir_results, "Cohorts/two_dose_crude.txt"), row.names = FALSE)

two_vaccine_cohort <- two_vaccine_cohort |>
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

attr(two_vaccine_cohort,"cohort_attrition") <- attr(two_vaccine_cohort,"cohort_attrition") |>
  union_all(tibble("N" = c(two_vaccine_cohort |> filter(immuneResponse == 1) |> tally() |> pull(),
                           two_vaccine_cohort |> filter(immuneResponse == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(two_vaccine_cohort, paste0(dir_results,'Cohorts/two_dose.txt'), row.names = FALSE)
write.table(attr(two_vaccine_cohort,"cohort_attrition"), paste0(dir_results,'Cohorts/two_dose_attrition.txt'), row.names = FALSE)

# write .phe files
print("load .phe files")
one_dose_cohort <- as.phe(one_vaccine_cohort, "FID", "IID")
two_dose_cohort <- as.phe(two_vaccine_cohort, "FID", "IID")

write.phe(paste0(dir_results,'/Cohorts/one_dose.phe'), one_dose_cohort)
write.phe(paste0(dir_results,'/Cohorts/two_dose.phe'), two_dose_cohort)

# ============================================================================ #
#                                   Validation                                 #
# ============================================================================ #
print("load validation")
# Main analysis ----------------------------------------------------------------
vaccine_cohort <- ukb_covid |>
  filter(!is.na(Antibody_test_result)) |>
  recordAttrition(reason = "participants in the SARS-CoV-2 coronavirus antibody seroprevalence study") |>
  left_join(ukb_covar, by = "eid") |>
  filter(is.na(Caucasian)) |>
  recordAttrition(reason = "participants with no European ancestry") |>
  filter(Sex == genetic_sex) |>
  recordAttrition(reason = "participants with the same sex and genetic sex registered") |>
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  recordAttrition(reason = "participants with no sex chromosome aneuploidy") |>
  # filter(kinship_to_other_participants == 0) |>
  # recordAttrition(reason = "no kinship to other participants") |>
  filter(Received_first_COVID19_vaccination == 1) %>%
  recordAttrition(reason = "participants that have received at least one dose of a COVID-19 vaccine") |>
  filter(Year_first_COVID19_vaccination >= 2021) |>
  recordAttrition(reason = "participants were vaccinated after 2021") |>
  filter(Covid_antibody_test_results == 0 | is.na(Covid_antibody_test_results)) |>
  recordAttrition(reason = "participants tested negative or did not test for COVID-19 infection") |>
  filter(Year_antibody_test_performed >= 2021) |>
  recordAttrition(reason = "participants COVID-19 test date is valid")

# One dose vaccine cohort ------------------------------------------------------
print("load one dose cohort")
one_vaccine_cohort <- vaccine_cohort |>
  filter(Received_second_COVID19_vaccination == 0) |>
  recordAttrition(reason = "participants have not received a second covid-19 vaccine dose") |> 
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_first_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine >= 8,
         Days_since_vaccine <= 56) |>
  # filter(Days_since_vaccine > 0, Days_since_vaccine <= 84) |>
  recordAttrition(reason = "participants have received the first COVID-19 vaccine between 8 and 56 days before the antibody test") |> 
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) |>
  mutate('FID' = eid, 'IID' = eid) |>
  relocate('FID', 'IID') |>
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         starts_with("PC"))

attr(one_vaccine_cohort,"cohort_attrition") <- attr(one_vaccine_cohort,"cohort_attrition") |>
  union_all(tibble("N" = c(one_vaccine_cohort |> filter(immuneResponse == 1) |> tally() |> pull(),
                           one_vaccine_cohort |> filter(immuneResponse == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(one_vaccine_cohort, paste0(dir_results,'Cohorts/one_dose_validation.txt'), row.names = FALSE)
write.table(attr(one_vaccine_cohort,"cohort_attrition"), paste0(dir_results,'Cohorts/one_dose_validation_attrition.txt'), row.names = FALSE)

# Two dose vaccine cohort ------------------------------------------------------
print("load two dose cohort")
two_vaccine_cohort <- vaccine_cohort |>
  filter(Received_second_COVID19_vaccination == 1) |>
  recordAttrition(reason = "participants have received a second COVID-19 vaccine dose") |> 
  filter(Year_second_COVID19_vaccination >= 2021) %>%
  recordAttrition(reason = "participants have received after 2021") |> 
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_second_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine <= 56 & Days_since_vaccine >= 8) %>%
  recordAttrition(reason = "participants have received the second vaccine dose between 8 and 56 days before antibody test") |> 
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID', 'IID') %>%
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

attr(two_vaccine_cohort,"cohort_attrition") <- attr(two_vaccine_cohort,"cohort_attrition") |>
  union_all(tibble("N" = c(two_vaccine_cohort |> filter(immuneResponse == 1) |> tally() |> pull(),
                           two_vaccine_cohort |> filter(immuneResponse == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(two_vaccine_cohort, paste0(dir_results,'Cohorts/two_dose_validation.txt'), row.names = FALSE)
write.table(attr(two_vaccine_cohort,"cohort_attrition"), paste0(dir_results,'Cohorts/two_dose_validation_attrition.txt'), row.names = FALSE)

# write .phe files
print("load .phe files")
one_dose_cohort <- as.phe(one_vaccine_cohort, "FID", "IID")
two_dose_cohort <- as.phe(two_vaccine_cohort, "FID", "IID")

write.phe(paste0(dir_results,'/Cohorts/one_dose_validation.phe'), one_dose_cohort)
write.phe(paste0(dir_results,'/Cohorts/two_dose_validation.phe'), two_dose_cohort)