# ============================================================================ #
#                           BREAKTHROUGH IMPUTED DATA                          #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
# Create cohort ----------------------------------------------------------------
print("load databases")
covid19_result_england <- loadCovidEnglandVaccination() |> as_tibble()
ukb_covariates   <- loadCovariates() |> as_tibble()
vaccination_data <- loadVaccinationData() |> as_tibble()

# Attrition --------------------------------------------------------------------
print("load attrition")
attr(vaccination_data, "cohort_attrition") <- tibble("N" = vaccination_data |> nrow(), 
                                              "Reason"   = "UK Biobank participants",
                                              "Excluded" = 0)

# Cohort -----------------------------------------------------------------------
print("load breakthrough cohort")
breakthrough <- vaccination_data |>
  # Restrict to vaccinated people
  select("eid" = "f.eid",
         "vaccination_status" = "f.32040.0.0",
         "date_first_dose"    = "f.32041.0.0") |>
  filter(!is.na(vaccination_status)) |>
  recordAttrition("Restrict to vaccinated people") |>
  left_join(
    covid19_result_england |>
      filter(result == 1) |>
      group_by(eid) |>
      filter(specdate == min(specdate)) |>
      select("eid", "specdate", "result", "origin") |>
      distinct() |> 
      group_by(eid) |>
      mutate(origin = max(origin)) |> 
      ungroup() |>
      distinct(),
    by = "eid"
  ) |>
  mutate(date_first_dose = date_first_dose + 15) |>
  filter(specdate >= date_first_dose | is.na(specdate)) |>
  recordAttrition("Restrict to participants with no infections prior than two weeks after the first vaccine") |>
  inner_join(
    ukb_covariates,
    by = "eid"
  ) |> 
  filter(Caucasian == 1) |>
  recordAttrition(reason = "participants with European ancestry") |>
  filter(Sex == genetic_sex) |>
  recordAttrition(reason = "participants with the same sex and genetic sex registered") |>
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  recordAttrition(reason = "participants with no sex chromosome aneuploidy") |>
  # filter(kinship_to_other_participants == 0) |>
  # recordAttrition(reason = "no kinship to other participants") |>
  mutate(Age = year(date_first_dose) - Year_of_birth) |>
  mutate(result = if_else(is.na(result), 0, result),
         origin = if_else(is.na(origin), 0, origin)) |>
  mutate("IID" = eid, "FID" = eid) |>
  select(FID, IID, bt_infection = result, severity_index = origin, Sex, Age, Genetic_batch, starts_with("PC")) 

attr(breakthrough,"cohort_attrition") <- attr(breakthrough,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough |> filter(bt_infection == 1) |> tally() |> pull(),
                           breakthrough |> filter(bt_infection == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))


# Breakthrough severity --------------------------------------------------------
print("load breakthrough severity")
breakthrough_severity <- breakthrough |>
  filter(bt_infection == 1) |>
  recordAttrition("Restrict to cases (infected)") |>
  select(FID, IID, severity_index, Sex, Age, Genetic_batch, starts_with("PC")) |>
  relocate('FID', 'IID')

attr(breakthrough_severity,"cohort_attrition") <- attr(breakthrough,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough_severity |> filter(severity_index == 1) |> tally() |> pull(),
                           breakthrough_severity |> filter(severity_index == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(breakthrough, paste0(dir_results,'/Cohorts/breakthrough_susceptibility.txt'), row.names = FALSE)
write.table(attr(breakthrough,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_susceptibility_attrition.txt'), row.names = FALSE)

write.table(breakthrough_severity, paste0(dir_results,'/Cohorts/breakthrough_severity.txt'), row.names = FALSE)
write.table(attr(breakthrough_severity,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_severity_attrition.txt'), row.names = FALSE)

# write .phe files
print("load .phe files")
breakthrough_susceptibility <- as.phe(breakthrough, "FID", "IID")
breakthrough_severity       <- as.phe(breakthrough_severity, "FID", "IID")

write.phe(paste0(dir_results,'/Cohorts/breakthrough_susceptibility.phe'), breakthrough_susceptibility)
write.phe(paste0(dir_results,'/Cohorts/breakthrough_severity.phe'), breakthrough_severity)


# # Validation -------------------------------------------------------------------
# # Breakthrough infection, general covid ----------------------------------------
# print("load breakthrough validation cohort")
print("load databases")
covid19_result_england <- loadCovidEnglandVaccination() |> as_tibble()
ukb_covariates   <- loadCovariates() |> as_tibble()
vaccination_data <- loadVaccinationData() |> as_tibble()

# Attrition --------------------------------------------------------------------
print("load attrition")
attr(vaccination_data, "cohort_attrition") <- tibble("N" = vaccination_data |> nrow(), 
                                                     "Reason"   = "UK Biobank participants",
                                                     "Excluded" = 0)

# Cohort -----------------------------------------------------------------------
print("load breakthrough cohort")
breakthrough <- vaccination_data |>
  # Restrict to vaccinated people
  select("eid" = "f.eid",
         "vaccination_status" = "f.32040.0.0",
         "date_first_dose"    = "f.32041.0.0") |>
  filter(!is.na(vaccination_status)) |>
  recordAttrition("Restrict to vaccinated people") |>
  left_join(
    covid19_result_england |>
      filter(result == 1) |>
      group_by(eid) |>
      filter(specdate == min(specdate)) |>
      select("eid", "specdate", "result", "origin") |>
      distinct() |> 
      group_by(eid) |>
      mutate(origin = max(origin)) |> 
      ungroup() |>
      distinct(),
    by = "eid"
  ) |>
  mutate(date_first_dose = date_first_dose + 15) |>
  filter(specdate >= date_first_dose | is.na(specdate)) |>
  recordAttrition("Restrict to participants with no infections prior than two weeks after the first vaccine") |>
  inner_join(
    ukb_covariates,
    by = "eid"
  ) |> 
  filter(is.na(Caucasian)) |>
  recordAttrition(reason = "participants with no European ancestry") |>
  filter(Sex == genetic_sex) |>
  recordAttrition(reason = "participants with the same sex and genetic sex registered") |>
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  recordAttrition(reason = "participants with no sex chromosome aneuploidy") |>
  # filter(kinship_to_other_participants == 0) |>
  # recordAttrition(reason = "no kinship to other participants") |>
  mutate(Age = year(date_first_dose) - Year_of_birth) |>
  mutate(result = if_else(is.na(result), 0, result),
         origin = if_else(is.na(origin), 0, origin)) |>
  mutate("IID" = eid, "FID" = eid) |>
  select(FID, IID, bt_infection = result, severity_index = origin, Sex, Age, Genetic_batch, starts_with("PC")) 

attr(breakthrough,"cohort_attrition") <- attr(breakthrough,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough |> filter(bt_infection == 1) |> tally() |> pull(),
                           breakthrough |> filter(bt_infection == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))


# Breakthrough severity --------------------------------------------------------
print("load breakthrough severity")
breakthrough_severity <- breakthrough |>
  filter(bt_infection == 1) |>
  recordAttrition("Restrict to cases (infected)") |>
  select(FID, IID, severity_index, Sex, Age, Genetic_batch, starts_with("PC")) |>
  relocate('FID', 'IID')

attr(breakthrough_severity,"cohort_attrition") <- attr(breakthrough,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough_severity |> filter(severity_index == 1) |> tally() |> pull(),
                           breakthrough_severity |> filter(severity_index == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(breakthrough, paste0(dir_results,'/Cohorts/breakthrough_susceptibility_validation.txt'), row.names = FALSE)
write.table(attr(breakthrough,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_susceptibility_validation_attrition.txt'), row.names = FALSE)

write.table(breakthrough_severity, paste0(dir_results,'/Cohorts/breakthrough_severity_validation.txt'), row.names = FALSE)
write.table(attr(breakthrough_severity,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_severity_validation_attrition.txt'), row.names = FALSE)

# write .phe files
print("load .phe files")
breakthrough_susceptibility <- as.phe(breakthrough, "FID", "IID")
breakthrough_severity       <- as.phe(breakthrough_severity, "FID", "IID")

write.phe(paste0(dir_results,'/Cohorts/breakthrough_susceptibility_validation.phe'), breakthrough_susceptibility)
write.phe(paste0(dir_results,'/Cohorts/breakthrough_severity_validation.phe'), breakthrough_severity)
