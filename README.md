# GeneticDeterminantsCovid19Vaxs

**Background**: The impact of genetics on the heterogeneous immune response and risks of breakthrough infections remains unclear for the different COVID-19 vaccines. Understanding the role of genetic variants on seroconversion and subsequent breakthrough outcomes is key to reveal mechanisms underlying vaccine effectiveness. 

**Objective**: To perform four genome wide association studies among vaccinated participants for COVID-19 vaccine seroconversion and breakthrough outcomes. 

For more information see the original publication [here](https://www.nature.com/articles/s41467-024-52890-6).

# Running the analysis
Previously, you must ensure that you have a directory "${folder directory}" with the following files:
 - breakthrough.rds [RDS file]: File with breakthrough information. Columns required include
   
        -> eid: individual identification
        -> event_dt: vaccination date
        -> case_control_index: if it was infected within 84 days after vaccination (Cases) // not infected within 84 days after vaccination (Control)
        -> breakthrough_infection_type: Omicron / No information (NA) / Pre_Omicron
        -> origin: if it was hospitalised (1) / not hospitalised (0) / no information (NA)
   
 - immune_response.tab [TAB FILE]: File with immune response information. Columns required include

        -> eid
        -> Antibody_test_result: Field 27981.0.0 from UK Biobank
        -> Date_antibody_test_performed: Field 27982.0.0 from UK Biobank
        -> Received_first_COVID19_vaccination: Field 27983.0.0 from UK Biobank
        -> Date_first_COVID19_vaccination: Field 27984.0.0 from UK Biobank
        -> Received_second_COVID19_vaccination: Field 27985.0.0 from UK Biobank
        -> Date_second_COVID19_vaccination: Field 27986.0.0 from UK Biobank
        -> Covid_antibody_test_results: Field 27990.0.0 from UK Biobank
   
 - ukb669864_covariates.csv [CSV file]: File with covariates information. Columns required include

        -> Sex: Field 31 from UK Biobank
        -> BMI: Field 21001.0.0 from UK Biobank
        -> PC1,...,PC10: Field 22009 from UK Biobank
        -> Year_of_birth: Field 34 from UK Biobank
        -> Genetic batch: Field 22000 from UK Biobank
        -> Caucasian: Field 22006 from UK Biobank
        -> sex_chromosome_aneuploidy: Field 22019 from UK Biobank
        -> kinship_to_other_participants: Field 22021 from UK Biobank
        -> genetic_sex: Field 22001 from UK Biobank

Also, the following variables must be specified:

    -> dir_data <- "...": The path to the data directory.
    -> dir_ukb  <- "...": The path where ukb files are

The ukb directory must contain the file "ukb65397.tab", with the following columns:

     -> f.21001.0.0: Body mass index
     -> f.26410.0.0: Index of multiple deprivation (England)
     -> f.26427.0.0: Index of multiple deprivation (Scotland)
     -> f.26426.0.0: Index of multiple deprivation (Wales)
     -> f.26411.0.0: Incomde score (England)
     -> f.26428.0.0: Income score (Scotland)
     -> f.26418.0.0: Income score (Wales)
     -> f.26412.0.0: Employment score (England)
     -> f.26419.0.0: Employment score (Wales)
     -> f.26413.0.0: Health score (England)
     -> f.26430.0.0: Health score (Scotland)
     -> f.26420.0.0: Health score (Wales)
     -> f.26414.0.0: Education score (England)
     -> f.26421.0.0: Education score (Wales)
     -> f.26431.0.0: Education score (Scotland)
     -> f.26415.0.0: Housing score (England)
     -> f.26432.0.0: Housing score (Scotland)
     -> f.26423.0.0: Housing score (Wales)
     -> f.26416.0.0: Crime score (England)
     -> f.26434.0.0: Crime score (Scotland)

 
Results will be saved into the dir_data directory, within a "Results" folder. For more instructions when runing the analysis, see CodeToRun.R
    
