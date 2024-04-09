# Path to files
pathUKB <- paste0(pathData,"UKBiobank/")

# Load data
ukb_data <- read.table(paste0(pathUKB,"alzheimer_data.tab"), header=TRUE, sep="\t") |> as_tibble()
ukb_data$f.42018.0.0 <- as.Date(ukb_data$f.42018.0.0)
ukb_data$f.42020.0.0 <- as.Date(ukb_data$f.42020.0.0)
ukb_data$f.42022.0.0 <- as.Date(ukb_data$f.42022.0.0)
ukb_data$f.42024.0.0 <- as.Date(ukb_data$f.42024.0.0)

# Rename
ukb_data <- ukb_data |>
  select("eid" = "f.eid",
         "sex" = "f.31.0.0",
         "year_of_birth" = "f.34.0.0",
         "diastolic_blood_pressure" = "f.4079.0.0",
         "date_of_dementia_report"  = "f.42018.0.0",
         "date_of_alzheimer_report" = "f.42020.0.0",
         "date_of_vascular_dementia_report" = "f.42022.0.0",
         "date_of_alzheimer_report" = "f.42020.0.0",
         "date_of_frontotemporal_dementia_report" = "f.42024.0.0"
  ) 

cohort <- ukb_data |> 
  filter(!is.na(diastolic_blood_pressure)) |>
  mutate(ad_status = if_else(is.na(date_of_alzheimer_report),0,1))

