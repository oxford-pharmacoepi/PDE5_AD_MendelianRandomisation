exp <- loadGwas(exposure_i, onlyInstruments = FALSE) %>%
  dplyr::filter(chr == 4, pos >= gene_start-window, pos <= gene_end+window)

exp |> select("SNP") |> write.table(paste0(pathResults,"Colocalisation/snp_list.txt"), row.names = FALSE)

# Upload the file to the UK Biobank RAP Platform

# From now on, the analysis was conducted in the UK Biobank R cloud ------------
# Once one has converted the inputed .bgen files to .bed files using plink 
# ("see ConvertBgenToBed.sh"), run in the R Cloud terminal the following line:
# dx download GWAS/Plink_files/ped_files_c4.map
# dx download GWAS/Plink_files/ped_files_c4.ped
# dx download GWAS/Plink_files/plink_files_c4.bim
# dx download GWAS/PDE5_AD/snp_list.txt
# dx download GWAS/PDE5_AD/data_participant.tsv
# dx download GWAS/PDE5_AD/coding22000.tsv
library(readr)
library(dplyr)

# general ----
map_path <- "ped_files_c4.map"
ped_path <- "ped_files_c4.ped"
bim_path <- "plink_files_c4.bim"
snp_path <- "snp_list.txt"

map <- readr::read_delim(map_path,col_names = FALSE) %>%
  dplyr::slice(rep(1:dplyr::n(), each = 2)) %>%
  dplyr::group_by(X2) %>%
  dplyr::mutate(SNP = paste0(X2,"_",dplyr::row_number())) %>%
  dplyr::pull(SNP)

ped <- readr::read_delim(ped_path,col_names = FALSE) %>%
  dplyr::select(-"X2":-"X6")
colnames(ped) <- c("eid", map)

alleles <- readr::read_delim(bim_path,col_names = FALSE) %>%
  dplyr::select("SNP" = "X2", "ALLELE0" = "X5") %>%
  tidyr::pivot_wider(names_from = "SNP", values_from = "ALLELE0") %>%
  dplyr::slice(rep(1:dplyr::n(), each = nrow(ped)))
snps <- colnames(alleles)
alleles <- alleles %>%
  dplyr::rename_with(~paste0(.x,"_0")) %>%
  dplyr::mutate(eid = ped$eid) %>%
  dplyr::inner_join(ped)

for (snp in snps) {
  alleles <- alleles %>%
    dplyr::mutate(
      !!snp := dplyr::if_else(
        .data[[paste0(snp, "_0")]] == .data[[paste0(snp, "_1")]], 1, 0
      ) + dplyr::if_else(
        .data[[paste0(snp, "_0")]] == .data[[paste0(snp, "_2")]], 1, 0
      ) + dplyr::if_else(
        .data[[paste0(snp, "_1")]] == "0", NA, 0
      )
    )
}

# Build cohort
ukb_data <- read.delim("data_participant.tsv") |>
  as_tibble() |>
  select("eid",
         "sex" = "X31.0.0",
         "age_when_assessment" = "X21003.0.0",
         "date_of_dementia_report"  = "X42018.0.0",
         "date_of_alzheimer_report" = "X42020.0.0",
         "date_of_vascular_dementia_report" = "X42022.0.0",
         "date_of_frontotemporal_dementia_report" = "X42024.0.0",
         "PC1" = "X22009.0.1", "PC2" = "X22009.0.2", "PC3" = "X22009.0.3",
         "PC4" = "X22009.0.4", "PC5" = "X22009.0.5", "PC6" = "X22009.0.6",
         "PC7" = "X22009.0.7", "PC8" = "X22009.0.8", "PC9" = "X22009.0.9", 
         "PC10" = "X22009.0.10",
         "genetic_batch" = "X22000.0.0") |>
  left_join(
    read.delim("coding22000.tsv") |>
      select("genetic_batch" = "meaning", "coding"),
    by = "genetic_batch"
  ) |>
  select(-"genetic_batch")|>
  rename("genetic_batch" = "coding")


ukb_data <- ukb_data |> 
  mutate(ad_status = if_else(date_of_alzheimer_report == "",0,1)) |>
  # Exclude people with vascular dementia and ad
  mutate(vascular_status = if_else(date_of_vascular_dementia_report  == "",0,1)) |>
  filter(vascular_status == 0)  |>
  # Exclude people with frontotemporal dementia and ad
  mutate(frontotemporal_status = if_else(date_of_frontotemporal_dementia_report == "",0,1)) |>
  filter(frontotemporal_status == 0) |>
  # Exclude people with other forms of dementia
  mutate(dementia_status = if_else(date_of_dementia_report == "",0,1)) |>
  filter((dementia_status == 0 & ad_status == 0) | (dementia_status == 1 & ad_status == 1))

ukb_data <- ukb_data |> filter(!is.na(PC1))
ukb_data |> group_by(sex, ad_status) |> tally()

# Colocalisation ----
genetics <- alleles %>%
  dplyr::select(eid, all_of(snps), all_of(paste0(snps,"_0"))) %>%
  dplyr::inner_join(ukb_data,by = "eid")

snp_outcome <- tibble(
  "snp" = as.character(),
  "effect_allele.outcome" = as.character(),
  "other_allele,outcome"  = as.character(),
  "samplesize.outcome" = as.numeric(),
  "beta.outcome" = as.numeric(),
  "se.outcome" = as.numeric(),
  "eaf.outcome" = as.numeric(),
  "pval.outcome" = as.numeric()
)

loadSnpOutcomeResults <- function(i, data, regression){
  tibble(
    "snp" = i,
    "effect_allele.outcome" = data |> select(all_of(paste0(i,"_0"))) |> distinct() |> pull(),
    "samplesize.outcome" = data |> tally() |> pull(),
    "beta.outcome" = regression[1],
    "se.outcome" = regression[2],
    "eaf.outcome" = (data |> select("snp") |> pull() |> sum())/(2*(data |> tally() |> pull())),
    "pval.outcome" = regression[4]
  )
}

for(i in snps){
  data <- genetics |>
    select("snp" = all_of(i), all_of(paste0(i,"_0")),"sex","age_when_assessment",
           starts_with("PC"), "genetic_batch","ad_status") |>
    filter(!is.na(snp)) |>
    filter(sex == "Female")
  regression_females <- glm(ad_status ~ snp + age_when_assessment + genetic_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                            data = data,
                            family = "binomial")
  regression_females <- coefficients(summary(regression_females))[row.names(coefficients(summary(regression_females))) == "snp",]
  
  snp_outcome <- snp_outcome |>
    rbind(loadSnpOutcomeResults(i, data = data, regression = regression_females))
}

readr::write_delim(snp_outcome,
                   paste0("Colocalisation_outcome.txt"), delim = "\t")



snps_exposure <- c("rs17355550","rs80223330","rs12646525","rs66887589","rs10050092")
snp_outcome_males <- tibble(
  "snp" = as.character(),
  "effect_allele.outcome" = as.character(),
  "other_allele,outcome"  = as.character(),
  "samplesize.outcome" = as.numeric(),
  "beta.outcome" = as.numeric(),
  "se.outcome" = as.numeric(),
  "eaf.outcome" = as.numeric(),
  "pval.outcome" = as.numeric()
)
snp_outcome_females <- snp_outcome_males

for(i in snps_exposure){
  # Select snp of interest
  data <- genetics |>
    select("snp" = all_of(i), all_of(paste0(i,"_0")),"sex","age_when_assessment",
           starts_with("PC"), "genetic_batch","ad_status") |>
    filter(!is.na(snp))
  
  data_males <- data |> filter(sex == "Male")
  
  data_females <- data |> filter(sex == "Female")
  
  regression_males <- glm(ad_status ~ snp + age_when_assessment + genetic_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                          data = data_males,
                          family = "binomial")
  
  regression_females <- glm(ad_status ~ snp + age_when_assessment + genetic_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                            data = data_females,
                            family = "binomial")
  
  regression_males <- coefficients(summary(regression_males))[row.names(coefficients(summary(regression_males))) == "snp",]
  regression_females <- coefficients(summary(regression_females))[row.names(coefficients(summary(regression_females))) == "snp",]
  
  snp_outcome_males <- snp_outcome_males |>
    rbind(loadSnpOutcomeResults(i, data = data_males, regression = regression_males))
  snp_outcome_females <- snp_outcome_females |>
    rbind(loadSnpOutcomeResults(i, data = data_females, regression = regression_females))
}


write_delim(snp_outcome_males |>
              mutate(type = "males") |>
              rbind(
                snp_outcome_females  |>
                  mutate(type = "females")),
            paste0("iv_outcome.tsv"), delim = "\t")