# Load data
ukb_data <- loadUkbData(pathUKB)

# Extract instruments ----
snps_exposure <- read.table(paste0(pathResults,"InstrumentSelection/iv_DBP.txt"))

# Calculate SNPS-outcome interaction ----
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

for(i in snps_exposure$SNP){
  # Select snp of interesta
  data <- ukb_data |>
    select("snp" = all_of(i), all_of(paste0(i,"_0")),"sex","age_when_assessment",
           starts_with("PC"), "genetic_batch","ad_status") |>
    filter(!is.na(snp))
  
  data_males <- data |> filter(sex == 1)
  
  data_females <- data |> filter(sex == 0)
  
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


write.csv(snp_outcome_males |>
            mutate(type = "males") |>
            rbind(
              snp_outcome_females  |>
                mutate(type = "females")),
          paste0(pathResults,"SexSpecific/iv_outcome.txt"), row.names = FALSE)



# Calculate MR estimates ----
ld_mat <- read.table(paste0(pathResults,"InstrumentSelection/ld_matrix_dbp.txt"))

pattern <- ld_mat %>%
  tibble::as_tibble() %>%
  colnames() %>%
  tibble::tibble() %>%
  tidyr::separate(col = ".", into = c("SNP", "other_allele","effect_allele"), sep = "_")

# Snp - exposure
exp <- read.table(paste0(pathResults,"InstrumentSelection/iv_dbp.txt"), header = TRUE) %>%
  dplyr::mutate(beta.exposure = -5.5*beta.exposure)

for(outcome_i in c("males","females")){
  # snp - outcome
  out <- read.csv(paste0(pathResults,"SexSpecific/iv_outcome.txt")) |> 
    filter(type == outcome_i) |>
    mutate(other_allele.outcome = c("T","G","T","T","T"), id.outcome = outcome_i, outcome = outcome_i) |>
    rename("SNP" = "snp")
  
  # harmonise data
  dat_harmonised <- pattern %>%
    dplyr::left_join(
      TwoSampleMR::harmonise_data(exposure = exp, outcome  = out),
      by = "SNP"
    ) %>%
    dplyr::mutate(
      beta.exposure = dplyr::if_else(effect_allele.exposure != effect_allele, -beta.exposure, beta.exposure),
      eaf.exposure  = dplyr::if_else(effect_allele.exposure != effect_allele, 1-eaf.exposure, eaf.exposure),
      effect_allele.exposure  = dplyr::if_else(effect_allele.exposure != effect_allele, effect_allele, effect_allele.exposure),
      other_allele.exposure   = dplyr::if_else(other_allele.exposure  != other_allele,  other_allele,  other_allele.exposure),
    ) %>%
    dplyr::mutate(
      beta.outcome = dplyr::if_else(effect_allele.outcome != effect_allele, -beta.outcome, beta.outcome),
      eaf.outcome  = dplyr::if_else(effect_allele.outcome != effect_allele, 1-eaf.outcome, eaf.outcome),
      effect_allele.outcome  = dplyr::if_else(effect_allele.outcome != effect_allele, effect_allele, effect_allele.outcome),
      other_allele.outcome   = dplyr::if_else(other_allele.outcome  != other_allele,  other_allele,  other_allele.outcome),
    )
  
  dat <- TwoSampleMR::harmonise_ld_dat(dat_harmonised,ld_mat)
  res <- TwoStepCisMR::IVWcorrel(betaYG  = dat$x$beta.outcome,
                                 sebetaYG = dat$x$se.outcome,
                                 betaXG   = dat$x$beta.exposure,
                                 sebetaXG = dat$x$se.exposure,
                                 rho = dat$ld) %>%
    as.data.frame() %>%
    dplyr::rename("beta" = "beta_IVWcorrel",
                  "se"   = "se_IVWcorrel.random",
                  "pval" = "p_IVWcorrel.random",
                  "instruments" = "n_snp") %>%
    dplyr::mutate("outcome" = outcome_i,
                  "OR" = exp(beta),
                  "cilow" = exp(beta-1.96*se),
                  "cihigh" = exp(beta+1.96*se))
  
  
  if(outcome_i == "males"){
    MR_result <- res
  }else{
    MR_result <- MR_result |> rbind(res)
  }
}

write.csv(MR_result,
          paste0(pathResults,"SexSpecific/MR_result.txt"), row.names = FALSE)
