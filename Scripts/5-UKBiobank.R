# Load data
snps_outcome <- as_tibble(read.delim(paste0(pathResults,"UK Biobank/iv_outcome.tsv")))

# Extract instruments ----
snps_exposure <- read.table(paste0(pathResults,"InstrumentSelection/iv_DBP.txt"))

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
  out <- snps_outcome |> 
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
