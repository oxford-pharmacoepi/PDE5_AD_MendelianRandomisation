exposure_i <- c("DBP")

# ld matrix
ld_mat <- read.table(paste0(pathResults,"InstrumentSelection/ld_matrix_",exposure_i,".txt"))

pattern <- ld_mat %>%
  tibble::as_tibble() %>%
  colnames() %>%
  tibble::tibble() %>%
  tidyr::separate(col = ".", into = c("SNP", "effect_allele","other_allele"), sep = "_")

# Snp - exposure
exp <- read.table(paste0(pathResults,"InstrumentSelection/iv_",exposure_i,".txt")) %>%
  dplyr::mutate(beta.exposure = -5.5*beta.exposure) |>
  dplyr::mutate(se.exposure = 5.5*se.exposure)
rm("MR_result")
for(outcome_i in c("Lambert","Wightman","deRojas","Bellenguez")){
  # snp - outcome
  out <- readr::read_table(paste0(pathResults,"InstrumentSelection/iv_",outcome_i,".txt")) %>%
    dplyr::inner_join(exp %>% dplyr::select("SNP"),
                      by = "SNP")
    
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
    
    for(i in c(pattern$SNP,"All")){
    dat_loo <- dat$x %>%
      dplyr::filter(SNP != i)
    ld_loo  <- dat$ld[!colnames(dat$ld) == i,]
    ld_loo  <- ld_loo[,!colnames(ld_loo) == i]

    res <- TwoStepCisMR::IVWcorrel(betaYG   = dat_loo$beta.outcome,
                                   sebetaYG = dat_loo$se.outcome,
                                   betaXG   = dat_loo$beta.exposure,
                                   sebetaXG = dat_loo$se.exposure,
                                   rho = ld_loo) %>%
      as.data.frame() %>%
      dplyr::rename("beta" = "beta_IVWcorrel",
                    "se"   = "se_IVWcorrel.random",
                    "pval" = "p_IVWcorrel.random",
                    "instruments" = "n_snp") %>%
      dplyr::mutate("outcome" = outcome_i,
                    "OR" = exp(beta),
                    "cilow" = exp(beta-1.96*se),
                    "cihigh" = exp(beta+1.96*se),
                    "SNP" = i)
    
    if(!"MR_result" %in% ls()){
      MR_result  <- res 
      
    }else{
      MR_result <- MR_result %>%
        dplyr::union_all(res)
    }
  }
}

readr::write_delim(MR_result, paste0(pathResults,"LeaveOneOutAnalysis/LeaveOneOut_Results_",exposure_i,".txt"))



