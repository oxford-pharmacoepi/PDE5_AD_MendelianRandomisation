
exposure <- c("DBP")
outcome  <- c("Lambert","Wightman","deRojas","Bellenguez","Kunkle")

exposure <- c("new")
outcome  <- c("new_Wightman","new_finngen")

rm(harmonised)
rm(MR_result)
for(exposure_i in exposure){
  # ld matrix
  ld_mat <- read.table(paste0(pathResults,"InstrumentSelection/ld_matrix_",exposure_i,".txt"))
  
  pattern <- ld_mat %>%
    tibble::as_tibble() %>%
    colnames() %>%
    tibble::tibble() %>%
    tidyr::separate(col = ".", into = c("SNP", "effect_allele","other_allele"), sep = "_")
  
  # Snp - exposure
  exp <- read.table(paste0(pathResults,"InstrumentSelection/iv_",exposure_i,".txt"), header = TRUE) %>%
    dplyr::mutate(beta.exposure = -5.5*beta.exposure)
  
  for(outcome_i in outcome){
    # snp - outcome
    out <- read.table(paste0(pathResults,"InstrumentSelection/iv_",outcome_i,".txt"), header = TRUE) %>%
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
    
    if(!"harmonised" %in% ls()){
      # harmonised <- dat$x %>%
      #   dplyr::mutate("outcome" = outcome_i) |>
      #   dplyr::select("SNP", starts_with(c("effect","other","beta","se","eaf",
      #                                      "chr","id","samplesize","pval")))
      harmonised <- 1
      MR_result  <- res 
        
    }else{
      # harmonised <- harmonised %>%
      #   dplyr::union_all(
      #     dat$x %>%
      #       dplyr::mutate("outcome" = outcome_i) |>
      #       dplyr::select("SNP", starts_with(c("effect","other","beta","se","eaf",
      #                                          "chr","id","samplesize","pval")))
      #   )

      MR_result <- MR_result %>%
        dplyr::union_all(res)
    }
  }

  readr::write_delim(harmonised, paste0(pathResults,"MR_Results/harmonised_",exposure_i,"_scaled.txt"))
  readr::write_delim(MR_result, paste0(pathResults,"MR_Results/MR_Results_",exposure_i,"_scaled.txt"))
}



