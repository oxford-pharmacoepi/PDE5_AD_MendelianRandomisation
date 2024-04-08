# Two step cis-MR
source(here::here("Functions/loadGwas.R"))
dir.create(pathResults,"TwoStepMR/")
confounders <- c("ukb-b-19953", # BMI - https://gwas.mrcieu.ac.uk/datasets/ukb-b-19953/, HG19/GRCh37
                 'ukb-b-7376', # Impedance of leg (right) - https://gwas.mrcieu.ac.uk/datasets/ukb-b-7376/, HG19/GRCh37
                 'ukb-b-14068', # Impedance of leg (left) - https://gwas.mrcieu.ac.uk/datasets/ukb-b-14068/, HG19/GRCh37
                 'ukb-b-7859', # Impedance of arm (right) - https://gwas.mrcieu.ac.uk/datasets/ukb-b-7859/, HG19/GRCh37
                 'ukb-b-19379', # Impedance of arm (left) - https://gwas.mrcieu.ac.uk/datasets/ukb-b-19379/, HG19/GRCh37
                 'ukb-b-19921', # Impedance of whole body - https://gwas.mrcieu.ac.uk/datasets/ukb-b-19921/, HG19/GRCh37
                 'ukb-b-10787', # Standing height - https://gwas.mrcieu.ac.uk/datasets/ukb-b-10787/, HG19/GRCh37
                 "ebi-a-GCST004607", # Plateletcrit - https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004607/, HG19/GRCh37
                 'ebi-a-GCST004626', # Myeloid white cell count - https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004626/, HG19/GRCh37
                 'ukb-d-30080_irnt', # Platelet count - https://gwas.mrcieu.ac.uk/datasets/ukb-d-30080_irnt/, HG19/GRCh37
                 'ieu-b-30', # white blood cell count - https://gwas.mrcieu.ac.uk/datasets/ieu-b-30/, 	HG19/GRCh37
                 'ebi-a-GCST005195', # Coronary artery disease - https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST005195/, HG19/GRCh37
                 'ebi-a-GCST004614', # Granulocyte count - https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004614/, 	HG19/GRCh37
                 'ebi-a-GCST004620') # Sum basophil neutrophil counts - https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004620/, HG19/GRCh37


for(exposure_i in "DBP"){
  # Harmonisation
  print("line 22")
  ld_mat <- read.table(paste0(pathResults,"InstrumentSelection/ld_matrix.txt"))
  print("line 25")
  pattern <- ld_mat %>%
    tibble::as_tibble() %>%
    colnames() %>%
    tibble::tibble() %>%
    tidyr::separate(col = ".", into = c("SNP", "effect_allele","other_allele"), sep = "_")
  
  for(confounder_i in confounders){
    # SNP - exposure
    print(paste0(confounder_i," line 34"))
    exp <- read.table(paste0(pathResults,"InstrumentSelection/iv_",exposure_i,".txt"))
    
    # SNP - confounder
    print("line 38")
    conf <- TwoSampleMR::extract_outcome_data(exp$SNP,confounder_i, proxies = FALSE)
    conf <- TwoSampleMR::harmonise_data(exp,conf)
    
    for(outcome_i in c("Lambert","Wightman","deRojas","Bellenguez")){
      # SNP - outcome
      print(paste0(outcome_i,"line 44"))
      out <- read.table(paste0(pathResults,"InstrumentSelection/iv_",outcome_i,".txt"), header = TRUE)
      out <- TwoSampleMR::harmonise_data(exp, out)
      
      # Confounder - outcome
      print("line 49")
      conf_out_conf <- TwoSampleMR::extract_instruments(outcomes = confounder_i, p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000, force_server = FALSE)
      print("line 51")
      if(outcome_i == "Wightman"){
        conf_out_out  <- loadGwas(outcome_i, onlyInstruments = FALSE) %>% dplyr::inner_join(conf_out_conf %>% 
                                                                                              dplyr::select("chr" = "chr.exposure", "pos" = "pos.exposure", "SNP") %>% 
                                                                                              dplyr::mutate(chr = as.numeric(chr), pos = as.numeric(pos)))
      }else{
        conf_out_out  <- loadGwas(outcome_i, onlyInstruments = FALSE, bellenguez = FALSE) %>%
          dplyr::filter(SNP %in% conf_out_conf$SNP)
      }
      print("line 60")
      conf_out    <- TwoSampleMR::harmonise_data(conf_out_conf, conf_out_out, action = 2)
      conf_out_mr <- TwoSampleMR::mr(conf_out, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
      print("line 62")
      conf_out_mr <- conf_out_mr %>%
        dplyr::full_join(
          TwoSampleMR::mr_heterogeneity(conf_out),
          by = c("id.outcome", "method", "id.exposure", "outcome", "exposure")
        )
      print("line 68")
      if(conf_out_mr$Q_pval[1] < 0.05){
        for(snp in out$SNP){
          a <- TwoStepCisMR::TSCMR(
            # SNP - outcome association
            Bgo  = out %>% dplyr::filter(SNP == snp) %>% dplyr::pull(beta.outcome),
            SEgo = out %>% dplyr::filter(SNP == snp) %>% dplyr::pull(se.outcome),
            # SNP - confounder association
            Bgc  = conf %>% dplyr::filter(SNP == snp) %>% dplyr::pull(beta.outcome),
            SEgc = conf %>% dplyr::filter(SNP == snp) %>% dplyr::pull(se.outcome),
            # Counder - outcome association
            Bco  = conf_out_mr %>% dplyr::filter(method == "Weighted median") %>% dplyr::pull(b),
            SEco = conf_out_mr %>% dplyr::filter(method == "Weighted median") %>% dplyr::pull(se)
          )
          out$beta.outcome[out$SNP == snp] <- a$Bgo
          out$se.outcome[out$SNP == snp]   <- a$BSSE
        }
      }else{
        a <- TwoStepCisMR::TSCMR(
          Bgo  = out %>% dplyr::filter(SNP == snp) %>% dplyr::pull(beta.outcome),
          SEgo = out %>% dplyr::filter(SNP == snp) %>% dplyr::pull(se.outcome),
          Bgc  = conf %>% dplyr::filter(SNP == snp) %>% dplyr::pull(beta.outcome),
          SEgc = conf %>% dplyr::filter(SNP == snp) %>% dplyr::pull(se.outcome),
          Bco  = conf_out_mr %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(b),
          SEco = conf_out_mr %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::pull(se)
        )
        out$beta.outcome[out$SNP == snp] <- a$Bgo
        out$se.outcome[out$SNP == snp]   <- a$BSSE
      }
      
      # Harmonised 
      print("line 100")
      dat_harmonised <- pattern %>%
        dplyr::left_join(out, by = "SNP") %>%
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
      
      # Two step cis-MR
      print("line 117")
      resTwoStepMR <- TwoStepCisMR::IVWcorrel(
        betaYG   = dat_harmonised$beta.outcome,
        sebetaYG = dat_harmonised$se.outcome,
        betaXG   = dat_harmonised$beta.exposure,
        sebetaXG = dat_harmonised$se.exposure,
        rho      = ld_mat
      ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          exposure   = exposure_i,
          outcome    = outcome_i,
          confounder = confounder_i
        )
      print("line 130")
      if("TwoStepMR_table" %in% ls()){
        TwoStepMR_table <- TwoStepMR_table %>% 
          dplyr::union_all(resTwoStepMR)
      }else{
        TwoStepMR_table <- resTwoStepMR
      }
      
    }
  }
}


readr::write_delim(TwoStepMR_table,paste0(pathResults,"TwoStepMR/twostepMR.txt"))

