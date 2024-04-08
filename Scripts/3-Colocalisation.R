# Colocalisation
source(here::here("Functions/loadGwas.R"))

chr <- 4
window <- 0 # Stringent window for the analysis
gene_start <- 120415550
gene_end   <- 120550146 # from ensembl

exposure <- "DBP"
outcome <- c("Lambert","Wightman","deRojas","Bellenguez")



for(exposure_i in exposure){
  exp <- loadGwas(exposure_i, onlyInstruments = FALSE) %>%
    dplyr::filter(chr == 4, pos >= gene_start-window, pos <= gene_end+window)
  
  
  for(outcome_i in outcome){
    out <- loadGwas(outcome_i, onlyInstruments = FALSE) %>%
      dplyr::filter(chr == 4, pos >= gene_start-window, pos <= gene_end+window)
    
    # harmonise data
    dat <- exp %>% 
      dplyr::inner_join(out, by = c("chr", "pos")) %>%
      dplyr::mutate(beta.outcome = dplyr::if_else(effect_allele.exposure == effect_allele.outcome, beta.outcome, -beta.outcome),
                    eaf.outcome  = dplyr::if_else(effect_allele.exposure == effect_allele.outcome, eaf.outcome, 1-eaf.outcome),
                    effect_allele.outcome = dplyr::if_else(effect_allele.exposure == effect_allele.outcome, effect_allele.outcome, other_allele.outcome),
                    other_allele.outcome  = dplyr::if_else(other_allele.exposure  == other_allele.outcome,  other_allele.outcome,  other_allele.exposure)) %>%
      dplyr::group_by(pos) %>% dplyr::filter(dplyr::n() == 1) %>% dplyr::ungroup()
    
    # List exposure
    exp_list <- list(
      beta = dat$beta.exposure,
      MAF  = dat$eaf.exposure,
      pvalues = dat$pval.exposure,
      varbeta = dat$se.exposure^2,
      N       = dat$samplesize.exposure,
      type    = "quant",
      pos     = dat$pos,
      chr     = dat$chr,
      snp     = paste0(dat$chr,":",dat$pos)
    )
    
    # List outcome
    if(is.na(dat$eaf.outcome) %>% unique()){
    
      out_list <- list(
        beta = dat$beta.outcome,
        pvalues = dat$pval.outcome,
        varbeta = dat$se.outcome^2,
        N       = dat$samplesize.outcome,
        type    = "cc",
        pos     = dat$pos,
        chr     = dat$chr,
        snp     = paste0(dat$chr,":",dat$pos)
      )
    }else{
      out_list <- list(
        beta = dat$beta.outcome,
        MAF  = dat$eaf.outcome,
        pvalues = dat$pval.outcome,
        varbeta = dat$se.outcome^2,
        N       = dat$samplesize.outcome,
        type    = "cc",
        pos     = dat$pos,
        chr     = dat$chr,
        snp     = paste0(dat$chr,":",dat$pos)
      )
    }
    
    res <- coloc::coloc.abf(exp_list, out_list)
    
    coloc_tab <- tibble::tibble(
      nsnps = res$summary[1],
      H0    = res$summary[2],
      H1    = res$summary[3],
      H2    = res$summary[4],
      H3    = res$summary[5],
      H4    = res$summary[6],
      P1    = 1e-4,
      P2    = 1e-4,
      P12   = 1e-5
    ) %>%
      dplyr::mutate(
        exposure = .env$exposure_i,
        outcome  = .env$outcome_i
      )
    
    
    if("coloc_table" %in% ls()){
     coloc_table <- coloc_table %>% 
       dplyr::union_all(coloc_tab)
    }else{
      coloc_table <- coloc_tab
    }
  }
}

readr::write_delim(coloc_table, paste0(pathResults,"Colocalisation/colocalisation.txt"))
