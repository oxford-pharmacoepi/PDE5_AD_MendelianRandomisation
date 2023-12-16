rm(list = ls())
pacman::p_load('readr','devtools','dplyr','TwoSampleMR','LDlinkR','forestplot','tidyr')
devtools::install_github("bar-woolf/TwoStepCisMR")
library(TwoStepCisMR)
pathToData <- 'C:/Users/marta/Desktop/PDE5_AD_MendelianRandomisation/'

# Mendelian randomisation ------------------------------------------------------
beta <- numeric()
beta_s <- numeric()
se   <- numeric()
cu   <- numeric()
cl   <- numeric()
j <- 1

dec <- c(-5.5,-8.4)
dec <- c(1,1)
for (i in c('DBP','SBP')){
  # Snp - exposure
  exp <- read.table(paste0('iv_',i,'.txt')) %>%
    mutate(id.exposure = 'exposure',
           exposure = 'exposure') %>%
    rename(beta.exposure = beta, se.exposure = se, effect_allele.exposure = effect_allele,
           other_allele.exposure = other_allele, eaf.exposure = eaf) %>%
    mutate(beta.exposure = beta.exposure)
  # Snp - outcome
  out <- read.table(paste0('iv_AlzD_',i,'.txt')) %>%
    mutate(id.outcome = 'outcome',
           outcome = 'outcome',
           eaf.outcome = NA) %>%
    rename(beta.outcome = beta, se.outcome = se, effect_allele.outcome = effect_allele,
           other_allele.outcome = other_allele)
  
  # Mendelian randomisation
  dat <- harmonise_data(exp,out)
  res <- as.data.frame(IVWcorrel(betaYG = dat$beta.outcome,  sebetaYG = dat$se.outcome,
                                 betaXG = dat$beta.exposure, sebetaXG = dat$se.exposure,
                                 rho = ld_matrix(dat$SNP)))
  

  beta[j]   <- exp(res$beta_IVWcorrel*dec[j]) 
  se[j] <- res$se_IVWcorrel.random*dec[j]
  cl[j] <- exp(dec[j]*(res$beta_IVWcorrel + 1.96*res$se_IVWcorrel.random))
  cu[j] <- exp(dec[j]*(res$beta_IVWcorrel - 1.96*res$se_IVWcorrel.random))
  j <- j + 1
}

t <- data.frame(OR = beta, SE = se, upper = cu, lower = cl, Exposure = c('Diastolic blood pressure','Systolic blood pressure'))
write.csv(t, 'mr_results_scaled.csv')


#Two-step MR -------------------------------------------------------------------
beta <- numeric()
se   <- numeric()
cu   <- numeric()
cl   <- numeric()
num <- 1
ex  <- 1
i <- c("ukb-b-19953", # BMI
       'ukb-b-7376', # Impedance of leg (right)
       'ukb-b-14068', # Impedance of leg (left)  
       'ukb-b-7859', # Impedance of arm (right)
       'ukb-b-19379', # Impedance of arm (left)
       'ukb-b-19921', # Impedance of whole body
       'ukb-b-10787', # Standing height
       "ebi-a-GCST004607", # Plateletcrit
       'ebi-a-GCST004626', # Myeloid white cell count
       'ukb-d-30080_irnt', # Platelet count
       'ieu-b-30', # white blood cell count       
       'ebi-a-GCST005195', # Coronary artery disease
       'ebi-a-GCST004614', # Granulocyte count
       'ebi-a-GCST004620') # Sum basophil neutrophil counts

for(j in c('DBP','SBP')){
  for (jj in i){
    # Snp - exposure 
    exp <- read.table(paste0('iv_',j,'.txt')) %>%
      mutate(id.exposure = 'exposure',
             exposure = 'exposure') %>%
      rename(beta.exposure = beta, se.exposure = se, effect_allele.exposure = effect_allele,
             other_allele.exposure = other_allele, eaf.exposure = eaf, pval.exposure = pval) %>%
      mutate(beta.exposure = beta.exposure)
    # Snp - confounder
    conf <- extract_outcome_data(exp$SNP, jj, proxies = F)
    conf <- harmonise_data(exp,conf)
    
    # Snp - outcome
    out <- read.table(paste0('iv_AlzD_',j,'.txt')) %>%
      mutate(id.outcome = 'outcome',
             outcome = 'outcome',
             eaf.outcome = NA) %>%
      rename(beta.outcome = beta, se.outcome = se, effect_allele.outcome = effect_allele,
             other_allele.outcome = other_allele, pval.outcome = pval)
    out <- harmonise_data(exp, out)
    
    # Confounder-outcome associations
    # - Confounder instruments
    conf_instruments <- extract_instruments(outcomes = jj, # Confounder ID
                                            p1 = 5e-08, # Significance threshold
                                            clump = T, # Clumping
                                            p2 = 5e-08, # Secondary clumping threshold
                                            r2 = 0.001, # Clumping r2 clut off
                                            kb = 10000, # Clumping distance cutoff
                                            access_token = ieugwasr::check_access_token(),
                                            force_server = FALSE)
    # - Confounder instruments on the outcome
    out_conf_instruments <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_Lambert_IGAP_2013/IGAP_summary_statistics/IGAP_stage_1.txt'))) 
    out_conf_instruments <- out_conf_instruments %>% 
      filter(MarkerName %in% conf_instruments$SNP) %>%
      rename(effect_allele.outcome = Effect_allele,
             other_allele.outcome = Non_Effect_allele,
             beta.outcome = Beta,
             se.outcome = SE,
             pvalue.outcome = Pvalue,
             SNP = MarkerName) %>%
      mutate(id.outcome = 'AD', outcome = 'AD', eaf.outcome = NA)
    # - Harmonise data
    conf_out <- harmonise_data(conf_instruments,out_conf_instruments, action = 2)
    # - Mendelian randomisation
    conf_out_mr <- TwoSampleMR::mr(conf_out,
                                   method_list=c("mr_ivw","mr_egger_regression",
                                                 "mr_weighted_median", "mr_weighted_mode"))
    
    # Two step mendelian randomisation
    for (snp in out$SNP){ # Each one of the instruments
      res <- TwoStepCisMR::TSCMR(Bgo  = out$beta.outcome[out$SNP == snp], # exp-outcome
                                 SEgo = out$se.outcome[out$SNP == snp],
                                 Bgc  = conf$beta.outcome[conf$SNP == snp], # exp - confounder
                                 SEgc = conf$se.outcome[conf$SNP == snp],
                                 Bco  = conf_out_mr$b[conf_out_mr$method == "Inverse variance weighted"], # Confounder - outcome
                                 SEco = conf_out_mr$se[conf_out_mr$method == "Inverse variance weighted"]) 
      
      out$beta.outcome[out$SNP == snp] <- res$Bgo
      out$se.outcome[out$SNP == snp]   <- res$BSSE
    } 
    
    IVW <- IVWcorrel(betaYG   = out$beta.outcome,
                     sebetaYG = out$se.outcome,
                     betaXG   = out$beta.exposure,
                     sebetaXG = out$se.exposure,
                     rho = ld_matrix(out$SNP))
    
    IVW <- data.frame(IVW)
    
    beta[num] <- exp(IVW$beta_IVWcorrel*dec[ex])
    se[num]   <- IVW$se_IVWcorrel.random
    cl[num]   <- exp(dec[ex]*(IVW$beta_IVWcorrel+1.96*IVW$se_IVWcorrel.random))
    cu[num]   <- exp(dec[ex]*(IVW$beta_IVWcorrel-1.96*IVW$se_IVWcorrel.random))
    
    num <- num + 1
  }
  
  t <- data.frame(i,OR = beta, SE = se, upper = cu, lower = cl)
  write.csv(t, paste0('mr_results_twoStep_',j,'_scaled.csv'))
  ex <- ex+1
}


