rm(list = ls())
pacman::p_load('readr','devtools','dplyr','TwoSampleMR','LDlinkR','forestplot','tidyr')
devtools::install_github("bar-woolf/TwoStepCisMR")
library(TwoStepCisMR)
pathToData <- 'C:/Users/martaa/Desktop/Projects/PDE5_AD_MendelianRandomisation/'

# Mendelian randomisation ------------------------------------------------------
beta <- numeric()
se   <- numeric()
cu   <- numeric()
cl   <- numeric()
j <- 1
for (i in c('DBP','SBP')){
  exp <- read.table(paste0('iv_',i,'.txt')) %>%
    mutate(id.exposure = 'exposure',
           exposure = 'exposure') %>%
    rename(beta.exposure = beta, se.exposure = se, effect_allele.exposure = effect_allele,
           other_allele.exposure = other_allele, eaf.exposure = eaf)
  
  out <- read.table(paste0('iv_AlzD_',i,'.txt')) %>%
    mutate(id.outcome = 'outcome',
           outcome = 'outcome',
           eaf.outcome = NA) %>%
    rename(beta.outcome = beta, se.outcome = se, effect_allele.outcome = effect_allele,
           other_allele.outcome = other_allele)
  
  dat   <- harmonise_data(exp,out)
  res <- as.data.frame(IVWcorrel(betaYG = dat$beta.outcome,  sebetaYG = dat$se.outcome,
                                 betaXG = dat$beta.exposure, sebetaXG = dat$se.exposure,
                                 rho = ld_matrix(dat$SNP)))
  
  beta[j] <- exp(res$beta_IVWcorrel)
  se[j] <- res$se_IVWcorrel.random
  cu[j] <- exp(res$beta_IVWcorrel + 1.96*res$se_IVWcorrel.random)
  cl[j] <- exp(res$beta_IVWcorrel - 1.96*res$se_IVWcorrel.random)
  j <- j + 1
}

t <- data.frame(OR = beta, SE = se, upper = cu, lower = cl, Exposure = c('Diastolic blood pressure','Systolic blood pressure'))
write.csv(t, 'mr_results.csv')



#Two-step MR -------------------------------------------------------------------
beta <- numeric()
se   <- numeric()
cu   <- numeric()
cl   <- numeric()
num <- 1
i <- "ukb-b-19953" # BMI
k <- 'ieu-a-297' # Alzheimer

for(j in c('DBP','SBP')){
  exp_dsbp <- read.table(paste0('iv_',j,'.txt')) %>%
    mutate(id.exposure = 'exposure',
           exposure = 'exposure') %>%
    rename(beta.exposure = beta, se.exposure = se, effect_allele.exposure = effect_allele,
           other_allele.exposure = other_allele, eaf.exposure = eaf)
  conf_dbp <- extract_outcome_data(exp_dsbp$SNP, i, proxies = F)
  conf_dbp <- harmonise_data(exp_dsbp,conf_dbp)
  
  kids_men_dbp <- extract_outcome_data(exp_dsbp$SNP, k, proxies = F)
  kids_men_dbp <- harmonise_data(exp_dsbp, kids_men_dbp)
  
  # Confounder-outcome associations
  # Search for GWAS significant SNPs. It performs LD-Clumping to return only 
  # independent significant associations.
  
  # Snps associated with confounder:
  conf <- extract_instruments(outcomes = i, # Confounder ID
                              p1 = 5e-08, # Significance threshold
                              clump = T, # Clumping
                              p2=5e-08, # Secondary clumping threshold
                              r2 = 0.001, # Clumping r2 clut off
                              kb = 10000, # Clumping distance cutoff
                              access_token = ieugwasr::check_access_token(),
                              force_server = FALSE)
  # Confounder instruments, effect size on the outcome
  conf_kids <- extract_outcome_data(conf$SNP, # Snps ID
                                    k, # Outcome
                                    proxies = T, # LD tags
                                    rsq = 0.8, # Minimum LD rsq value
                                    align_alleles = 1, # Align tag alleles to target alleles
                                    palindromes = 1, # Allow palindromic snps
                                    maf_threshold = 0.3, # MAF threshold to try to infern palindromic SNPs
                                    access_token = ieugwasr::check_access_token(), 
                                    splitsize = 10000,
                                    proxy_splitsize = 500)
  conf_kids <- harmonise_data(conf, conf_kids, action = 2)
  conf_kids_mr <- TwoSampleMR::mr(conf_kids, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
  conf_kids_mr<-merge(conf_kids_mr,mr_heterogeneity(conf_kids), by=c("id.outcome", "method", "id.exposure", "outcome","exposure") , all=TRUE )
  
  dat<-kids_men_dbp # exp - outcome association
  if (conf_kids_mr$Q_pval[1] < 0.05){ # if there is no heterogeneity
    for (snp in dat$SNP){
      a <- TSCMR(Bgo=dat$beta.outcome[dat$SNP==snp], # Byz: Beta instrument exp-outcome (SNP-)
                 SEgo=dat$se.outcome[dat$SNP==snp],  # SEyz: Se instrument exp-outcome
                 Bgc=conf_dbp$beta.outcome[conf_dbp$SNP==snp], # Bkz: Beta instrument exp-outcome in exp-confunder
                 SEgc=conf_dbp$se.outcome[conf_dbp$SNP==snp], # SEkz: Se instrument exp-outcome in dbp-confounder
                 Bco=conf_kids_mr$b[conf_kids_mr$method=="Weighted median"], # MR result of weighted median
                 SEco=conf_kids_mr$se[conf_kids_mr$method=="Weighted median"]) # MR result of weighted median
      dat$beta.outcome[dat$SNP==snp]<-a$Bgo # New Byz
      dat$se.outcome[dat$SNP==snp]<-a$BSSE # New SEyz
    }
  }else {
    for (snp in dat$SNP){
      a<-TSCMR(Bgo=dat$beta.outcome[dat$SNP==snp], 
               SEgo=dat$se.outcome[dat$SNP==snp],
               Bgc=conf_dbp$beta.outcome[conf_dbp$SNP==snp], 
               SEgc=conf_dbp$se.outcome[conf_dbp$SNP==snp],
               Bco=conf_kids_mr$b[conf_kids_mr$method=="Inverse variance weighted"],
               SEco=conf_kids_mr$se[conf_kids_mr$method=="Inverse variance weighted"])
      dat$beta.outcome[dat$SNP==snp]<-a$Bgo #Inverse variance weighted
      dat$se.outcome[dat$SNP==snp]<-a$BSSE
    }
  }
  
  IVW <- data.frame(IVWcorrel(betaYG=dat$beta.outcome,  
                 sebetaYG=dat$se.outcome,  
                 betaXG=dat$beta.exposure, 
                 sebetaXG=dat$se.exposure, 
                 rho=ld_matrix(dat$SNP)))
  beta[num] <- exp(IVW$beta_IVWcorrel)
  se[num]   <- IVW$se_IVWcorrel.random
  cu[num]   <- exp(IVW$beta_IVWcorrel+1.96*IVW$se_IVWcorrel.random)
  cl[num]   <- exp(IVW$beta_IVWcorrel-1.96*IVW$se_IVWcorrel.random)
  num <- num + 1
}

t <- data.frame(OR = beta, SE = se, upper = cu, lower = cl, Exposure = c('Diastolic blood pressure','Systolic blood pressure'))
write.csv(t, 'mr_results_twoStep.csv')


# Colocalization
chr <- 4
window <- 0 #using a more stringent window for the anlysis. 
gene_start <- 120415550#119494395#
gene_end <- 120550146#119628991#
chrpos <- paste0(chr, ":", gene_start - window, "-", gene_end + window)

dbp <- as_tibble(read_table(paste0(pathToData,'GWAS/BP_Evangelou_2018/Evangelou_30224653_DBP.txt.gz'))) %>%
    separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
    mutate(pos = as.numeric(pos)) %>%
  filter(pos > gene_start-window & pos < gene_end+window & chr == chr) 

sbp <- as_tibble(read_table(paste0(pathToData,'GWAS/BP_Evangelou_2018/Evangelou_30224653_SBP.txt.gz'))) %>%
  separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
  mutate(pos = as.numeric(pos)) %>%
  right_join(dbp %>% select(chr, pos), by = c('chr','pos'))
  
alz <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_Lambert_IGAP_2013/IGAP_summary_statistics/IGAP_stage_1.txt'))) %>%
  rename('chr' = 'Chromosome', 'pos' = 'Position') %>%
  right_join(dbp %>% mutate(chr = as.numeric(chr)) %>% select(chr, pos), by = c('chr','pos'))

# Prove that all gwas have same length -> same number of snps
nrow(dbp) == nrow(sbp)
nrow(sbp) == nrow(alz)

library(hyprcoloc)
hyper <- dbp %>% select(chr, pos, Beta_dbp = Effect, se_dbp = StdErr, p_dbp=DBP, maf_dbp = Freq1) %>%
  left_join(sbp %>% select(chr, pos, beta_sbp = Effect, se_sbp = StdErr, p_sbp=DBP, maf_sbp = Freq1), by = c('chr','pos')) %>%
  left_join(alz %>% select(chr,pos,snp=MarkerName, beta_alz = Beta, se_alz = SE, p_alz = Pvalue) %>% mutate(maf = 0.5), by = c('chr','pos'))

hyper<-merge(dsbp %>% 
               mutate(se_dbp = varbeta) %>%
               select("snp","beta_dbp" = "beta","se_dbp","pvalues_dbp" = pvalues),
             sbp %>% mutate(se_sbp = varbeta) %>%
               select("snp","beta_sbp" = "beta","se_sbp","pvalues_sbp" = pvalues),by="snp")
hyper<-merge(hyper,kids %>%
               mutate(se_kids= varbeta) %>%
               select("snp","beta_kids" = "beta","se_kids","pvalues_kids" = pvalues),by="snp")

betas <- as.matrix(hyper[,grepl("beta_", names(hyper))])
ses <- as.matrix(hyper[,grepl("se_", names(hyper))])
hyprcoloc_results <- hyprcoloc::hyprcoloc(effect.est = betas,
                                          effect.se = ses,
                                          trait.names = c("sbp", "dsbp","number of childern fathered" ),
                                          snp.id = hyper$snp,
                                          binary.outcomes = c(0, 1), 
                                          prior.1 = 1e-04, prior.c = 0.02)
hyprcoloc_results$results
