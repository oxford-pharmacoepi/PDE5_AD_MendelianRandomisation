#Two-step MR
rm(list = ls())
pacman::p_load('readr','devtools','dplyr','TwoSampleMR')
devtools::install_github("bar-woolf/TwoStepCisMR")
library(TwoStepCisMR)
pathToData <- 'C:/Users/martaa/Desktop/Projects/PDE5_AD_MendelianRandomisation/'

i <- "ukb-b-19953" # BMI


exp <- read.table('iv_DBP.txt') %>%
  mutate(id.exposure = 'exposure',
         exposure = 'exposure') %>%
  rename(beta.exposure = beta, se.exposure = se, effect_allele.exposure = effect_allele,
         other_allele.exposure = other_allele, eaf.exposure = eaf)

out <- read.table('iv_AlzD_DBP.txt') %>%
  mutate(id.outcome = 'outcome',
         outcome = 'outcome',
         eaf.outcome = NA) %>%
  rename(beta.outcome = beta, se.outcome = se, effect_allele.outcome = effect_allele,
         other_allele.outcome = other_allele)

conf <- extract_outcome_data(exp$SNP,"ukb-b-19953", proxies = F)

conf_dbp <- harmonise_data(exp,conf)
kinds_men_dbp  <- harmonise_data(exp,out)

## Confounder outcome assoicaiton
conf <- extract_instruments(outcomes= i, p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000, access_token = ieugwasr::check_access_token(), force_server = FALSE)

conf_kids<- extract_outcome_data(conf$SNP, 'ieu-a-297', proxies = T, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = ieugwasr::check_access_token(), splitsize = 10000, proxy_splitsize = 500)
conf_kids<- harmonise_data(conf, conf_kids, action = 2)
conf_kids_mr <- TwoSampleMR::mr(conf_kids, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
conf_kids_mr<- merge(conf_kids_mr,mr_heterogeneity(conf_kids), by=c("id.outcome", "method", "id.exposure", "outcome","exposure") , all=TRUE )


# ## Confounder outcome association
# conf1 <- extract_instruments(outcomes = "ukb-b-19953", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000,
#                              access_token = ieugwasr::check_access_token(), force_server = FALSE)
# 
# 
# out1 <- extract_outcome_data(conf1$SNP, 'ieu-a-297', proxies = T, rsq = 0.8, align_alleles = 1, palindromes = 1,
#                                  maf_threshold = 0.3, access_token = ieugwasr::check_access_token(), splitsize = 10000,
#                                  proxy_splitsize = 500)
# 
# conf1_out1 <- harmonise_data(conf1, out1, action=2)
# 
# conf1_out1_mr <- TwoSampleMR::mr(conf1_out1, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
# conf1_out1_mr <- merge(conf1_out1_mr,mr_heterogeneity(conf1_out1), by=c("id.outcome", "method", "id.exposure", "outcome","exposure") , all=TRUE)
# 
# 
# dat <- out_exp
# if (conf1_out1_mr$Q_pval[1] < 0.05){
#   for (sni in dat:SNP){
#     a <- TSCMR(Bgo=dat$beta.outcome[dat$SNP==snp],
#                SEgo = dat$se.outcome[dat$SNP==snp],               Bgc = conf_exp$beta.outcome[conf_exp$SNP==snp],
#                SEgc = dat$conf_exp$se.outcome[conf_exp$SNP==snp], Bco = conf1_out1_mr$b[conf1_kids1_mr$method=="Weighted median"],)
#   }
# }
