library(dplyr)

gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/BP_Evangelou_2018/Evangelou_30224653_SBP.txt/Evangelou_30224653_SBP.txt'))) %>%
  tidyr::separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
  dplyr::mutate(pos = as.numeric(pos),
                chr = as.numeric(chr)) %>%
  dplyr::filter(Type == "SNP") %>%
  dplyr::select(
    "chr", "pos",
    "effect_allele.exposure" = "Allele1",
    "other_allele.exposure"  = "Allele2",
    "pval.exposure"          = "P",
    "samplesize.exposure"    = "TotalSampleSize",
    "beta.exposure"          = "Effect",
    "se.exposure"            = "StdErr",
    "eaf.exposure"           = "Freq1"
  ) %>% 
  dplyr::mutate(
    "exposure" = "SBP",
    "id.exposure" = "SBP"
  ) %>%
  dplyr::mutate(
    effect_allele.exposure = toupper(effect_allele.exposure),
    other_allele.exposure  = toupper(other_allele.exposure)
  )
ED_finngen <- TwoSampleMR::extract_outcome_data(gwas1$SNP,
                                                c("finn-b-ERECTILE_DYSFUNCTION"), 
                                                proxies = F, 
                                                rsq = 0.8, 
                                                align_alleles = 1, 
                                                palindromes = 1,
                                                maf_threshold = 0.3, 
                                                access_token = ieugwasr::check_access_token(),splitsize = 10000,  proxy_splitsize = 500)

gwas3 <- extract_instruments(
  "ieu-b-38",
  p1 = 5e-08,
  clump = FALSE,
  p2 = 5e-08,
  r2 = 1,
  kb = 1,
  access_token = ieugwasr::check_access_token(),
  force_server = FALSE
) |> filter(chr.exposure == 4)

gwas31 <- gwas3 |> filter(pos.exposure >= 120415558-150000 & pos.exposure <= 120549959+150000)

gwas32 <- clump_data(
  gwas31 |> rename(chr_name = chr.exposure, chrom_start = pos.exposure),
  clump_kb = 250000,
  clump_r2 = 0.9,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)


ld <- TwoSampleMR::ld_matrix(gwas31$SNP, with_alleles = TRUE, pop = "EUR") 

ED_finngen <- TwoSampleMR::extract_outcome_data(gwas1$SNP,
                                        c("finn-b-ERECTILE_DYSFUNCTION"), 
                                        proxies = F, 
                                        rsq = 0.8, 
                                        align_alleles = 1, 
                                        palindromes = 1,
                                        maf_threshold = 0.3, 
                                        access_token = ieugwasr::check_access_token(),splitsize = 10000,  proxy_splitsize = 500)

alz <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/AlzD_wightman_2021_excluding_23andme/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz'))) %>%
  dplyr::select(chr           = "chromosome",
                pos           = "base_pair_location",
                effect_allele.outcome = "effect_allele",
                other_allele.outcome  = "other_allele",
                eaf.outcome           = "effect_allele_frequency",
                pval.outcome          = "p_value",
                beta.outcome          = "beta",
                se.outcome            = "standard_error",
                samplesize.outcome    = "N") %>%
  dplyr::mutate(id.outcome = "Wightman",
                outcome    = "Wightman",
                chr = as.numeric(chr))
alz <- alz |>
  inner_join(
    tibble(
      SNP = c("rs4643791", "rs2389873"),
      pos = c(120265619,120554714),
      chr = c(4,4)
    )
  )


# mr analysis-
for(outcome_i in outcome){
  # snp - outcome
  dat <- TwoSampleMR::harmonise_data(gwas31,alz1)

  r1 <- TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, b_out = dat$beta.outcome,se_exp = dat$se.exposure, se_out = dat$se.outcome)
  r2 <- TwoSampleMR::mr_ivw_mre(b_exp = dat$beta.exposure, b_out = dat$beta.outcome,se_exp = dat$se.exposure, se_out = dat$se.outcome)
  
  result <- tibble(
    Outcome = "ad",
    Method = c("mr_ivw","mr_ivw_mre"),
    beta = c(r1$b, r2$b),
    se = c(r1$se, r2$se),
    OR = exp(-beta),
    low = exp(-beta-1.96*se),
    high = exp(-beta+1.96*se),
    pval = c(r1$pval, r2$pval)
    )
  

}

i <- 119494397
f <- 119628804
snp1 <- 119344464 
snp2 <- 119633559 
table <- tibble(x = c(i-150*10^3,f+150*10^3),
                y = 1)
square <- tibble(x = c(i, f, i,f),
                 y = c(0.9,0.9,1.1,1.1))
ggplot(table, aes(x = x, y = y)) +
  geom_segment(aes(x = i-150*10^3, xend = f+150*10^3, y = 1, yend = 1)) +
  geom_segment(aes(x = x, y = y-0.1, xend = x, yend = y+0.1)) +
  scale_y_continuous(limits = c(0.5,1.5)) + 
  theme_classic() +
  scale_x_continuous(limits = c(i-150*10^3, f+150*10^3)) +
  geom_rect(aes(xmin = i, xmax = f, ymin = 0.9, ymax = 1.1), fill = "orange") +
  geom_rect(aes(xmin = snp1-0.1*10^4, xmax = snp1+0.25*10^4, ymin = 0.95, ymax = 1.05), fill = "brown2") + 
  geom_rect(aes(xmin = snp2-0.25*10^4, xmax = snp2+0.25*10^4, ymin = 0.95, ymax = 1.05), fill = "brown2")



# RACER





  