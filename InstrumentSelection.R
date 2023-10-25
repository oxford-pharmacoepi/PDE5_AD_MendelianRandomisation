rm(list=ls())
pacman::p_load('dplyr','readr','here','stringr','tidyr','TwoSampleMR',
               'flextable','ftExtra')

## Read gwas -------------------------------------------------------------------
pathToData <- 'C:/Users/marta/Desktop/PDE5_AD_MendelianRandomisation/'

# Diastolic blood pressure  ====================================================
gwasDBP <- as_tibble(read_table(paste0(pathToData,'GWAS/BP_Evangelou_2018/Evangelou_30224653_DBP.txt.gz'))) %>%
  separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
  mutate(pos = as.numeric(pos))

snpsDBP <- gwasDBP %>%
  filter(chr == 4  & Type =='SNP') %>%
  filter(pos %in% c(120423094,120502461,120416096,120532085,120509279)) %>%
  mutate(SNP = case_when(
    pos == 120423094 ~ "rs80223330",
    pos == 120502461 ~ "rs12646525",
    pos == 120416096 ~ "rs17355550",
    pos == 120532085 ~ "rs10050092",
    pos == 120509279 ~ "rs66887589"
  )) %>%
  select(SNP, chr, pos,
         effect_allele.exposure = "Allele1",
         other_allele.exposure  = "Allele2",
         pval.exposure          = "P",
         samplesize.exposure    = "N_effective",
         beta.exposure          = "Effect",
         se.exposure            = "StdErr",
         eaf.exposure           = 'Freq1')
write.table(snpsDBP,here('Results','Instruments','iv_DBP.txt'))

# Alzheimer disease gwas - Lambert
gwasAlzD_Lambert <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_Lambert_IGAP_2013/IGAP_summary_statistics/IGAP_stage_1.txt')))

snps_Lambert_DBP <- gwasAlzD_Lambert %>%
  filter(MarkerName %in% snpsDBP$SNP) %>%
  select(SNP           = "MarkerName",
         chr           = "Chromosome",
         pos           = "Position", #hg19 
         effect_allele.outcome = "Effect_allele",
         other_allele.outcome  = "Non_Effect_allele",
         pval.outcome          = "Pvalue",
         beta.outcome          = "Beta",
         se.outcome            = "SE") %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome")
write.table(snps_Lambert_DBP,here('Results','Instruments','iv_Lambert_DBP.txt'))

# Alzheimer disease gwas - Wightman
gwasAlzD_Wightman <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_wightman_2021_excluding_23andme/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz')))

snps_Wightman_DBP <- gwasAlzD_Wightman %>%
  filter(chromosome == 4) %>%
  mutate(SNP = case_when(
    base_pair_location == 120423094 ~ "rs80223330",
    base_pair_location == 120502461 ~ "rs12646525",
    base_pair_location == 120416096 ~ "rs17355550",
    base_pair_location == 120532085 ~ "rs10050092",
    base_pair_location == 120509279 ~ "rs66887589"
  )) %>% 
  filter(!is.na(SNP)) %>%
  select(SNP,
         chr           = "chromosome",
         pos           = "base_pair_location",
         effect_allele.outcome = "effect_allele",
         other_allele.outcome  = "other_allele",
         eaf.outcome           = "effect_allele_frequency",
         pval.outcome          = "p_value",
         beta.outcome          = "beta",
         se.outcome            = "standard_error",
         samplesize.outcome    = "N") %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome")

write.table(snps_Wightman_DBP,here('Results','Instruments','iv_Wightman_DBP.txt'))

# Alzheimer disease GWAS - de Rojas
gwasAlzD_deRojas <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_deRojas_2021/Sumstats_SPIGAPUK2_20190625/Sumstats_SPIGAPUK2_20190625.txt')))

snps_deRojas_DBP <- gwasAlzD_deRojas %>%
  filter(RS %in% c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589")) %>%
  select(SNP = RS,
         chr = CHR,
         effect_allele.outcome = A1, 
         other_allele.outcome  = A2,
         beta.outcome          = Beta,
         se.outcome            = SE,
         pval.outcome          = P) %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome",
         eaf.outcome = NA)
write.table(snps_deRojas_DBP,here('Results','Instruments','iv_deRojas_DBP.txt'))

# Alzheimer disease GWAS - Bellenguez
gwasAlz_Bellenguez <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_Bellenguez_2022/GCST90027158_buildGRCh38.tsv.gz')))

snps_Bellenguez_DBP <- gwasAlz_Bellenguez %>%
  filter(variant_id %in% c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589")) %>%
  select(SNP = variant_id,
         chr = chromosome,
         effect_allele.outcome = effect_allele,
         other_allele.outcome  = other_allele,
         beta.outcome          = beta,
         se.outcome            = standard_error,
         eaf.outcome           = effect_allele_frequency,
         pval.outcome          = p_value) %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome")
write.table(snps_Bellenguez_DBP,here('Results','Instruments','iv_Bellenguez_DBP.txt'))

# Sistolic blood pressure ======================================================
gwasSBP <- as_tibble(read_table(paste0(pathToData,'GWAS/BP_Evangelou_2018/Evangelou_30224653_SBP.txt.gz'))) %>%
  separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
  mutate(pos = as.numeric(pos))

snpsSBP <- gwasSBP %>%
  filter(chr == 4  & Type =='SNP') %>%
  mutate(SNP = case_when(
    pos == 120423094 ~ "rs80223330",
    pos == 120502461 ~ "rs12646525",
    pos == 120416096 ~ "rs17355550",
    pos == 120544112 ~ "rs7672519"
  )) %>%
  select(SNP, chr, pos,
         effect_allele.exposure = "Allele1",
         other_allele.exposure  = "Allele2",
         pval.exposure          = "P",
         samplesize.exposure    = "N_effective",
         beta.exposure          = "Effect",
         se.exposure            = "StdErr",
         eaf.exposure           = 'Freq1')
write.table(snpsSBP,here('Results','Instruments','iv_SBP.txt'))

# Alzheimer disease gwas - Lambert 
snps_lambert_SBP <- gwasAlzD_lambert %>%
  filter(MarkerName %in% snpsSBP$SNP) %>%
  select(SNP           = "MarkerName",
         chr           = "Chromosome",
         pos           = "Position", #hg19 
         effect_allele.outcome = "Effect_allele",
         other_allele.outcome  = "Non_Effect_allele",
         pval.outcome          = "Pvalue",
         beta.outcome          = "Beta",
         se.outcome            = "SE") %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome")
write.table(snps_lambert_SBP,here('Results','Instruments','iv_Lambert_SBP.txt'))

# Alzheimer disease gwas - Wightman
snps_Wightman_SBP <- gwasAlzD_Wightman %>%
  filter(chromosome == 4) %>%
  mutate(SNP = case_when(
    base_pair_location == 120423094 ~ "rs80223330",
    base_pair_location == 120502461 ~ "rs12646525",
    base_pair_location == 120416096 ~ "rs17355550",
    base_pair_location == 120544112 ~ "rs7672519"
  )) %>% 
  filter(!is.na(SNP)) %>%
  select(SNP,
         chr           = "chromosome",
         pos           = "base_pair_location",
         effect_allele.outcome = "effect_allele",
         other_allele.outcome  = "other_allele",
         eaf.outcome           = "effect_allele_frequency",
         pval.outcome          = "p_value",
         beta.outcome          = "beta",
         se.outcome            = "standard_error",
         samplesize.outcome    = "N") %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome")

write.table(snps_Wightman_SBP,here('Results','Instruments','iv_Wightman_SBP.txt'))

# Alzheimer disease gwas - deRojas
snps_deRojas_SBP <- gwasAlzD_deRojas %>%
  filter(RS %in% snpsSBP$SNP) %>%
  select(SNP = RS,
         chr = CHR,
         effect_allele.outcome = A1, 
         other_allele.outcome  = A2,
         beta.outcome          = Beta,
         se.outcome            = SE,
         pval.outcome          = P) %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome",
         eaf.outcome = NA)
write.table(snps_deRojas_SBP,here('Results','Instruments','iv_deRojas_SBP.txt'))

# Alzheimer disease GWAS - Bellenguez
snps_Bellenguez_SBP <- gwasAlz_Bellenguez %>%
  filter(variant_id %in% snpsSBP$SNP) %>%
  select(SNP = variant_id,
         chr = chromosome,
         effect_allele.outcome = effect_allele,
         other_allele.outcome  = other_allele,
         beta.outcome          = beta,
         se.outcome            = standard_error,
         eaf.outcome           = effect_allele_frequency,
         pval.outcome          = p_value) %>%
  mutate(id.outcome = "outcome",
         outcome    = "outcome")
write.table(snps_Bellenguez_SBP,here('Results','Instruments','iv_deRojas_SBP.txt'))




# a <- harmonise_data(iv_DBP,iv_Alz_DBP) %>%
#   summarise(SNP = SNP,
#          CHR = chr.outcome, 
#          POS = pos.outcome,
#          'EA_exposure' = effect_allele.exposure,
#          'EA_outcome'  = effect_allele.outcome,
#          'OA_exposure' = other_allele.exposure,
#          'OA_outcome' = other_allele.outcome,
#          'EAF_exposure' = round(eaf.exposure, digits = 2),
#          'EAF_outcome'  = round(eaf.outcome, digits = 2),
#          'Beta_exposure' = round(beta.exposure, digits = 2),
#          'Beta_outcome' = round(beta.outcome, digits = 2),
#          'SE_exposure' = round(se.exposure, digits = 2),
#          'SE_outcome'  = round(se.outcome, digits = 2),
#          'PVAL_exposure' = formatC(pval.exposure, format = "e", digits = 1),
#          'PVAL_outcome' = pval.outcome) 
# write.table(a, 'DBP_AlzD_harmonised.txt')
# a <- a %>%
#   flextable() %>%
#   span_header(sep = "_") %>%
#   align(i = 1, align = 'center', part = "header") %>%
#   bg(j = "EA_exposure", bg = "#EFEFEF", part = "all") %>%
#   bg(j = "OA_exposure", bg = "#EFEFEF", part = "body") %>%
#   bg(i = 2, j = 6, bg = "#EFEFEF",part = "header") %>%
#   bg(j = "EAF_exposure",bg = "#EFEFEF", part = "all") %>%
#   bg(j = "Beta_exposure", bg = "#EFEFEF", part = "body") %>%
#   bg(j = "SE_exposure",bg = "#EFEFEF", part = "all") %>%
#   bg(i = 2, j = 10, bg = "#EFEFEF", part = "header") %>%
#   bg(i = 2, j = 14,  bg = "#EFEFEF", part = "header") %>%
#   bg(j = "PVAL_exposure",bg = "#EFEFEF", part = "body") %>%
#   vline(i = 1, j = c(3:13), part = "header") %>%
#   vline(i = 2, j = seq(3,13,2), part = "header") %>%
#   vline(j = seq(3,14,2), part = "body") %>%
#   set_caption("Exposure = DBP, Outcome = AlzD (Lambert2013)")
# # save_as_docx(a, path = here('iv_DBP.docx'))
# 
# 
# # SBP --------------------------------------------------------------------------
# iv_SBP <- as_tibble(read.table('iv_SBP.txt')) %>%
#   format_data(type = "exposure") %>%
#   mutate(id.exposure = "exposure",
#          exposure = "exposure")
# 
# iv_Alz_SBP <- as_tibble(read.table('iv_AlzD_SBP.txt')) %>% 
#   format_data(type = "outcome") %>%
#   mutate(id.outcome = "outcome",
#          outcome = "outcome")
# 
# b <- harmonise_data(iv_SBP,iv_Alz_SBP) %>%
#   summarise(SNP = SNP,
#             CHR = chr.outcome, 
#             POS = pos.outcome,
#             'EA_exposure' = effect_allele.exposure,
#             'EA_outcome'  = effect_allele.outcome,
#             'OA_exposure' = other_allele.exposure,
#             'OA_outcome' = other_allele.outcome,
#             'EAF_exposure' = round(eaf.exposure, digits = 2),
#             'EAF_outcome'  = round(eaf.outcome, digits = 2),
#             'Beta_exposure' = round(beta.exposure, digits = 2),
#             'Beta_outcome' = round(beta.outcome, digits = 2),
#             'SE_exposure' = round(se.exposure, digits = 2),
#             'SE_outcome'  = round(se.outcome, digits = 2),
#             'PVAL_exposure' = formatC(pval.exposure, format = "e", digits = 1),
#             'PVAL_outcome' = pval.outcome) 
# write.table(b, 'SBP_AlzD_harmonised.txt')
# b <- b %>%
#   flextable() %>%
#   span_header(sep = "_") %>%
#   align(i = 1, align = 'center', part = "header") %>%
#   bg(j = "EA_exposure", bg = "#EFEFEF", part = "all") %>%
#   bg(j = "OA_exposure", bg = "#EFEFEF", part = "body") %>%
#   bg(i = 2, j = 6, bg = "#EFEFEF",part = "header") %>%
#   bg(j = "EAF_exposure",bg = "#EFEFEF", part = "all") %>%
#   bg(j = "Beta_exposure", bg = "#EFEFEF", part = "body") %>%
#   bg(j = "SE_exposure",bg = "#EFEFEF", part = "all") %>%
#   bg(i = 2, j = 10, bg = "#EFEFEF", part = "header") %>%
#   bg(i = 2, j = 14,  bg = "#EFEFEF", part = "header") %>%
#   bg(j = "PVAL_exposure",bg = "#EFEFEF", part = "body") %>%
#   vline(i = 1, j = c(3:13), part = "header") %>%
#   vline(i = 2, j = seq(3,13,2), part = "header") %>%
#   vline(j = seq(3,14,2), part = "body") %>%
#   set_caption("Exposure = SBP, Outcome = AlzD (Lambert2013)")
# # save_as_image(a, here('iv_SBP.png'), expand = 10, res = 300)
# 
