install.packages("devtools") 
library(devtools) 
install_github("oliviasabik/RACER") 
library(RACER)
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
  dplyr::filter(chr == 4) |>
  dplyr::filter(pos >= (120415558-150000) & pos <= (120549959+150000)) |>
  dplyr::mutate(
    effect_allele.exposure = toupper(effect_allele.exposure),
    other_allele.exposure  = toupper(other_allele.exposure)
  )


# submit the query
gwas_lambert <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/AlzD_Lambert_IGAP_2013/IGAP_summary_statistics/IGAP_stage_1.txt'))) %>%
  dplyr::select(SNP           = "MarkerName",
                chr           = "Chromosome",
                pos           = "Position", #hg19 
                effect_allele.outcome = "Effect_allele",
                other_allele.outcome  = "Non_Effect_allele",
                pval.outcome          = "Pvalue",
                beta.outcome          = "Beta",
                se.outcome            = "SE") %>%
  dplyr::mutate(id.outcome = "Lambert",
                outcome    = "Lambert",
                eaf.outcome = NA,
                chr = as.numeric(chr),
                samplesize.outcome = 54162)

gwas <- gwas |> inner_join(gwas_lambert |> dplyr::select(chr, pos, SNP))

# Erectile disfunction
ed <- read.delim("D:/Projects/PDE5_AD_MendelianRandomisation/SummaryStatistics/summary_stats_finngen_R10_ERECTILE_DYSFUNCTION/summary_stats_finngen_R10_ERECTILE_DYSFUNCTION")
ed <- ed |> filter(X.chrom == 4) 
ed <- ed |> rename(chr = X.chrom, SNP = rsids) |> inner_join(gwas |> dplyr::select(SNP))

# Alzheimer's disease
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
  

gwas_racer <- RACER::formatRACER(assoc_data = gwas, chr_col = 1, pos_col = 2, p_col = 5) |> as_tibble()
ed        <- RACER::formatRACER(assoc_data = ed, chr_col = 1, pos_col = 2, p_col = 7) |> as_tibble()
alz_racer <- RACER::formatRACER(assoc_data = alz|> left_join(gwas_lambert |> dplyr::select(chr, pos, SNP)), chr_col = 1, pos_col = 2, p_col = 6) |> as_tibble()


gwas_racer1 <- RACER::ldRACER(assoc_data = gwas_racer , rs_col = 12, pops = "EUR", lead_snp = "rs2389873") 
ed_racer1 <- RACER::ldRACER(assoc_data = ed, rs_col = 5, pops = "EUR", lead_snp = "rs2389873") 
alz_racer1 <- RACER::ldRACER(assoc_data = alz_racer, rs_col = 5, pops = "EUR", lead_snp = "rs2389873") 

ed_racer2 <- ed_racer1 |> dplyr::select(-POS) |> inner_join(gwas_racer1 |> dplyr::select(RS_ID, POS))
ed_racer2 <- ed_racer2 |> relocate(POS, .after = RS_ID)

RACER::singlePlotRACER(assoc_data = alz_racer1 |> dplyr::select(-"LABEL"), chr = 4, build = "hg19", plotby = "gene", gene_plot = "PDE5A")
RACER::singlePlotRACER(assoc_data = ed_racer2  |> dplyr::select(-"LABEL"), chr = 4, build = "hg19", plotby = "gene", gene_plot = "PDE5A")

mirrorPlotRACER(assoc_data1 = gwas_racer1, assoc_data2 = ed_racer2, chr = 4, plotby = "gene", gene_plot = "PDE5A")
mirrorPlotRACER(assoc_data1 = gwas_racer1, assoc_data2 = alz_racer1, chr = 4, plotby = "gene", gene_plot = "PDE5A")




data("mark3_bmd_gwas")
data("mark3_eqtl")