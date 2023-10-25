rm(list = ls())
pacman::p_load('readr','devtools','dplyr','TwoSampleMR','LDlinkR','forestplot','tidyr','coloc',
               'here','ggplot2')
devtools::install_github("bar-woolf/TwoStepCisMR")
library(TwoStepCisMR)
pathToData <- 'C:/Users/marta/Desktop/PDE5_AD_MendelianRandomisation/'

# Functions --------------------------------------------------------------------
mh <- function(outcome, dbp, sbp, col = "black"){
  p1 <- ggplot(outcome,aes(x = pos, y = -log10(pval))) +
    geom_point(show.legend = FALSE, colour = col) +
    labs(x = "", y = bquote(-log[10](italic(p))),
         title = "Alzheimer disease") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(120420000, 120450000, 120480000, 120520000,120550000),
                       labels = c(120420, 120450, 120480,120520,120550)) +
    labs(x = "Chromosome 4 Position (Kb)", y = bquote(-log[10](italic(p))))
  
  p3 <-ggplot(dbp,aes(x = pos, y = -log10(pval))) +
    geom_point(show.legend = FALSE) +
    labs(x = "", y = bquote(-log[10](italic(p))),
         title = "Diastolic blood pressure") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(120420000, 120450000, 120480000, 120520000,120550000),
                       labels = c(120420, 120450, 120480,120520,120550)) +
    labs(x = "Chromosome 4 Position (Kb)", y = bquote(-log[10](italic(p))))
  
  p4 <-ggplot(sbp,aes(x = pos, y = -log10(pval))) +
    geom_point(show.legend = FALSE) +
    labs(x = "", y = bquote(-log[10](italic(p))),
         title = "Systolic blood pressure") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(120420000, 120450000, 120480000, 120520000,120550000),
                       labels = c(120420, 120450, 120480,120520,120550)) + 
    labs(x = "Chromosome 4 Position (Kb)", y = bquote(-log[10](italic(p))))
  
  return(list(p1,p3,p4))
}


createList <- function(d,tip){
  x <- list(
    beta = d$beta,
    MAF  = d$eaf,
    pvalues = d$pval,
    varbeta = d$se^2,
    N = d$n,
    type = tip,
    pos = d$pos,
    chr = d$chr,
    snp = d$snp
  )
}



# Colocalization ---------------------------------------------------------------
chr <- 4
window <- 0 #using a more stringent window for the anlysis. 
gene_start <- 120415550#119494395#
gene_end <- 120550146#119628991#
chrpos <- paste0(chr, ":", gene_start - window, "-", gene_end + window)

# Colocalisation datasets for diastolic blood pressure and systolic blood pressure

# dbp <- as_tibble(read_table(paste0(pathToData,'GWAS/BP_Evangelou_2018/Evangelou_30224653_DBP.txt.gz'))) %>%
#   separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
#   mutate(pos = as.numeric(pos)) %>%
#   mutate(chr = as.numeric(chr)) %>%
#   select(chr,pos,eaf = Freq1, beta = Effect, se = StdErr, pval = P, n = TotalSampleSize)
# write.table(dbp,'dbp_colocalisation.txt')
# sbp <- as_tibble(read_table(paste0(pathToData,'GWAS/BP_Evangelou_2018/Evangelou_30224653_SBP.txt.gz'))) %>%
#   separate(MarkerName, into=c('chr','pos','Type'), sep =':') %>%
#   mutate(pos = as.numeric(pos)) %>%
#   mutate(chr = as.numeric(chr)) %>%
#   select(chr,pos,eaf = Freq1, beta = Effect, se = StdErr, pval = P, n = TotalSampleSize)
# write.table(sbp,'sbp_colocalisation.txt')
dbp <- read.table('dbp_colocalisation.txt')
sbp <- read.table('sbp_colocalisation.txt') 

# WIGHTMAN ---------------------------------------------------------------------
# alz <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_wightman_2021_excluding_23andme/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz')))
# alz <- alz %>% filter(chromosome == 4) %>%
#    select(chr = chromosome, pos = base_pair_location, beta, se = standard_error, eaf = effect_allele_frequency,
#    pval = p_value, n = N)
# write.table(alz,'wightman_colocalisation.txt')
alz <- read.table('wightman_colocalisation.txt') 

dbp <-  dbp %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)
sbp <- sbp %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)
alz <- alz %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)

snp_col <- dbp %>% select(chr, pos) %>%
  inner_join(sbp %>% select(chr,pos)) %>%
  inner_join(alz %>% select(chr,pos)) %>%
  mutate(snp = paste0(chr,':',pos)) %>%
  unique()

dbp_list <- createList(dbp %>% right_join(snp_col), tip = 'quant')
sbp_list <- createList(sbp %>% right_join(snp_col), tip = 'quant')
alz_list <- createList(alz %>% right_join(snp_col) %>% filter(pos != 120455648 & n != 74004), tip = 'cc')

res3 <- coloc::coloc.abf(dbp_list,alz_list)
res4 <- coloc::coloc.abf(sbp_list,alz_list)

p <- mh(alz, dbp, sbp, col = "#EFBA00")
ggsave(here('Figures','mh_Wightman.png'),p[[1]],width = 15, height = 7, units = 'cm', dpi = 300)
ggsave(here('Figures','mh_Wightman_dbp.png'),p[[2]],width = 15, height = 7, units = 'cm', dpi = 300)
ggsave(here('Figures','mh_Wightman_sbp.png'),p[[3]],width = 15, height = 7, units = 'cm', dpi = 300)

# DE ROJAS  --------------------------------------------------------------------
# alz <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_deRojas_2021/Sumstats_SPIGAPUK2_20190625/Sumstats_SPIGAPUK2_20190625.txt')))
# alz <- alz %>% filter(CHR == 4) %>%
#   select(snp = RS, chr = CHR, pos = BP, beta = Beta, se = SE, pval = P)
# write.table(alz,'deRojas_colocalisation.txt')
dbp <- read.table('dbp_colocalisation.txt')
sbp <- read.table('sbp_colocalisation.txt') 
alz <- as_tibble(read.table('deRojas_colocalisation.txt'))

dbp <-  dbp %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)
sbp <- sbp %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)
alz <- alz %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)

snp_col <- dbp %>% select(chr, pos) %>%
  inner_join(sbp %>% select(chr,pos)) %>%
  inner_join(alz %>% select(chr,pos,snp)) %>%
  unique()

dbp_list <- createList(dbp %>% right_join(snp_col), tip = 'quant')
sbp_list <- createList(sbp %>% right_join(snp_col), tip = 'quant')
alz_list <- createList(alz %>% right_join(snp_col), tip = 'cc')
alz_list <- alz_list[names(alz_list) != "MAF"]    

res3 <- coloc::coloc.abf(dbp_list,alz_list)
res4 <- coloc::coloc.abf(sbp_list,alz_list)

p <- mh(alz, dbp, sbp, "#D57100")
ggsave(here('Figures','mh_deRojas.png'),p[[1]],width = 15, height = 7, units = 'cm', dpi = 300)
ggsave(here('Figures','mh_deRojas_dbp.png'),p[[2]],width = 15, height = 7, units = 'cm', dpi = 300)
ggsave(here('Figures','mh_deRojas_sbp.png'),p[[3]],width = 15, height = 7, units = 'cm', dpi = 300)

# Lambert ----------------------------------------------------------------------
# alz <- as_tibble(read_table(paste0(pathToData,'GWAS/AlzD_Lambert_IGAP_2013/IGAP_summary_statistics/IGAP_stage_1.txt')))
# alz <- alz %>% filter(Chromosome == 4) %>%
#   select(snp = MarkerName, chr = Chromosome, pos = Position, beta = Beta, se = SE, pval = Pvalue)
# write.table(alz,'Lambert_colocalisation.txt')

dbp <- read.table('dbp_colocalisation.txt')
sbp <- read.table('sbp_colocalisation.txt') 
alz <- read.table('Lambert_colocalisation.txt')

dbp <-  dbp %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)
sbp <- sbp %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)
alz <- alz %>%
  filter(chr == 4 & pos >= gene_start-window & pos <= gene_end+window)

snp_col <- dbp %>% select(chr, pos) %>%
  inner_join(sbp %>% select(chr,pos)) %>%
  inner_join(alz %>% select(chr,pos,snp)) %>%
  unique()

dbp_list <- createList(dbp %>% right_join(snp_col), tip = 'quant')
sbp_list <- createList(sbp %>% right_join(snp_col), tip = 'quant')
alz_list <- createList(alz %>% right_join(snp_col), tip = 'cc')
alz_list <- alz_list[names(alz_list) != "MAF"]    

res3 <- coloc::coloc.abf(dbp_list,alz_list)
res4 <- coloc::coloc.abf(sbp_list,alz_list)

p <- mh(alz, dbp, sbp, "#71ABB8")
ggsave(here('Figures','mh_Lambert.png'),p[[1]],width = 15, height = 7, units = 'cm', dpi = 300)
ggsave(here('Figures','mh_Lambert_dbp.png'),p[[2]],width = 15, height = 7, units = 'cm', dpi = 300)
ggsave(here('Figures','mh_Lambert_sbp.png'),p[[3]],width = 15, height = 7, units = 'cm', dpi = 300)
