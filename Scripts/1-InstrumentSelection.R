
source(here::here("Functions/loadGwas.R"))
dir.create(paste0(pathResults,"InstrumentSelection/"))

# Diastolic blood pressure -----------------------------------------------------
# Snp position: hg19
snps <- loadGwas("DBP", onlyInstruments = TRUE)
write.table(snps, paste0(pathResults,"InstrumentSelection/iv_DBP.txt"))

# LD matrix --------------------------------------------------------------------
ld_mat <- TwoSampleMR::ld_matrix(snps$SNP)
write.table(ld_mat,paste0(pathResults,"InstrumentSelection/ld_matrix_DBP.txt"))

# Instruments from outcome -----------------------------------------------------
outcome_i  <- c("Lambert","deRojas","Wightman","Bellenguez","Kunkle")
for(outcome in outcome_i){
  gwas <- loadGwas(outcome, onlyInstruments = TRUE)
  print(gwas)
  readr::write_delim(gwas, paste0(pathResults,"InstrumentSelection/iv_",outcome,".txt"))
}

# Check individuals SNPs effect
exp <- read.table(paste0(pathResults,"InstrumentSelection/iv_DBP.txt"))
out <- read.table(paste0(pathResults,"InstrumentSelection/iv_Lambert.txt"), header = TRUE)

p <- exp |>
  inner_join(out, by = c("SNP", "chr","pos")) |>
  as_tibble() |>
  ggplot(aes(x = beta.exposure, y = beta.outcome, label = SNP,
             ymin = beta.outcome - 1.96*se.outcome,
             ymax = beta.outcome + 1.96*se.outcome,
             xmin = beta.exposure - 1.96*se.outcome,
             xmax = beta.exposure + 1.96*se.outcome)) +
  geom_point(size = 2) +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, colour = "darkred") +
  xlim(c(-0.2,0.2)) +
  ylim(c(-0.2,0.2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab("Beta effect of the SNPs on the exposure") +
  ylab("Beta effect of the SNPs on the outcome") +
  theme_bw() +
  geom_text(size = 4, hjust = c(0,0,0,0.5,0)-0.05, vjust = -c(0.5,0.5,0.5,0.5,0.75)) +
  geom_errorbar()
ggsave(paste0(pathResults, "InstrumentSelection/InstrumentsEffect.png"), dpi = 400)  
  
  
