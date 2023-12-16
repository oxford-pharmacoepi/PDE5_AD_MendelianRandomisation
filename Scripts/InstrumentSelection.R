
source(here::here("Functions/loadGwas.R"))
dir.create(paste0(pathResults,"InstrumentSelection/"))

# Diastolic blood pressure -----------------------------------------------------
# Snp position: hg19
snps <- loadGwas("DBP", onlyInstruments = TRUE)
write.table(snps, paste0(pathResults,"InstrumentSelection/iv_DBP.txt"))

# LD matrix --------------------------------------------------------------------
ld_mat <- TwoSampleMR::ld_matrix(snps$SNP)
write.table(ld_mat,paste0(pathResults,"InstrumentSelection/ld_matrix_DBP.txt"))

# Systolic blood pressure  -----------------------------------------------------
# Snp position: hg19
snps <- loadGwas("SBP", onlyInstruments = TRUE)
write.table(snps, paste0(pathResults,"InstrumentSelection/iv_SBP.txt"))

# LD matrix --------------------------------------------------------------------
ld_mat <- TwoSampleMR::ld_matrix(snps$SNP)
write.table(ld_mat,paste0(pathResults,"InstrumentSelection/ld_matrix_SBP.txt"))

# Instruments from outcome -----------------------------------------------------
outcome_i  <- c("Lambert","deRojas","Wightman","Bellenguez")

for(outcome in outcome_i){
  gwas <- loadGwas(outcome, onlyInstruments = TRUE)
  print(gwas)
  readr::write_delim(gwas, paste0(pathResults,"InstrumentSelection/iv_",outcome,".txt"))
}