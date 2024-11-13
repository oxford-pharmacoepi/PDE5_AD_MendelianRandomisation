# Colocalisation
source(here::here("Functions/Functions.R"))

chr <- 4
window <- 0
gene_start <- 120415550
gene_end   <- 120550146 # from ensembl

exp <- loadGwas("DBP", onlyInstruments = FALSE) %>%
  dplyr::filter(chr == 4, pos >= gene_start-window, pos <= gene_end+window)
out <- as_tibble(read.delim(paste0(pathResults, "UK Biobank/Colocalisation_outcome.tsv"))) |>
  rename("SNP" = "snp") |>
  left_join(loadGwas("Lambert", onlyInstruments = FALSE) |> select("SNP", "chr", "pos"), by = "SNP")

# harmonise data
dat <- exp %>% 
  dplyr::inner_join(out, by = c("chr", "pos")) %>%
  dplyr::mutate(beta.outcome = dplyr::if_else(effect_allele.exposure == effect_allele.outcome, beta.outcome, -beta.outcome),
                eaf.outcome  = dplyr::if_else(effect_allele.exposure == effect_allele.outcome, eaf.outcome, 1-eaf.outcome),
                other_allele.outcome  = other_allele.exposure,
                effect_allele.outcome = effect_allele.exposure) %>%
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

coloc_table <- coloc_tab  

readr::write_delim(coloc_table, paste0(pathResults,"Colocalisation/colocalisation.txt"))


# Create a plot
table <- tibble(
  pval = -log10(exp_list$pvalues),
  pos = exp_list$pos,
  type = "Diastolic blood pressure") |>
  rbind(
    tibble(
      pval = -log10(out_list$pvalues),
      pos = out_list$pos,
      type = "Alzheimer's disease"
    )
  )

p1 <- ggplot(table, aes(x = pos, y = pval)) +
  geom_point(colour = "#6A99D0") +
  facet_grid(rows = vars(type)) +
  xlab("Position (Chromosome 4)") +
  ylab(bquote(-log[10](PValue))) +
  theme_bw()

ggsave(p1, filename = paste0(pathResults,"Figures/Colocalisation_UKB.png"),width = 20, height = 15, units = "cm", dpi = 600)
