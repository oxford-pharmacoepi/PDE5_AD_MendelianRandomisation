windowsFonts("Calibri" = windowsFont("Calibri"))
library(forestplot)
tables <- list()

# Main figure ----
res <- read.table(paste0(pathResults,"MR_Results/MR_Results_DBP_scaled.txt"), header = TRUE) |>
  filter(outcome != "Kunkle") |>
  mutate(samplesize = c("54,162","398,108", "409,435", "487,511"),
         cases = c("17,008 (31)", "39,968 (10)", "81,611 (20)", "85,934 (18)"),
         outcome = c("Lambert et al. (2013)", "Wightman et al. (2021)", "De Rojas et al. (2021)", "Bellenguez et al. (2022)")) |>
  add_row(outcome = "Main analysis",.before = 1) |>
  add_row(outcome = "Secondary analyses", .before = 3) |>
  mutate(orci = ifelse(row_number() %in% c(1,3)," ",paste0(round(OR,2)," (",round(cilow,2),", ",round(cihigh,2),")")),
         pval = round(pval, 2)) |>
  rename(mean = OR, lower = cilow, upper = cihigh) 

res |>
  forestplot(labeltext = c(outcome, samplesize, cases, orci, pval),
             clip = c(0.9,1.1),
             xlog = TRUE,
             boxsize = 0.1,
             line.margin = 0.22,
             xlab = "Odds Ratio",
             xticks = log(c(0.9, 1, 1.1)),
             graphwidth = unit(7, "cm"),
             graph.pos = 4,
             is.summary = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)) |>
  fp_add_lines(h_2 = gpar(lty = 1)) |>
  fp_add_header(outcome = "Outcome study",
                samplesize = "Sample size",
                cases = "Cases (%)",
                orci = "OR (95% CI)",
                pval = "P-Value") |>
  # fp_set_zebra_style("white","#F1F5F5","white","#F1F5F5","white") |>
  fp_set_style(txt_gp = fpTxtGp(label = gpar(fontfamily = "Calibri", cex = 1),
                                ticks = gpar(cex = 0.65),
                                xlab = gpar(cex = 0.75),
                                legend = gpar(cex = 1),
                                title = gpar(fontfamily = "Calibri"),
                                summary = gpar(fontface = "bold")),
               box = "#4472C4", lines = "#4472C4")



# Sex specific analysis ----
res <- 
  read.csv(paste0(pathResults,"SexSpecific/MR_Result.txt")) |>
  mutate(variable = stringr::str_to_sentence(outcome)) |>
  mutate(outcome = if_else(row_number() == 1, "UK Biobank", "")) |>
  mutate(samplesize = c("220,352","262,037"),
         cases = c("1,722 (0.78)","2,030 (0.77)")) |>
  mutate(orci = paste0(round(OR,2)," (",round(cilow,2),", ", round(cihigh,2),")"),
         pval = round(pval, 2)) |>
  rename(mean = OR, lower = cilow, upper = cihigh)

res |>
  forestplot(labeltext = c(variable, samplesize, cases, orci, pval),
             clip = c(0.8,1.2),
             xlog = TRUE,
             boxsize = 0.1,
             line.margin = 0.22,
             xlab = "Odds Ratio",
             xticks = c(.8,.9,1, 1.1, 1.2),
             graphwidth = unit(7, "cm"),
             graph.pos = 4) |>
  fp_add_lines(h_2 = gpar(lty = 1)) |>
  fp_add_header(variable = "Sex", samplesize = "Sample size",
                cases = "Cases (%)",
                orci = "OR (95% CI)",
                pval = "P-Value") |>
  forestplot::fp_set_zebra_style("#F1F5F5","white") |>
  fp_set_style(txt_gp = fpTxtGp(label = gpar(fontfamily = "Calibri", cex = 1),
                                ticks = gpar(cex = 0.65),
                                xlab = gpar(cex = 0.75),
                                legend = gpar(cex = 1),
                                title = gpar(fontfamily = "Calibri"),
                                summary = gpar(fontface = "bold")),
               box = "#FEAE00", lines = "#FEAE00")
