tables <- list()

# Table 1: Instrumental variables ----
tables[["t1"]] <- read.table(paste0(pathResults,"MR_Results/harmonised_DBP.txt"), header = TRUE) |>
  select("SNP", ends_with("exposure")) |> distinct() |> arrange(desc(pval.exposure)) |>
  select("SNP","Effect allele" = "effect_allele.exposure", "Other allele" = "other_allele.exposure",
         "EAF" = "eaf.exposure", "Beta" = "beta.exposure", "SE" = "se.exposure", "P-Value" = "pval.exposure", "Sample size" = "samplesize.exposure"
         )  |>
  arrange(`P-Value`) |>
  mutate(Beta = round(Beta,2),
         SE = round(SE, 2),
         EAF = round(EAF,2)) |>
  mutate("P-Value" = formatC(`P-Value`,digits = 0)) |>
  flextable() |>
  bold(part = "header") |>
  align(j = c(2,3,4,5,6,7,8), align = "center")

# ST2. Linkage disequilibrium ----
# Linkage disequilibrium matrix of the instruments used for the Mendelian randomization analysis. Notice that values correspond to r (not r2).
iv <- read.table(paste0(pathResults,"InstrumentSelection/ld_matrix.txt")) %>%
  as.data.frame() %>%
  dplyr::rename("rs10050092 (C/T)" = "rs10050092_C_T",
                "rs12646525 (T/C)" = "rs12646525_T_C",
                "rs17355550 (C/T)" = "rs17355550_C_T",
                "rs66887589 (C/T)" = "rs66887589_C_T",
                "rs80223330 (A/G)" = "rs80223330_A_G") %>%
  dplyr::mutate(`rs10050092 (C/T)` = round(`rs10050092 (C/T)`, digits = 2),
                `rs12646525 (T/C)` = round(`rs12646525 (T/C)`, digits = 2),
                `rs17355550 (C/T)` = round(`rs17355550 (C/T)`, digits = 2),
                `rs66887589 (C/T)` = round(`rs66887589 (C/T)`, digits = 2),
                `rs80223330 (A/G)` = round(`rs80223330 (A/G)`, digits = 2)) 
tables[["ST2"]] <- iv %>%
  dplyr::mutate(SNP = colnames(iv)) %>%
  dplyr::relocate(SNP) %>%
  flextable::flextable() %>%
  flextable::bold(i = 1, bold = TRUE, part = "header") %>%
  flextable::bold(j = 1, bold = TRUE, part = "body") %>%
  flextable::align(part = "all", align = "center") %>%
  flextable::vline(j = 1)



# ST3. Harmonised variants ----
# Harmonized variants respect the alleles of the linkage disequilibrium matrix. 
# *Note:* EA = Effect allele, OA = Other allele, SE = Standard error, EAF = Effect allele frequency.
tab <- read.table(paste0(pathResults,"MR_Results/harmonised_DBP.txt"), header = TRUE) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    beta.exposure = round(beta.exposure, digits = 2),
    beta.outcome  = round(beta.outcome, digits = 2),
    se.exposure   = round(se.exposure, digits = 2),
    se.outcome    = round(se.outcome, digits = 2),
    pval.exposure = formatC(pval.exposure, digits = 1, format = "e"),
    pval.outcome  = formatC(pval.outcome, digits = 1, format = "e"),
    eaf.exposure  = round(eaf.exposure, digits = 2),
    eaf.outcome   = round(eaf.outcome, digits = 2)
  ) %>%
  mutate(outcome = case_when(
    outcome == "Lambert" ~"Lambert et al. (2013)",
    outcome == "Wightman" ~"Wightman et al. (2021)", 
    outcome == "deRojas" ~ "De Rojas et al. (2021)",
    outcome == "Bellenguez" ~"Bellenguez et al. (2022)"))   |>
  rename("Study outcome" = "outcome") |>
  dplyr::select("Study outcome", "SNP","EA" = "effect_allele", "OA" = "other_allele",
                "EAF_Exposure"  = "eaf.exposure",
                "EAF_Outcome"   = "eaf.outcome",
                "Beta_Exposure" = "beta.exposure",
                "Beta_Outcome"  = "beta.outcome",
                "SE_Exposure"   = "se.exposure",
                "SE_Outcome"    = "se.outcome",
                "PValue_Exposure" = "pval.exposure",
                "PValue_Outcome"  = "pval.outcome")
tables[["ST3"]]  <- tab %>%
  flextable::flextable() %>%
  ftExtra::span_header(sep = "_") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::bold(i = c(1,2), bold = TRUE, part = "header") %>%
  flextable::vline(j = seq(4,11,2)) |>
  flextable::merge_at(j = 1, c(1:5)) |>
  flextable::merge_at(j = 1, c(6:10)) |>
  flextable::merge_at(j = 1, c(11:15)) |>
  flextable::merge_at(j = 1, c(16:20)) |>
  flextable::hline(i = c(5,10,15,20)) |>
  flextable::valign(j = 1, valign = "top") |>
  flextable::align(j = 1, align = "left") |>
  flextable::width(j = 1, width = 8, unit = "cm")



# ST4. Baseline characteristics
ukb_data <- loadUkbData(pathUKB) |>
  mutate(index_of_multiple_deprivation = if_else(is.na(index_of_multiple_deprivation_england),
                                                 index_of_multiple_deprivation_wales,
                                                 index_of_multiple_deprivation_england),
         index_of_multiple_deprivation = if_else(is.na(index_of_multiple_deprivation),
                                                 index_of_multiple_deprivation_scotland,
                                                index_of_multiple_deprivation)) |>
  mutate(ethnic_background = case_when(
    ethnic_background %in% c(1,1001,2001,3001,4001) ~ "White",
    ethnic_background %in% c(-1,-3) ~ "Unknown",
    .default = "Non-white"
  )) |>
  mutate(ethnic_background = as.factor(ethnic_background)) 
  
males <- ukb_data |> filter(sex == 1)
females <- ukb_data |> filter(sex == 0)

t1_males <- tableone::CreateTableOne(data = males, 
                         strata = "ad_status", 
                         vars = c("age_when_assessment",
                                  "index_of_multiple_deprivation",
                                  "ethnic_background",
                                  "body_mass_index"),
                         test = FALSE
                         ) |> 
  print(smd = TRUE) |> 
  as.data.frame() 
tables[["ST4_males"]] <- t1_males |>
  mutate("Variable" = row.names(t1_males)) |>
  relocate(Variable) |>
  flextable() 

t1_females <- tableone::CreateTableOne(data = females, 
                                     strata = "ad_status", 
                                     vars = c("age_when_assessment",
                                              "index_of_multiple_deprivation",
                                              "ethnic_background",
                                              "body_mass_index"),
                                     test = FALSE
) |> 
  print(smd = TRUE)  |> 
  as.data.frame() 
tables[["ST5_females"]] <- t1_females |>
  mutate("Variable" = row.names(t1_females)) |>
  relocate(Variable) |>
  flextable() 


# ST5. Colocalisation results ----
tables[["ST5"]] <- read.table(paste0(pathResults,"Colocalisation/colocalisation.txt"), header = TRUE) |>
  mutate(outcome = c("Lambert et al. (2013)", "Wightman et al. (2021)", "De Rojas et al. (2021)", "Bellenguez et al. (2022)")) |>
  rename("Outcome study" = "outcome", "N SNPs" = "nsnps") |>
  relocate("Outcome study") |>
  relocate(P1, .after = `N SNPs`) |>
  relocate(P2, .after = P1) |>
  relocate(P12, .after = P2) |>
  mutate_at(vars(starts_with("H")), funs(. * 100)) |>
  mutate_at(vars(starts_with("H")), funs(round(.,digits = 2))) |>
  mutate_at(vars(starts_with("H")), funs(as.character)) |>
  mutate_at(vars(starts_with("H")), funs(if_else(.=="0","<0.01",.))) |>
  select(-"exposure") |>
  rename_at(vars(starts_with("P")), funs(paste0(.," (%)"))) |>
  rename_at(vars(starts_with("H")), funs(paste0(.," (%)"))) |>
  flextable() |>
  ftExtra::span_header(sep = "_") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::bold(i = c(1), bold = TRUE, part = "header") %>%
  flextable::vline(j = c(2,5)) |>
  flextable::width(j = 1, width = 5, unit = "cm")

# ST6. Two-step cis-MR ----
tab_twostep <- read.table(paste0(pathResults,"TwoStepMR/twostepMR.txt"), header = TRUE) |>
  mutate(outcome = case_when(
    outcome == "Lambert" ~"Lambert et al. (2013)",
    outcome == "Wightman" ~"Wightman et al. (2021)", 
    outcome == "deRojas" ~ "De Rojas et al. (2021)",
    outcome == "Bellenguez" ~"Bellenguez et al. (2022)"))   |>
  dplyr::distinct() %>%
  dplyr::relocate(outcome) %>%
  dplyr::mutate("Confounder" = dplyr::case_when(
    confounder == "ukb-b-19953" ~ "BMI",
    confounder == "ukb-b-7376"  ~ "Impedance of leg (right)",
    confounder == "ukb-b-14068" ~ "Impedance of leg (left)",
    confounder == "ukb-b-7859"  ~ "Impedance of arm (right)",
    confounder == "ukb-b-19379" ~ "Impedance of arm (left)",
    confounder == "ukb-b-19921" ~ "Impedance of whole body",
    confounder == "ukb-b-10787" ~ "Standing height",
    confounder == "ebi-a-GCST004607" ~ "Plateletcrit",
    confounder == "ebi-a-GCST004626" ~ "Myeloid white cell count",
    confounder == "ukb-d-30080_irnt" ~ "Platelet count",
    confounder == "ieu-b-30" ~ "White blood cell count",
    confounder == "ebi-a-GCST005195" ~ "Coronary artery disease",
    confounder == "ebi-a-GCST004614" ~ "Granulocyte count",
    confounder == "ebi-a-GCST004620" ~ "Sum basophil neutrophil counts"
  )) %>%
  dplyr::arrange(Confounder) %>%
  dplyr::rename("GWAS ID" = "confounder",
                "Study outcome" = "outcome") %>%
  dplyr::mutate("Associated SNPs" = dplyr::case_when(
    Confounder == "BMI" ~ " ",
    Confounder == "Impedance of leg (right)"       ~ "rs10050092, rs12646525",
    Confounder == "Impedance of leg (left)"        ~ "rs10050092, rs12646525",
    Confounder == "Impedance of arm (right)"       ~ "rs10050092",
    Confounder == "Impedance of arm (left)"        ~ "rs10050092",
    Confounder == "Impedance of whole body"        ~ "rs10050092",
    Confounder == "Standing height"                ~ "rs10050092",
    Confounder == "Plateletcrit"                   ~ "rs66887589",
    Confounder == "Myeloid white cell count"       ~ "rs10050092, rs66887589",
    Confounder == "Platelet count"                 ~ "rs10050092",
    Confounder == "White blood cell count"         ~ "rs66887589",
    Confounder == "Coronary artery disease"        ~ "rs66887589, rs80223330",
    Confounder == "Granulocyte count"              ~ "rs66887589",
    Confounder == "Sum basophil neutrophil counts" ~ "rs66887589")
  ) %>%
  dplyr::relocate("Confounder", .after = "Study outcome") %>%
  dplyr::relocate("GWAS ID",    .after = "Confounder") %>%
  dplyr::relocate("Associated SNPs", .after = "GWAS ID") %>%
  dplyr::mutate(OR = round(exp(beta_IVWcorrel), digits = 2)) %>%
  dplyr::mutate(SE = round(se_IVWcorrel.random, digits = 2)) %>%
  dplyr::mutate(`P-Value` = round(p_IVWcorrel.random, digits = 2)) %>%
  dplyr::select(-tidyselect::contains("IVWcorrel"), -n_snp, -F_stat, -exposure)

tables[["ST6"]] <- tab_twostep |>
  flextable::flextable() %>%
  flextable::bold(bold = TRUE, part = "header") %>%
  flextable::align(align = c("center"), part = "all") |>
  flextable::width(j = 1, width = 5, unit = "cm") |>
  flextable::align(align = "left", j = 1) |>
  flextable::hline(i = seq(4,nrow(tab_twostep),4))

# ST7. Leave one out analysis ----
tables[["ST7"]] <- tibble(
  "outcome" = c("Lambert","Lambert", "Lambert", "Lambert", "Lambert", "Lambert", 
                "Wightman", "Wightman","Wightman","Wightman","Wightman","Wightman",
                "deRojas", "deRojas", "deRojas", "deRojas", "deRojas", "deRojas", 
                "Bellenguez","Bellenguez","Bellenguez","Bellenguez","Bellenguez","Bellenguez"),
  "SNP" = c("All","rs10050092","rs66887589", "rs12646525","rs80223330","rs17355550",
            "All","rs10050092","rs66887589", "rs12646525","rs80223330","rs17355550",
            "All","rs10050092","rs66887589", "rs12646525","rs80223330","rs17355550",
            "All","rs10050092","rs66887589", "rs12646525","rs80223330","rs17355550")
) |>
  left_join(
    read.table(paste0(pathResults,"LeaveOneOutAnalysis/LeaveOneOut_Results_DBP.txt"),
               header = TRUE) |>
      select( "outcome", "SNP", "N Instruments" = "instruments", "OR", "SE" = "se","P-value" = "pval") |>
      mutate_at(vars(c("SE","P-value")), ~round(.,digits = 2)) |>
      mutate_at(vars(c("OR","P-value")), ~round(.,digits = 3)),
    by = c("outcome","SNP")
  ) |>
  mutate(outcome = case_when(
    outcome == "Lambert" ~"Lambert et al. (2013)",
    outcome == "Wightman" ~"Wightman et al. (2021)", 
    outcome == "deRojas" ~ "De Rojas et al. (2021)",
    outcome == "Bellenguez" ~"Bellenguez et al. (2022)"))   |>
  rename("Study outcome" = "outcome") |>
  flextable() |>
  ftExtra::span_header(sep = "_") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::bold(i = c(1), bold = TRUE, part = "header") |>
  flextable::width(j = 1, width = 5, unit = "cm") |>
  flextable::align(align = "left", j = 1) |>
  flextable::hline(i = seq(6,6*4,6))

save_as_docx(values = tables, path = "C:/Users/martaa/OneDrive - Nexus365/1_PDE5_AlzD/ManuscriptVersions/1-FirstVersion/Tables.docx")
save_as_docx(values = tables, path = paste0(pathResults, "Tables.docx"))
