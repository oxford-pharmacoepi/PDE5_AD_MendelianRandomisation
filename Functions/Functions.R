loadGwas <- function(gwas, onlyInstruments = TRUE, bellenguez = TRUE){
  
  gwas <- switch(
    gwas,
    
    "DBP" = {
      gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/BP_Evangelou_2018/Evangelou_30224653_DBP.txt/Evangelou_30224653_DBP.txt'))) %>%
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
          "exposure" = "DBP",
          "id.exposure" = "DBP"
        ) %>%
        dplyr::mutate(
          effect_allele.exposure = toupper(effect_allele.exposure),
          other_allele.exposure  = toupper(other_allele.exposure)
        )
      
      if(onlyInstruments == TRUE){
        gwas <- gwas %>%
          dplyr::inner_join(
            tibble::tibble("pos" = c(120423094,120502461,120416096,120532085,120509279),
                           "SNP" = c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589")),
            by = "pos"
          )
      }
      return(gwas)
    },
    
    
    "SBP" = {
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
      
      if(onlyInstruments == TRUE){
        gwas <- gwas %>%
          dplyr::inner_join(
            tibble::tibble("pos" = c(120423094,120502461,120416096,120544112),
                           "SNP" = c("rs80223330","rs12646525","rs17355550","rs7672519")),
            by = "pos"
          ) 
      }
      return(gwas)
    },
    
    
    "Lambert" = {
      # (Build 37, Assembly Hg19)
      gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/AlzD_Lambert_IGAP_2013/IGAP_summary_statistics/IGAP_stage_1.txt'))) %>%
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
      
      if(onlyInstruments == TRUE){
        gwas <- gwas %>%
          dplyr::inner_join(tibble::tibble(
            "SNP" =  c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589","rs7672519")),
            by = "SNP"
          )
      }
      return(gwas)
    },
    
    "Wightman" = {
      # (Build 37, Assembly Hg19)
      gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/AlzD_wightman_2021_excluding_23andme/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz'))) %>%
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
        
        if(onlyInstruments == TRUE){
          gwas <- gwas %>%
            dplyr::inner_join(tibble::tibble(
              "SNP" = c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589","rs7672519"),
              "pos" = c(120423094,120502461,120416096,120532085,120509279,120544112)),
              by = "pos")
        }
      
      return(gwas)
    },
    
    "deRojas" = {
      gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/AlzD_deRojas_2021/Sumstats_SPIGAPUK2_20190625/Sumstats_SPIGAPUK2_20190625.txt'))) %>%
        dplyr::select(SNP = "RS",
                      chr = "CHR",
                      pos = "BP",
                      effect_allele.outcome = "A1", 
                      other_allele.outcome  = "A2",
                      beta.outcome          = "Beta",
                      se.outcome            = "SE",
                      pval.outcome          = "P") %>%
        dplyr::mutate(id.outcome = "deRojas",
                      outcome    = "deRojas",
                      eaf.outcome = NA,
                      chr = as.numeric(chr),
                      samplesize.outcome = 409435)
      
      if(onlyInstruments == TRUE){
        gwas <- gwas %>%
          dplyr::inner_join(tibble::tibble(
            "SNP" = c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589","rs7672519"),
            "pos" = c(120423094,120502461,120416096,120532085,120509279,120544112)),
            by = c("SNP","pos")
          )
      }
      return(gwas)
    },
    
    "Kunkle" = {
      gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,"SummaryStatistics/AlzD_kunkle_2019/Kunkle_etal_Stage1_results.txt"))) |>
        dplyr::select("SNP" = "MarkerName",
                      "chr" = "Chromosome",
                      "pos" = "Position",
                      "effect_allele.outcome" = "Effect_allele",
                      "other_allele.outcome"  = "Non_Effect_allele",
                      "beta.outcome" = "Beta",
                      "se.outcome" = "SE",
                      "pval.outcome" = "Pvalue") |>
        dplyr::mutate(id.outcome = "Kunkle",
                      outcome = "Kunkle",
                      samplesize.outcome = 63926)
      if(onlyInstruments == TRUE){
        gwas <- gwas %>%
          dplyr::inner_join(tibble::tibble(
            "SNP" = c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589"),
            "pos" = c(120423094,120502461,120416096,120532085,120509279),
            "eaf.outcome" = c(0.9674,0.1411,0.2154,0.5221,0.6616)),
            by = c("SNP","pos")
          )
      }
    },
    
    "Bellenguez" = {
      gwas <- tibble::as_tibble(readr::read_table(paste0(pathData,'SummaryStatistics/AlzD_Bellenguez_2022/GCST90027158_buildGRCh38.tsv.gz'))) %>%
        dplyr::select(SNP = "variant_id",
                      chr = "chromosome",
                      pos = "base_pair_location",
                      effect_allele.outcome = "effect_allele",
                      other_allele.outcome  = "other_allele",
                      beta.outcome          = "beta",
                      se.outcome            = "standard_error",
                      eaf.outcome           = "effect_allele_frequency",
                      pval.outcome          = "p_value") %>%
        dplyr::mutate(id.outcome = "Bellenguez",
                      outcome    = "Bellenguez",
                      chr = as.numeric(chr),
                      samplesize.outcome = 487511)
      
      # Liftover  1200080:1250080
      # Prepare file for liftOver:
    
      if(bellenguez == TRUE){
      gwas <- gwas %>% dplyr::filter(chr == 4)
      
      readr::write_delim(gwas %>%
                           dplyr::select("chr", "pos") %>%
                           dplyr::mutate(pos = as.character(formatC(pos,format = "d"))) %>%
                           dplyr::mutate(chr = paste0("chr",chr,":",pos,"-",pos)) %>%
                           dplyr::select(-"pos"),
                         paste0(pathData,"SummaryStatistics/AlzD_Bellenguez_2022/mapping.txt")) 
      
      liftOver <- readr::read_delim(paste0(pathData,"SummaryStatistics/AlzD_Bellenguez_2022/transformed.bed"), col_names = FALSE, delim = "\t") %>%
        dplyr::mutate(X2 = gsub(x = X1, pattern = ".*:", replacement = ""),
                      X1 = gsub(x = X1, pattern = ":.*", replacement = "")) %>%
        dplyr::mutate(X2 = gsub(x = X2, pattern = "-.*", replacement = ""),
                      X1 = gsub(x = X1, pattern = ".*chr", replacement = "")) %>%
        dplyr::rename("chr" = "X1", "pos" = "X2") %>%
        dplyr::mutate(chr = as.numeric(chr),
                      pos = as.numeric(pos))
      
      gwas1 <- gwas %>%
        dplyr::anti_join(
          readr::read_delim(paste0(pathData,"SummaryStatistics/AlzD_Bellenguez_2022/failed.txt"), col_names = FALSE, delim = "\t") %>%
            dplyr::filter(X1 != "#Deleted in new") %>%
            dplyr::mutate(X2 = gsub(x = X1, pattern = ".*:", replacement = ""),
                          X1 = gsub(x = X1, pattern = ":.*", replacement = "")) %>%
            dplyr::mutate(X2 = gsub(x = X2, pattern = "-.*", replacement = ""),
                          X1 = gsub(x = X1, pattern = ".*chr", replacement = "")) %>%
            dplyr::rename("chr" = "X1", "pos" = "X2") %>%
            dplyr::mutate(chr = as.numeric(chr),
                          pos = as.numeric(pos)),
          by = c("chr", "pos")
        )  %>%
        dplyr::mutate(pos = liftOver$pos)
        
      

      if(onlyInstruments == TRUE){
        gwas1 <- gwas1 %>%
          dplyr::inner_join(tibble::tibble(
            "SNP" = c("rs80223330","rs12646525","rs17355550","rs10050092","rs66887589","rs7672519"),
            "pos" = c(120423094,120502461,120416096,120532085,120509279,120544112)),
            by = c("SNP","pos")
          )
      }
      gwas <- gwas1
      }
      return(gwas)
    }
  )
}

loadUkbData <- function(pathUKB){
  ukb_data <- read.table(paste0(pathUKB,"alzheimer_data.tab"), header=TRUE, sep="\t") |> as_tibble()
  ukb_data$f.42018.0.0 <- as.Date(ukb_data$f.42018.0.0)
  ukb_data$f.42020.0.0 <- as.Date(ukb_data$f.42020.0.0)
  ukb_data$f.42022.0.0 <- as.Date(ukb_data$f.42022.0.0)
  ukb_data$f.42024.0.0 <- as.Date(ukb_data$f.42024.0.0)
  
  
  # Rename
  ukb_data <- ukb_data |>
    select("eid" = "f.eid",
           "sex" = "f.31.0.0",
           "year_of_birth" = "f.34.0.0",
           "age_when_assessment" = "f.21003.0.0",
           "diastolic_blood_pressure" = "f.4079.0.0",
           "date_of_dementia_report"  = "f.42018.0.0",
           "date_of_alzheimer_report" = "f.42020.0.0",
           "date_of_vascular_dementia_report" = "f.42022.0.0",
           "date_of_alzheimer_report" = "f.42020.0.0",
           "date_of_frontotemporal_dementia_report" = "f.42024.0.0",
           "ethnic_background" = "f.21000.0.0",
           "body_mass_index" = "f.21001.0.0",
           "PC1" = "f.22009.0.1", "PC2" = "f.22009.0.2", "PC3" = "f.22009.0.3",
           "PC4" = "f.22009.0.4", "PC5" = "f.22009.0.5", "PC6" = "f.22009.0.6",
           "PC7" = "f.22009.0.7", "PC8" = "f.22009.0.8", "PC9" = "f.22009.0.9", "PC10" = "f.22009.0.10",
           "genetic_ethnic_grouping" = "f.22006.0.0",
           "sex_chromosome_aneuploidy" = "f.22019.0.0",
           "genetic_kinship_to_other_participants" = "f.22021.0.0",
           "heterozygosity" = "f.22027.0.0",
           "index_of_multiple_deprivation_england" = "f.26410.0.0",
           "index_of_multiple_deprivation_wales" = "f.26426.0.0",
           "index_of_multiple_deprivation_scotland" = "f.26427.0.0",
           "age_at_death" = "f.40007.0.0",
           "genetic_batch" = "f.22000.0.0"
    ) 
  
  # Restriction criteria ---
  ukb_data <- ukb_data |> 
    mutate(ad_status = if_else(is.na(date_of_alzheimer_report),0,1)) |>
    # Exclude people with vascular dementia and ad
    mutate(vascular_status = if_else(is.na(date_of_vascular_dementia_report),0,1)) |>
    filter(vascular_status == 0)  |>
    # Exclude people with frontotemporal dementia and ad
    mutate(frontotemporal_status = if_else(is.na(date_of_frontotemporal_dementia_report),0,1)) |>
    filter(frontotemporal_status == 0) |>
    # Exclude people with other forms of dementia
    mutate(dementia_status = if_else(is.na(date_of_dementia_report),0,1)) |>
    filter((dementia_status == 0 & ad_status == 0) | (dementia_status == 1 & ad_status == 1))
  
  # Genetic data ----
  map <- readr::read_delim(paste0(pathData,"UKBiobank/selected_snps_c4.map"),col_names = FALSE) %>%
    dplyr::slice(rep(1:dplyr::n(), each = 2)) %>%
    dplyr::group_by(X2) %>%
    dplyr::mutate(SNP = paste0(X2,"_",dplyr::row_number())) %>%
    dplyr::pull(SNP)
  
  ped <- readr::read_delim(paste0(pathData,"UKBiobank/selected_snps_c4.ped"),col_names = FALSE) %>%
    dplyr::select(-"X2":-"X6")
  colnames(ped) <- c("eid", map)
  
  alleles <- readr::read_delim(paste0(pathData,"UKBiobank/selected_snps_c4.bim"),col_names = FALSE) %>%
    dplyr::select("SNP" = "X2", "ALLELE0" = "X5") %>%
    tidyr::pivot_wider(names_from = "SNP", values_from = "ALLELE0") %>%
    dplyr::slice(rep(1:dplyr::n(), each = nrow(ped)))
  snps <- colnames(alleles)
  alleles <- alleles %>%
    dplyr::rename_with(~paste0(.x,"_0")) %>%
    dplyr::mutate(eid = ped$eid) %>%
    dplyr::inner_join(ped)
  
  for (snp in snps) {
    alleles <- alleles %>%
      dplyr::mutate(
        !!snp := dplyr::if_else(
          .data[[paste0(snp, "_0")]] == .data[[paste0(snp, "_1")]], 1, 0
        ) + dplyr::if_else(
          .data[[paste0(snp, "_0")]] == .data[[paste0(snp, "_2")]], 1, 0
        ) + dplyr::if_else(
          .data[[paste0(snp, "_1")]] == "0", NA, 0
        )
      )
  }
  
  alleles %>%
    dplyr::select(eid, all_of(snps), all_of(paste0(snps,"_0"))) %>%
    dplyr::inner_join(ukb_data,by = "eid")
}

loadSnpOutcomeResults <- function(i, data, regression){
  tibble(
    "snp" = i,
    "effect_allele.outcome" = data |> select(all_of(paste0(i,"_0"))) |> distinct() |> pull(),
    "samplesize.outcome" = data |> tally() |> pull(),
    "beta.outcome" = regression[1],
    "se.outcome" = regression[2],
    "eaf.outcome" = (data |> select("snp") |> pull() |> sum())/(2*(data |> tally() |> pull())),
    "pval.outcome" = regression[4]
  )
}

addCoding <- function(bd, variable, dir_data){
  coding <- tibble(read.table(paste0(dir_data,"UKBiobank/",variable,"_coding.tsv"), header = TRUE, sep = "\t")) |>
    select(!!variable := "coding", "meaning")
  
  bd <- bd |>
    left_join(
      coding,
      by = variable
    ) |>
    select(-!!variable) |>
    rename(!!variable := "meaning") |>
    mutate(!!variable := if_else(.data[[variable]] %in% c("Prefer not to answer", "Do not know"), NA, .data[[variable]]))
  return(bd)
}
