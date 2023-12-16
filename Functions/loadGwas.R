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