#' Kinase activity analysis based on known and predicted kinase-substrate relationships
#'
#' @param ptypes_data_ratio_in_single_exp A quantification vector from a single experiment.
#' @param ID A phosporylation ID vector like VIM_S56 (GeneSymbol_psite).
#' @param kinase_substrate_regulation_relationship A data frame contanning kinase-substrate relationships that consists of "kinase", "substrate", "site", "sequence" and "predicted" columns.
#' @param ksea_activity_i_pvalue A cutoff used for filtering significant activities computed from KSEA.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#'
#' @return A list containing results from ksea.
#'
#' @export
#'
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_ksea_result_list.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_ksea_result_list.RData")
#' load("get_ksea_result_list.RData")
#'
#' ksea_result_list_i <- get_ksea_result_list(
#'  ptypes_data_ratio_in_single_exp, ID,
#'  kinase_substrate_regulation_relationship,
#'  ksea_activity_i_pvalue = 0.05
#' )
#' head(ksea_result_list_i)
#' }


get_ksea_result_list <- function(ptypes_data_ratio_in_single_exp, ID, kinase_substrate_regulation_relationship, ksea_activity_i_pvalue = 0.05){
  symbol <- apply(data.frame(ID), 1, function(x){
    x <- strsplit(x, split = '_')[[1]]
    x[1]
  })
  site <- apply(data.frame(ID), 1, function(x){
    x <- strsplit(x, split = '_')[[1]]
    x[2]
  })
  ptypes_data_ID <- data.frame(site, symbol)

  sites_id <- paste(site, symbol, sep = '|')
  names(ptypes_data_ratio_in_single_exp) <- sites_id
  sites_id_count <- length(sites_id)

  index_of_ptypes_data_ratio_in_single_exp_desc <- order(ptypes_data_ratio_in_single_exp, decreasing=TRUE)
  ptypes_data_ratio_in_single_exp_desc <- ptypes_data_ratio_in_single_exp[index_of_ptypes_data_ratio_in_single_exp_desc]
  ptypes_data_ratio_in_single_exp_desc_names <- names(ptypes_data_ratio_in_single_exp_desc)

  kinases_i <- NULL
  kinases_site_substrate_i <- NULL
  site_substate_i <- NULL
  site_quant_ratio_i <- NULL
  kinase_substrate <- kinase_substrate_regulation_relationship
  # j <- index_of_sites_id
  for(j in seq_len(sites_id_count)){
    substrate_site <- as.vector(ptypes_data_ID[j,1])
    substrate_symbol <- as.vector(ptypes_data_ID[j,2])
    # extract kinase from a table called relationship of kinase-substrate
    index_of_kinases_i_j <-  which(kinase_substrate[,3]==substrate_site & kinase_substrate[,2]==substrate_symbol)

    kinase_substrate_i_j_df <- kinase_substrate[index_of_kinases_i_j, c(1,3,2)] # kinase, site, substrate
    kinases_i_j <- as.vector(kinase_substrate_i_j_df[,1])
    kinases_site_substrate_i_j <- paste(kinase_substrate_i_j_df[,1], kinase_substrate_i_j_df[,2], kinase_substrate_i_j_df[,3], sep = '|')
    site_substate_i_j <- paste(kinase_substrate_i_j_df[,2], kinase_substrate_i_j_df[,3], sep = '|')
    site_quant_ratio_i_j <- rep(ptypes_data_ratio_in_single_exp[j], nrow(kinase_substrate_i_j_df))

    kinases_i <- c(kinases_i, kinases_i_j)
    kinases_site_substrate_i <- c(kinases_site_substrate_i, kinases_site_substrate_i_j)
    site_substate_i <- c(site_substate_i, site_substate_i_j)
    site_quant_ratio_i <- c(site_quant_ratio_i, site_quant_ratio_i_j)
  }


  kinases_i_unique <- unique(kinases_i)
  kinases_i_unique_count <- length(kinases_i_unique)

  regulons_i <- list()
  kinase_regulation_direction_i <- NULL
  for(k in seq_len(kinases_i_unique_count)){
    kinases_i_unique_k <- kinases_i_unique[k]
    index_regulons_i_k <- which(grepl(kinases_i_unique_k, kinases_site_substrate_i))
    regulons_i_k <- site_substate_i[index_regulons_i_k]
    regulons_i[[k]] <-  regulons_i_k
    site_quant_ratio_i_k <- site_quant_ratio_i[index_regulons_i_k]
    if(length(site_quant_ratio_i_k)>0){
      site_quant_ratio_i_k_mean <- mean(site_quant_ratio_i_k)
      if(site_quant_ratio_i_k_mean > 0){
        direction_i_k <- 1
      }else{
        direction_i_k <- -1
      }
    }else{
      direction_i_k <- NA
    }
    kinase_regulation_direction_i <- c(kinase_regulation_direction_i, direction_i_k)
  }

  regulons_i <- lapply(regulons_i, function(x){gsub(' ','',x)})
  names(regulons_i) <- kinases_i_unique

  ksea_es_i <- NULL
  ksea_pvalue_i <- NULL
  ksea_regulons_i <- NULL
  for(l in seq_len(kinases_i_unique_count)){
    regulons_i_l <- regulons_i[[l]]
    if(TRUE){
      ksea_result_i_l <- get_kses(
        ptypes_data_ratio_in_single_exp_desc,
        regulons_i_l,
        1000
      )
      if(length(ksea_result_i_l) > 1){
        es <- as.vector(ksea_result_i_l$expected_enrichment_score)
        pvalue <- as.vector(ksea_result_i_l$pvalue)
        id <- paste(names(regulons_i[l]), '|', NULL)
      }else{
        es <- NA
        pvalue <- NA
        id <- paste(names(regulons_i[l]), '|', NA)
      }
    }
    ksea_es_i <- c(ksea_es_i, es)
    ksea_pvalue_i <- c(ksea_pvalue_i, pvalue)
    ksea_regulons_i <- c(ksea_regulons_i, id)
  }

  index_of_non_NA <- which(!is.na(ksea_es_i))

  ksea_es_i_non_NA <- ksea_es_i[index_of_non_NA]

  ksea_pvalue_i_non_NA <- ksea_pvalue_i[index_of_non_NA]
  index_of_zero <- which(ksea_pvalue_i_non_NA==0)
  if(length(index_of_zero)>0){
    min_value_in_ksea_pvalue_i_non_NA <- min(ksea_pvalue_i_non_NA[-index_of_zero])
    ksea_pvalue_i_non_NA[index_of_zero] <- min_value_in_ksea_pvalue_i_non_NA*0.1
  }

  ksea_regulons_i_non_NA <- ksea_regulons_i[index_of_non_NA]
  ksea_regulons_i_non_NA <- apply(data.frame(ksea_regulons_i_non_NA), 1, function(x){
    x <- gsub('\\|','',x)
    x <- gsub(' ','',x)
    x
  })

  kinase_regulation_direction_i_non_NA <- kinase_regulation_direction_i[index_of_non_NA]

  ksea_activity_i <- (-log10(ksea_pvalue_i_non_NA)*kinase_regulation_direction_i_non_NA)
  names(ksea_activity_i) <- ksea_regulons_i_non_NA

  # ksea_activity_i_pvalue <- 0.05
  ksea_activity_i_cutoff <- (-log10(ksea_activity_i_pvalue))
  ksea_trans_i <- apply(data.frame(ksea_activity_i), 1, function(x, cutoff){
    if(x > cutoff){
      return(1)
    }else if(x < (-cutoff)){
      return(-1)
    }else{
      return(0)
    }
  }, cutoff = ksea_activity_i_cutoff)

  names(ksea_trans_i) <- ksea_regulons_i_non_NA

  result_list <- list(
    ksea_es_i_non_NA = ksea_es_i_non_NA, # enrichment score from ksea
    ksea_pvalue_i_non_NA = ksea_pvalue_i_non_NA, # pvalue from ksea
    ksea_regulons_i_non_NA = ksea_regulons_i_non_NA, # regulons (kinase) from ksea
    ksea_activity_i = ksea_activity_i, # kinase activity based on pvalue and enrichment score computed by ksea
    ksea_trans_i = ksea_trans_i # regulation direction: 1 = activate, 0 = no work, -1 = supress
  )
  return(result_list)
}
