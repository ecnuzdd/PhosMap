#' Get a data frame only containing inforamtion of kinase inferred by KSEA
#'
#' @param ptypes_data A data frame of phosphorylation data after normalization.
#' @param species A string representing the species of imported data, the options are human, mouse and rat.
#' @param log2_label A boolean value representing whether data is logarithmetics, the default is FALSE.
#' @param ratio_cutoff A cutoff that depicts quantification changes at phosphorylation level relative to profiling level, the default is 3.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame that consists of kinase, psite, substrate, counting byond ratio_cutoff and corresponding original value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summary_df_list_from_ksea <- get_summary_from_ksea(
#'   ptypes_data,
#'   species = 'human',
#'   log2_label = TRUE,
#'   ratio_cutoff = 3
#' )
#' }
#'
#'
#'

get_summary_from_ksea <- function(
  ptypes_data,
  species = 'human',
  log2_label = TRUE,
  ratio_cutoff = 3
){
  requireNamespace('utils')
  # read relationship of kinase-substrate provided by PhosMap
  # KSRR: kinase substrate regulation relationship
  # A data frame contanning relationship of kinase-substrate that consists of "kinase", "substrate", "site", "sequence" and "predicted" columns.
  KSRR_FILE_NAME <- paste(species, 'ksrr.csv', sep = '_')
  KSRR_FILE_PATH <- normalizePath(
    system.file(
      'extdata',
      'kinase_substrate_regulation_relationship_table', species, KSRR_FILE_NAME,
      package = "PhosMap"
    ),
    mustWork = FALSE
  )
  if(!file.exists(KSRR_FILE_PATH)){
    cat(KSRR_FILE_PATH, ' -> ', 'No the file')
    stop('')
  }
  kinase_substrate_regulation_relationship <- utils::read.csv(KSRR_FILE_PATH, header = TRUE, sep= ",", stringsAsFactors = NA)

  ID <- as.vector(ptypes_data[,1])
  ptypes_data_ratio <- ptypes_data[,-1]
  if(!log2_label){

    ptypes_data_ratio <- log2(ptypes_data_ratio)
  }
  ptypes_data_ratio_colnames <- colnames(ptypes_data_ratio)



  ksea_es_list <- list()
  ksea_pvalue_list <- list()
  ksea_regulons_list <- list()
  ksea_activity_list <- list()
  ksea_trans_list <- list()
  ptypes_data_exp_count <- ncol(ptypes_data_ratio)
  cat('\n Starting KSEA')
  for(i in seq_len(ptypes_data_exp_count)){
    cat('\n completing: ', i, '/', ptypes_data_exp_count)
    ptypes_data_ratio_in_sigle_exp <- as.numeric(unlist(ptypes_data_ratio[,i]))
    ksea_result_list_i <- get_ksea_result_list(
      ptypes_data_ratio_in_sigle_exp,
      ID,
      kinase_substrate_regulation_relationship,
      ksea_activity_i_pvalue = 0.05
    )
    ksea_es_list[[i]] <- ksea_result_list_i$ksea_es_i_non_NA
    ksea_pvalue_list[[i]] <- ksea_result_list_i$ksea_pvalue_i_non_NA
    ksea_regulons_list[[i]] <- ksea_result_list_i$ksea_regulons_i_non_NA
    ksea_activity_list[[i]] <- ksea_result_list_i$ksea_activity_i
    ksea_trans_list[[i]] <- ksea_result_list_i$ksea_trans_i
    cat('\n completed: ', i, '/', ptypes_data_exp_count)
  }
  cat('\n Ending KSEA')

  cat('\n Extracting information data frame derived from KSEA')
  cat('\n ********** Regulation direction from KSEA **********')
  cat('\n ********** Pvalue from KSEA **********')
  cat('\n ********** Activity from KSEA **********')
  cat('\n ********** Kinase_site_substrate quantification matrix after KSEA **********')
  cat('\n')

  ksea_regulons <- unique(unlist(ksea_regulons_list))
  ksea_regulons_count <- length(ksea_regulons)
  # enrichment score from ksea
  # pvalue from ksea
  # regulons (kinase) from ksea
  # kinase activity based on pvalue and enrichment score computed by ksea
  # regulation direction: 1 = activate, 0 = no work, -1 = supress
  ksea_regulons_regulation_direction_df <- get_ksea_regulons_info(ksea_regulons, ksea_trans_list, ksea_trans_list,
                                                  ptypes_data_ratio_colnames)
  ksea_regulons_pvalue_df <- get_ksea_regulons_info(ksea_regulons, ksea_trans_list, ksea_pvalue_list,
                                                   ptypes_data_ratio_colnames)
  ksea_regulons_activity_df <- get_ksea_regulons_info(ksea_regulons, ksea_trans_list, ksea_activity_list,
                                                     ptypes_data_ratio_colnames)

  ksea_kinase_site_substrate_original_ratio_df <- get_substrate_expr_df(ID,
                                                      kinase_substrate_regulation_relationship,
                                                      ksea_regulons,
                                                      ptypes_data_ratio,
                                                      ratio_cutoff)
  summary_df_list_from_ksea <- list(
    ksea_regulons_regulation_direction_df = ksea_regulons_regulation_direction_df, # regulation direction: 1 = activate, 0 = no work, -1 = supress
    ksea_regulons_pvalue_df = ksea_regulons_pvalue_df, # pvalue from ksea
    ksea_regulons_activity_df = ksea_regulons_activity_df, # kinase activity based on pvalue and enrichment score computed by ksea
    ksea_kinase_site_substrate_original_ratio_df = ksea_kinase_site_substrate_original_ratio_df #
  )

  cat('\n KSEA OK! ^_^')

  return(summary_df_list_from_ksea)

}
