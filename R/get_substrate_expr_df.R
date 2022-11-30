#' Get a data frame only containing kinase inferred by KSEA
#'
#' @param ID A phosporylation ID vector like VIM_S56 (GeneSymbol_psite).
#' @param kinase_substrate_regulation_relationship A data frame contanning relationship of
#' kinase-substrate that consists of "kinase", "substrate", "site", "sequence" and "predicted" columns.
#' @param ksea_regulons A kinase vector from ksea
#' @param ptypes_data_ratio A data frame that the ratio of phosphorylation and profiling data
#' @param ratio_cutoff A cutoff that depicts quantification changes
#' at phosphorylation level relative to profiling level, the default is 3.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @return A data frame that consists of kinase, psite, substrate,
#' counting byond ratio_cutoff and corresponding original value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_substrate_expr_df.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_substrate_expr_df.RData")
#' load("get_substrate_expr_df.RData")
#'
#' kinase_site_substrate_original_ratio_df <- get_substrate_expr_df(
#'   ID,
#'   kinase_substrate_regulation_relationship,
#'   ksea_regulons,
#'   ptypes_data_ratio,
#'   ratio_cutoff = 3
#' )
#' head(kinase_site_substrate_original_ratio_df)
#' }

get_substrate_expr_df <- function(
  ID,
  kinase_substrate_regulation_relationship,
  ksea_regulons,
  ptypes_data_ratio,
  ratio_cutoff = 3
){

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

  kinase_substrate <- kinase_substrate_regulation_relationship
  kinases_site_substrate <- NULL
  site_substrate <- NULL
  sites_id_count <- length(sites_id)
  for(i in seq_len(sites_id_count)){
    # cat('\n complete: ', i, '/', sites_id_count)
    substrate_site <- as.vector(ptypes_data_ID[i,1])
    substrate_symbol <- as.vector(ptypes_data_ID[i,2])
    # extract kinase from a table called relationship of kinase-substrate
    index_of_kinases_i <- which(kinase_substrate[,3]==substrate_site & kinase_substrate[,2]==substrate_symbol)
    kinase_substrate_i_df <- kinase_substrate[index_of_kinases_i, c(1,3,2)] # kinase, site, substrate
    kinases_site_substrate_i <- paste(kinase_substrate_i_df[,1], kinase_substrate_i_df[,2], kinase_substrate_i_df[,3], sep = '|')
    site_substrate_i <-  paste(kinase_substrate_i_df[,2], kinase_substrate_i_df[,3], sep = '|')
    kinases_site_substrate <- c(kinases_site_substrate, kinases_site_substrate_i)
    site_substrate <- c(site_substrate, site_substrate_i)
  }

  regulons_list <- list()
  ksea_regulons_count <- length(ksea_regulons)
  regulons_count <- NULL
  for(i in seq_len(ksea_regulons_count)){
    ksea_regulon <- ksea_regulons[i]
    index_of_match <- which(grepl(ksea_regulon, kinases_site_substrate))
    regulons <- site_substrate[index_of_match]
    regulons_list[[i]] <- regulons
    regulons_count <- c(regulons_count, length(index_of_match))
  }
  regulons_list <- lapply(regulons_list, function(x){gsub(' ', '', x)})
  names(regulons_list) <- ksea_regulons

  ksea_regulons_with_count_greater_than_zero <- ksea_regulons[regulons_count>0]
  ksea_regulons_with_count_greater_than_zero_count <- length(ksea_regulons_with_count_greater_than_zero)

  ptypes_data_ratio_e2 <- 2^ptypes_data_ratio
  ptypes_data_ratio_e2 <- data.frame(ptypes_data_ratio_e2)
  RATIO_CUTOFF <- ratio_cutoff
  kinase_site_substrate_ratio_df <- c()

  for(i in seq_len(ksea_regulons_with_count_greater_than_zero_count)){
    ksea_regulon_i <- ksea_regulons_with_count_greater_than_zero[i]
    index_of_substrate_i <- which(ksea_regulons==ksea_regulon_i)
    substrate_i <- regulons_list[[index_of_substrate_i]]
    substrate_count_i <- length(substrate_i)
    index_of_match_i <- NULL
    count_of_greater_than_ratio_cutoff <- NULL
    count_of_less_than_ratio_cutoff <- NULL
    for(j in seq_len(substrate_count_i)){

      substrate_j <- substrate_i[j]
      index_of_match_j <- which(sites_id==substrate_j)
      if(length(index_of_match_j)>1){
        index_of_match_j <- index_of_match_j[1] # keep only one
      }

      ptypes_data_ratio_e2_j <- ptypes_data_ratio_e2[index_of_match_j,]


      count_of_greater_than_ratio_cutoff_j <- length(which(ptypes_data_ratio_e2_j > RATIO_CUTOFF))
      count_of_less_than_ratio_cutoff_j <- length(which(ptypes_data_ratio_e2_j < 1/RATIO_CUTOFF))
      count_of_greater_than_ratio_cutoff <- c(count_of_greater_than_ratio_cutoff, count_of_greater_than_ratio_cutoff_j)
      count_of_less_than_ratio_cutoff <- c(count_of_less_than_ratio_cutoff, count_of_less_than_ratio_cutoff_j)

      index_of_match_i <- c(index_of_match_i, index_of_match_j)
    }
    if(length(index_of_match_i)==1){
      substrate_ratio_all_exps <- ptypes_data_ratio_e2[c(index_of_match_i),]
    }else{
      substrate_ratio_all_exps <- ptypes_data_ratio_e2[index_of_match_i,] # The original value, not log2
    }

    site <- apply(data.frame(substrate_i), 1, function(x){
      x <- strsplit(x, '[|]')[[1]][1] # strsplit(x, '\\|')
      x
    })
    substrate <- apply(data.frame(substrate_i), 1, function(x){
      x <- strsplit(x, '[|]')[[1]][2] # strsplit(x, '\\|')
      x
    })

    kinase <- rep(ksea_regulon_i, substrate_count_i)
    df <- cbind(kinase, site, substrate,
               count_of_greater_than_ratio_cutoff, count_of_less_than_ratio_cutoff,
               substrate_ratio_all_exps)
    kinase_site_substrate_ratio_df <- rbind(kinase_site_substrate_ratio_df, df)
  }
  kinase_site_substrate_original_ratio_df <- kinase_site_substrate_ratio_df
  return(kinase_site_substrate_original_ratio_df)
}
