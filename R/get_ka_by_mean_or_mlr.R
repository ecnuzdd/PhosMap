#' Computing kinase activity using mean value and multiple linear regression (ridge regression) except KSEA
#'
#' @param ptypes_data A data frame of phosphorylation data after normalization.
#' @param species A string representing the species of imported data, the options are human, mouse and rat.
#' @param log2_label A boolean value representing whether data is logarithmetics, the default is FALSE.
#' @param method A string for the method to compute kinase activity, the options are 'mean' and 'mlr' (multiple linear regression),
#' the default is mean.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame that consists of kinase, psite, substrate, counting byond ratio_cutoff and corresponding original value.
#'
#' @export
#'
#' @examples
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_ka_by_mean_or_mlr.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_ka_by_mean_or_mlr.RData")
#' load("get_ka_by_mean_or_mlr.RData")
#'
#' kinase_activity_df <- get_ka_by_mean_or_mlr(
#'   cluster_df,
#'   species = 'human',
#'   log2_label = TRUE,
#'   method = 'mean'
#' )
#' head(kinase_activity_df)
#'

get_ka_by_mean_or_mlr <- function(
  ptypes_data,
  species = 'human',
  log2_label = FALSE,
  method = 'mean'
){
  requireNamespace('utils')
  requireNamespace('stats')
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

  ksea_regulons <- unique(
    as.vector(
      unlist(kinase_substrate_regulation_relationship[,1])
    )
  )

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
    if(length(index_of_kinases_i) > 0){
      kinase_substrate_i_df <- kinase_substrate[index_of_kinases_i, c(1,3,2)] # kinase, site, substrate
      kinases_site_substrate_i <- paste(kinase_substrate_i_df[,1], kinase_substrate_i_df[,2], kinase_substrate_i_df[,3], sep = '|')
      site_substrate_i <-  paste(kinase_substrate_i_df[,2], kinase_substrate_i_df[,3], sep = '|')
      kinases_site_substrate <- c(kinases_site_substrate, kinases_site_substrate_i)
      site_substrate <- c(site_substrate, site_substrate_i)
    }
  }

  regulons_list <- list()
  regulons_list_index <- 0
  ksea_regulons_count <- length(ksea_regulons)
  regulons_count <- NULL
  regulons_list_names <- NULL
  for(i in seq_len(ksea_regulons_count)){
    ksea_regulon <- ksea_regulons[i]
    index_of_match <- which(grepl(ksea_regulon, kinases_site_substrate))
    if(length(index_of_match)>0){
      regulons_list_index <- regulons_list_index + 1
      regulons <- site_substrate[index_of_match]
      regulons_list[[regulons_list_index]] <- regulons
      regulons_count <- c(regulons_count, length(index_of_match))
      regulons_list_names <- c(regulons_list_names, ksea_regulon)
    }
  }

  regulons_list <- lapply(regulons_list, function(x){gsub(' ', '', x)})
  names(regulons_list) <- regulons_list_names

  if(method == 'mean'){
    regulons_list_count <- length(regulons_list)
    kinase_site_substrate_activity_df <- NULL
    for(i in seq_len(regulons_list_count)){
      ksea_regulon_i <- regulons_list_names[i]
      substrate_i <- regulons_list[[i]]
      substrate_count_i <- length(substrate_i)
      index_of_match_i <- NULL
      for(j in seq_len(substrate_count_i)){
        substrate_j <- substrate_i[j]
        index_of_match_j <- which(sites_id==substrate_j)
        if(length(index_of_match_j)>1){
          index_of_match_j <- index_of_match_j[1] # keep only one
        }
        ptypes_data_ratio_j <- ptypes_data_ratio[index_of_match_j,]
        index_of_match_i <- c(index_of_match_i, index_of_match_j)
      }
      if(length(index_of_match_i)==1){
        substrate_ratio_all_exps <- t(ptypes_data_ratio[c(index_of_match_i),])
      }else{
        substrate_ratio_all_exps <- ptypes_data_ratio[index_of_match_i,]
      }
      ksea_regulon_i_activity <- 2^colMeans(substrate_ratio_all_exps) # The original value, not log2
      kinase_site_substrate_activity_df <- rbind(kinase_site_substrate_activity_df, ksea_regulon_i_activity)
    }
    result_df <- data.frame(regulons_list_names, kinase_site_substrate_activity_df)
    colnames(result_df) <- c('Kinase', colnames(kinase_site_substrate_activity_df))
    return(result_df)
  }else if(method == 'mlr'){
    requireNamespace("glmnet")
    x_vector <- sort(regulons_list_names)
    x_vector_count <- length(x_vector)
    y_vector <- sort(unique(as.vector(unlist(regulons_list))))
    y_vector_count <- length(y_vector)
    y_vector_assign_value_count <- rep(0, y_vector_count)

    xy_mat1 <- matrix(0, y_vector_count, x_vector_count, dimnames = list(y_vector, x_vector))
    xy_mat2 <- matrix(0, y_vector_count, ncol(ptypes_data_ratio), dimnames = list(y_vector, colnames(ptypes_data_ratio)))

    for(i in seq_len(x_vector_count)){
      regulons <- regulons_list_names[i]
      i_substrates <- regulons_list[[i]]
      xy_mat1_match_index <- which(y_vector == ij_substrate)
      xy_mat1[xy_mat1_match_index, regulons] <- 1
      i_substrates_count <- length(i_substrates)
      for(j in seq_len(i_substrates_count)){
        ij_substrate <- i_substrates[j]
        ij_match_index <- which(sites_id == ij_substrate)
        ij_match_value <- as.vector(unlist(ptypes_data_ratio[ij_match_index,]))
        xy_mat2_match_index = which(y_vector == ij_substrate)
        if(y_vector_assign_value_count[xy_mat2_match_index] == 0){
          xy_mat2[xy_mat2_match_index,] <- ij_match_value
        }else{
          xy_mat2[xy_mat2_match_index,] <- xy_mat2[xy_mat2_match_index,] + ij_match_value
        }
        y_vector_assign_value_count[xy_mat2_match_index] <- y_vector_assign_value_count[xy_mat2_match_index] + 1
      }
    }
    xy_mat2 <- xy_mat2/y_vector_assign_value_count
    # xy_mat2 <- 2^xy_mat2
    x <- xy_mat1
    coefV_df <- NULL
    for(i in seq_len(ncol(xy_mat2))){
      y <- as.vector(unlist(xy_mat2[,i]))
      # alpha = 0 -> ridge regression
      cv_fit_ridge <- glmnet::cv.glmnet(
        x, y, alpha = 0,
        intercept = FALSE, grouped=FALSE,
        thresh = 0.001,  family="gaussian",
        type.measure = "mse", standardize = TRUE,
        standardize.response = TRUE
      )
      # plot(cv_fit_ridge)
      coefV <- as.vector(stats::coef(cv_fit_ridge))
      coefV_df <- cbind(coefV_df, coefV)
    }
    result_df <- data.frame(regulons_list_names, coefV_df[-1,])
    colnames(result_df) <- c('Kinase', colnames(xy_mat2))
    return(result_df)

  }else{
    stop('The input parameters may be wrong.')
  }


}
