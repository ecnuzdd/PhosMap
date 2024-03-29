% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ksea_regulons_info.R
\name{get_ksea_regulons_info}
\alias{get_ksea_regulons_info}
\title{Get informational data frame by combining results from all experiments}
\usage{
get_ksea_regulons_info(
  ksea_regulons,
  ksea_trans_list,
  ksea_x_list,
  ptypes_data_ratio_colnames
)
}
\arguments{
\item{ksea_regulons}{A kinase vector from all experiments.}

\item{ksea_trans_list}{A list that consits of regulation direction of kinase from each experiment by ksea.}

\item{ksea_x_list}{A list that consits of sepecific information from each experiment by ksea, like regulation direction, p value and activity etc..}

\item{ptypes_data_ratio_colnames}{A vector that consists of column names from experiments.}
}
\value{
A data frame containing sepecific information of all experiments from ksea results, like regulation direction, pvalue and activity etc..
}
\description{
Get informational data frame by combining results from all experiments
}
\examples{
\dontrun{
ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_ksea_regulons_info.RData"
load_data <- load_data_with_ftp(ftp_url, 'RData')
writeBin(load_data, "get_ksea_regulons_info.RData")
load("get_ksea_regulons_info.RData")

ksea_regulons_activity_df <- get_ksea_regulons_info(
  ksea_regulons,
  ksea_trans_list,
  ksea_activity_list,
  ptypes_data_ratio_colnames
)
ksea_regulons_activity_df
}
}
\author{
Dongdong Zhan and Mengsha Tong
}
