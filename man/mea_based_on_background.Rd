% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mea_based_on_background.R
\name{mea_based_on_background}
\alias{mea_based_on_background}
\title{Motif enrichment based on global background (fasta library from Refseq).}
\usage{
mea_based_on_background(foreground, AA_in_protein, background, motifx_pvalue)
}
\arguments{
\item{foreground}{A vector for aligned sequence of foreground.}

\item{AA_in_protein}{A vector for the location of S/T/Y in sequence of protein.}

\item{background}{A vector for aligned sequence of background.}

\item{motifx_pvalue}{A numeric value for selecting motifs that meets the minimum cutoff.}
}
\value{
A list containing motifs and the corresponding sequences
}
\description{
Motif enrichment based on global background (fasta library from Refseq).
}
\examples{
## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
## It may take a few minutes.
if(FALSE){
    ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/mea_based_on_background.RData"
    load_data <- load_data_with_ftp(ftp_url, 'RData')
    writeBin(load_data, "mea_based_on_background.RData")
    load("mea_based_on_background.RData")

    motifs_list <- mea_based_on_background(
      foreground[1:100],
      AA_in_protein[1:100],
      background[1:1000],
      motifx_pvalue
    )
}

}
\author{
Dongdong Zhan and Mengsha Tong
}
