% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_aligned_seq_for_mea.R
\name{get_aligned_seq_for_mea}
\alias{get_aligned_seq_for_mea}
\title{Taking S/T/Y as the center, align sequence to fasta library by the specific length.}
\usage{
get_aligned_seq_for_mea(
  ID,
  Sequence,
  AA_in_protein,
  fixed_length,
  species = "human",
  fasta_type = "refseq"
)
}
\arguments{
\item{ID}{A vector for gi number of proteins.}

\item{Sequence}{A vector for sequence of peptides.}

\item{AA_in_protein}{A vector for the locations of S/T/Y in sequence of proteins.}

\item{fixed_length}{Length of aligned sequence,the default is 15.}

\item{species}{A string for the library of species, the options are human, mouse and rattus, the default is human.}

\item{fasta_type, }{A string for fasta source, the options are refseq and uniprot, the default is refseq}
}
\value{
A data frame containing ID, Sequence, AA_in_protein, aligned_seq.
}
\description{
Taking S/T/Y as the center, align sequence to fasta library by the specific length.
}
\examples{
## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
## It may take a few minutes.
if(FALSE){
ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_aligned_seq_for_mea.RData"
load_data <- load_data_with_ftp(ftp_url, 'RData')
writeBin(load_data, "get_aligned_seq_for_mea.RData")
load("get_aligned_seq_for_mea.RData")

foreground_df <- get_aligned_seq_for_mea(
  ID[1:100], Sequence[1:100], AA_in_protein[1:100],
  fixed_length, species = 'human',
  fasta_type = 'refseq'
)
head(foreground_df)

}
}
\references{
Hadley Wickham (2018). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.3.0.\
https://CRAN.R-project.org/package=stringr.
}
\author{
Dongdong Zhan and Mengsha Tong
}
