% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{dRNASb_pipeline}
\alias{dRNASb_pipeline}
\title{dRNASb pipeline}
\usage{
dRNASb_pipeline(
  data_file_path,
  phenotype_file_path,
  annotation_function_file_path,
  ppi_file_path,
  result_file_prefix = "",
  de_method = "limma",
  norm_method = "log_TMM",
  perform_filter = TRUE,
  hours_in_data = c("0h", "2h", "4h", "8h", "16h", "24h"),
  replicates_in_data = 3,
  logFC_cutoff = 1,
  num_of_clust = 10,
  hub_gene_cutoff = 10,
  betweenness_cutoff = 100
)
}
\arguments{
\item{data_file_path}{full file path of input data file containing read counts in (gene x samples) format
should have samples corresponding to different hours and replicates for each hour}

\item{phenotype_file_path}{full file path of phenotype file - contains mapping of sample to groups with groups column}

\item{annotation_function_file_path}{full file path of annotation function file}

\item{ppi_file_path}{full file path of ppi file}

\item{result_file_prefix}{a prefix string to be added to all result files}

\item{de_method}{differential expression method to be used}

\item{norm_method}{normalization method to be used - use show_allowed_norm_methods() for available norm methods}

\item{perform_filter}{Should filter preprocessing step be performed}

\item{hours_in_data}{Vector of different hours information present in the data}

\item{replicates_in_data}{Number of replicates for each hour in the data}

\item{logFC_cutoff}{log fold change cutoff used to select differentially expressed genes}

\item{num_of_clust}{number of clusters in mfuzz clustering}

\item{hub_gene_cutoff}{cutoff value in network analysis to determine hub genes}

\item{betweenness_cutoff}{cutoff value in network analysis to determine network betweenness}
}
\description{
dRNASb pipeline
}
