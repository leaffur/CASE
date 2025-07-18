% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CASE_models.R
\name{get_credible_sets}
\alias{get_credible_sets}
\title{CASE Obtain Credible Sets}
\usage{
get_credible_sets(
  pips,
  R,
  verbose = TRUE,
  cor.min = 0.5,
  coverage_thres = 0.95,
  ruled_out = 1e-04
)
}
\arguments{
\item{pips}{(M * C),The pips of SNPs.}

\item{R}{M * M matrix of LD.}

\item{verbose}{(optional) logical, whether to print logging information. Default = TRUE.}

\item{cor.min}{minimum correlation in the credible sets}

\item{coverage_thres}{threshold for the sum of PIPs.}

\item{ruled_out}{excluding SNPs with PIPs less than the threshold.}
}
\value{
a length C list of credible sets.
}
\description{
Obtain credible sets for any multi-trait fine-mapping results.
}
