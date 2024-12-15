#' Example Data
#'
#' An example data of three traits for the CASE fine-mapping illustration.
#'
#' @format ## `example_data`
#' A list contains 3 elements:
#' \describe{
#'   \item{Y}{500 by 3 matrix of phenotype.}
#'   \item{X}{500 by 1000 matrix of genotype. This is created by R package `CorBin`}
#'   \item{B}{1000 by 3 matrix of eQTL effects. Only the 10th and 950th SNPs have eQTL effects.}
#' }
"example_data"