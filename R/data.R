#' Example datasets
#'
#' This data list, \code{example.data}, includes three datasets generated based on Saturation mutagenesis results (M. Kircher, et al.,2019) and the DeepSEA features (Zhou & Troyanskaya, 2015).
#' The training and testing datasets in the data list include binary response vectors, which are truncations of the P values of tissue K562 from the Saturation mutagenesis results,
#' and reduced versions of the DeepSEA features for a faster computational demonstration. The \code{full.data} dataset includes the original P values, chromosome and allelic information, and the complete DeepSEA features.
#'
#'
#' @format The \code{example.data$train.data} and \code{example.data$test.data} are dataframes with 220 and 1574 observations and 146 variables.
#' \describe{
#'   \item{response}{A binary response vector}
#'   \item{features}{Standardized 145 DeepSEA features}
#' }
#'
#' @format The \code{example.data$full.data} is a dataframe with 1794 observations and 924 variables, i.e., including all 919 DeepSEA features.
#' \describe{
#'   \item{chr}{The chromosome of SNPs}
#'   \item{pos}{The position of SNPs}
#'   \item{ref.alt}{The reference and alternative alleles of SNPs}
#'   \item{p.value}{The P value of SNPs}
#'   \item{features}{The original 919 DeepSEA features}
#' }
"example.data"

