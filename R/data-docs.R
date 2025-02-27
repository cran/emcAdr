#' ATC Tree Upper Bound 2024
#'
#' Example dataset representing the ATC tree structure, sourced from the WHO website (2024-02-23).
#' This dataset is provided for demonstration and testing purposes with the package.
#'
#' @format A data frame with 4 variables:
#' \describe{
#'   \item{ATCCode}{The code of ATC nodes}
#'   \item{Name}{The name of ATC nodes}
#'   \item{ATC_length}{The number of characters in the ATCCode}
#'   \item{upperBound}{The index of the last child node in the tree}
#' }
#' @source World Health Organization, ATC classification register
"ATC_Tree_UpperBound_2024"

#' FAERS Myopathy Dataset
#'
#' Example dataset representing drug intake and adverse event reports from FAERS.
#' This dataset is provided to demonstrate the functionality of genetic and MCMC algorithms in the package.
#'
#' @format A data frame with 2 columns:
#' \describe{
#'   \item{patientATC}{Drug intake for each patient as a vector of ATC tree indices}
#'   \item{patientADR}{Indicates if the patient experienced myopathy as an adverse event}
#' }
#' @source Food & Drug Administration Event Reporting System (FAERS)
"FAERS_myopathy"
