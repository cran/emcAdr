#' Firth Penalized Logistic Regression for Drug Cocktails
#'
#' This function prepares a specific "cocktail" (a set of drugs) 
#' and performs a Firth's penalized logistic regression to estimate the interaction 
#' effect between the drug present in the combination detected as "at risk".
#'
#' @param cocktail An integer vector representing the ATC indices of drugs in the combination.
#' @param upper_bound A list or vector defining the hierarchy/bounds (upper_bound column of the ATC_tree).
#' @param patient_data A data frame containing patient-level data, including the ADR outcome.
#' @param adr_column A string specifying the column name in \code{patient_data} used as the dependent variable (Y). 
#'   Defaults to "patientADR".
#'
#' @return An object of class \code{logistf} containing the regression results, including 
#'   coefficients, p-values, and confidence intervals.
#'
#' @details 
#' Firth's method is preferred here as it handles 
#' "separation" issues common in sparse clinical data (where a drug combination 
#' might perfectly predict an ADR).
#'
#' @examples
#' \dontrun{
#' # Example using indices for drugs 888, 659
#' results <- run_firth_regression(
#'   cocktail = c(888, 659),
#'   upper_bound = ATC_Tree_UpperBound_2024$upperBound,
#'   patient_data = FAERS_myopathy
#' )
#' summary(results)
#' }
#' 
#' @importFrom logistf logistf
#' @export
run_firth_regression <- function(cocktail, 
                                 upper_bound, 
                                 patient_data, 
                                 adr_column = "patientADR") {
  
  df_model <- combination_data_frame(
    cocktail, 
    upper_bound, 
    patient_data
  )
  
  if (!adr_column %in% names(patient_data)) {
    stop(paste("Column", adr_column, "not found in patient_data"))
  }
  
  df_model$Y <- patient_data[[adr_column]]
  
  model_firth <- logistf::logistf(Y ~ ., data = df_model)
  
  return(model_firth)
}