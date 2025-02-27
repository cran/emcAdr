#' Calculate p-value of sampled value
#'
#' @param empirical_distribution A numeric vector of values representing the empirical distribution (return value of DistributionAproximation function)
#' @param sampled_values A scalar or a vector of real valued number representing the sampled value (score to be tested)
#' @param isFiltered A boolean representing if we want to use the filtered distribution or the distribution as is (False by default)
#' @param includeZeroValue A boolean that indicate if you want to take into account the null score (False by default)
#' @return A numeric value representing the empirical p-value
#' @examples
#' \donttest{
#' data("ATC_Tree_UpperBound_2024")
#' data("FAERS_myopathy")
#' 
#' cocktails = list(c(561, 904),
#'                c(1902, 4585))
#'                
#' estimated_score_distribution = DistributionApproximation(epochs = 10,
#'             ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy)
#'             
#' Hypergeom_of_cocktails = compute_hypergeom_on_list(cocktails = cocktails,
#'                               ATCtree = ATC_Tree_UpperBound_2024, 
#'                               observations = FAERS_myopathy)
#'             
#' p_value = p_value_on_sampled(empirical_distribution = estimated_score_distribution,
#'       sampled_values = Hypergeom_of_cocktails)
#'}
#' @export
p_value_on_sampled <- function(empirical_distribution, sampled_values, isFiltered = FALSE,
                               includeZeroValue = FALSE) {
  # Sort empirical distribution in ascending order (if the distribution comes. from 
  # the histogramToDitribution function it should already be sorted)
  if(isFiltered){
    if(includeZeroValue){
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$Filtered_score_distribution)
      }else{
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$Filtered_score_distribution[2:length(
        empirical_distribution$Filtered_score_distribution)])
      }
    empirical_distribution_array <- append(empirical_distribution_array, empirical_distribution$OutstandingRR)
  }
  else{
    if(includeZeroValue){
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$ScoreDistribution)
    }else{
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$ScoreDistribution[2:length(
        empirical_distribution$ScoreDistribution)])
    }
    empirical_distribution_array <- append(empirical_distribution_array, empirical_distribution$OutstandingRR)
  }
  
  empirical_distribution_array <- sort(empirical_distribution_array)
  
  p_values <- numeric(length(sampled_values))
  
  # Iterate over each sampled value
  for (i in seq_along(sampled_values)) {
    # Calculate ECDF value for the current sampled value
    ecdf_value <- sum(empirical_distribution_array <= sampled_values[i]) / length(empirical_distribution_array)
    
    # Calculate p-value
    p_values[i] <- 1 - ecdf_value
  }
  
  return(p_values) 
}

#' Calculate the divergence between 2 distributions (the true Distribution and the learned one)
#' @param empirical_distribution A numeric vector of values representing the empirical distribution (return value of DistributionAproximation function)
#' @param true_distribution A numeric vector of values representing the true distribution computed by the trueDistributionSizeTwoCocktail function
#' @param method A string, either "TV" or "KL" to respectively use the total variation distance or the Kullback-Leibler divergence. (default = "TV")
#' @param Filtered Should we use the filtered distribution or the normal one
#' @return A numeric value representing the divergence of the 2 distributions
#' @examples
#' \donttest{
#' data("ATC_Tree_UpperBound_2024")
#' data("FAERS_myopathy")
#' 
#' estimated_score_distribution = DistributionApproximation(epochs = 10,
#'             ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy[1:100,], Smax =2)
#'             
#' true_score_distribution = trueDistributionSizeTwoCocktail(ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy[1:100,], beta = 4)
#' 
#' divergence <- calculate_divergence(empirical_distribution = estimated_score_distribution,
#'                 true_distribution = true_score_distribution)
#'}
#' @export
calculate_divergence <- function(empirical_distribution, true_distribution, method = "TV", Filtered = FALSE){
  RRmax <- (length(empirical_distribution$ScoreDistribution) - 1) / 10
  
  if(RRmax > 30){
    if(Filtered){
      dist_ouststandingRR <- OutsandingScoreToDistribution(true_distribution$Outstanding_score, RRmax)
      length(true_distribution$Filtered_score_distribution) <- length(dist_ouststandingRR)
      true_distribution$Filtered_score_distribution[is.na(true_distribution$Filtered_score_distribution)] <- 0
      true_distribution$Filtered_score_distribution <- true_distribution$Filtered_score_distribution + dist_ouststandingRR
    }
    else{
      dist_ouststandingRR <- OutsandingScoreToDistribution(true_distribution$Outstanding_score, RRmax)
      length(true_distribution$ScoreDistribution) <- length(dist_ouststandingRR)
      true_distribution$ScoreDistribution[is.na(true_distribution$ScoreDistribution)] <- 0
      true_distribution$ScoreDistribution <- true_distribution$ScoreDistribution + dist_ouststandingRR
    }
  }
  else if(RRmax < 30){
    if(Filtered){
      true_distribution$Filtered_score_distribution[((RRmax*10)+1)] <- sum(true_distribution$Filtered_score_distribution[((RRmax*10)+1):length(true_distribution$Filtered_score_distribution)]) + length(true_distribution$OutstandingRR)
      length(true_distribution$Filtered_score_distribution) <- ((RRmax*10)+1)
    }else{
      true_distribution$ScoreDistribution[((RRmax*10)+1)] <- sum(true_distribution$ScoreDistribution[((RRmax*10)+1):length(true_distribution$ScoreDistribution)]) + length(true_distribution$OutstandingRR)
      length(true_distribution$ScoreDistribution) <- ((RRmax*10)+1)
    }
  }
  else{
    if (Filtered) {
      true_distribution$Filtered_score_distribution[((RRmax*10)+1)] <- length(true_distribution$Outstanding_score)
    }else{
      true_distribution$ScoreDistribution[((RRmax*10)+1)] <- length(true_distribution$Outstanding_score)
    }
  }
  
  #empirical_distribution$ScoreDistribution <- empirical_distribution$ScoreDistribution / sum(empirical_distribution$ScoreDistribution)
  #true_distribution$ScoreDistribution <- true_distribution$ScoreDistribution / sum(true_distribution$ScoreDistribution)
  adjusted_empirical <- empirical_distribution$ScoreDistribution[2:length(empirical_distribution$ScoreDistribution)] / sum(empirical_distribution$ScoreDistribution[2:length(empirical_distribution$ScoreDistribution)])
  adjusted_true <- true_distribution$ScoreDistribution[2:length(true_distribution$ScoreDistribution)] / sum(true_distribution$ScoreDistribution[2:length(true_distribution$ScoreDistribution)])
  
  if(method == "TV"){
    #return(sum(abs(empirical_distribution$Distribution - true_distribution$Distribution)))
    return (sum(abs(adjusted_empirical - adjusted_true)))
  }
  else if(method == "KL"){
    #return(sum(true_distribution$Distribution * log(true_distribution$Distribution / empirical_distribution$Distribution)))
    return(sum(adjusted_true * log(adjusted_true / adjusted_empirical)))
  }
  else{
    stop("The method should be either \"TV\" for Total variation distance or \"KL\" for Kullback-Leibler divergence")
  }
}