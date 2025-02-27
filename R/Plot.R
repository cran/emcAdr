#' Plot the evolution of the mean and the best value of the population used by the GeneticAlgorithm
#'
#' @param list A list with 2 elements returned by the GeneticAlgorithm: "mean" and "best", containing the numeric vectors representing the mean and best fitness of the population
#' @param mean_color A string specifying the color of the mean values
#' @param best_color A string specifying the color of the best values
#' @param xlab A string specifying the label for the x-axis
#' @param ylab A string specifying the label for the y-axis
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom grid units
#' @importFrom rlang .data
#' 
#' @examples
#' \donttest{
#' data("ATC_Tree_UpperBound_2024")
#' data("FAERS_myopathy")
#' 
#' results = GeneticAlgorithm(epochs = 10, nbIndividuals = 10, 
#'             ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy)
#' 
#' plot_evolution(list = results)
#'}
#' @return no returned value, should plot the evolution of the genetic algorithm
#' results (mean/max score for each epoch).
#' @export
#'
plot_evolution <- function(list, mean_color = "#F2A900", best_color = "#008080", xlab = "Epochs", ylab = "Score") {
  requireNamespace("dplyr") 
  requireNamespace("ggplot2")
  
  epochs <- seq_along(list$meanFitnesses)
  
  data <- data.frame(epochs = epochs, mean = list$meanFitnesses, best = list$BestFitnesses)  
  
  ggplot2::ggplot(data, ggplot2::aes(x = .data$epochs)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$best, color = "Best"), size = 2, shape = 16, fill = "white") +
    ggplot2::geom_point(ggplot2::aes(y = .data$mean, color = "Mean"), size = 2, shape = 16, fill = "white") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$epochs, 
                                       xend = dplyr::lead(.data$epochs, default = dplyr::last(.data$epochs)), 
                                       y = .data$mean, 
                                       yend = dplyr::lead(.data$mean, default = dplyr::last(.data$mean))), 
                          color = mean_color, size = 0.5, linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$epochs, 
                                       xend = dplyr::lead(.data$epochs, default = dplyr::last(.data$epochs)), 
                                       y = .data$best, 
                                       yend = dplyr::lead(.data$best, default = dplyr::last(.data$best))), 
                          color = best_color, size = 0.5, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c(best_color, mean_color)) +
    ggplot2::labs(title = "Evolution of the population", x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::unit(16, "pt"), margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::unit(14, "pt")),
      axis.text = ggplot2::element_text(size = ggplot2::unit(12, "pt")),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(face = "bold", size = ggplot2::unit(12, "pt")),
      legend.position = "top"
    )
}

#' Make a Quantile-Quantile diagram from the output of the MCMC algorithm (DistributionAproximation)
#' and the algorithm that exhaustively calculates the distribution
#' 
#' @param estimated Outputed object of DistributionApproximation function
#' @param true Outputed object of either DistributionApproximation function or True distribution
#' computation function
#' @param filtered Make use of the classic distributuion estimation or of the filtred one
#' (number of patient taking the cocktail > beta)
#' @param color The color of the dashed line of the qq-plot
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_minimal labs
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
#' qq_plot_output(estimated = estimated_score_distribution,
#'                 true = true_score_distribution)
#'}
#' @return no returned value, should plot the quantile-quantile plot
#' of the estimated distribution (estimated) vs the true distribution (true).
#' @export
qq_plot_output <- function(estimated, true, filtered = FALSE, color = "steelblue") {
  
  requireNamespace("ggplot2")
  
  if (filtered) {
    estimated_distribution <- histogramToDitribution(estimated$Filtered_score_distribution[2:length(estimated$Filtered_score_distribution)])
    true_distribution <- histogramToDitribution(true$Filtered_score_distribution[2:length(true$Filtered_score_distribution)])
  } else {
    estimated_distribution <- histogramToDitribution(estimated$ScoreDistribution[2:length(estimated$ScoreDistribution)])
    true_distribution <- histogramToDitribution(true$ScoreDistribution[2:length(true$ScoreDistribution)])
  }
  
  num_quantiles <- min(length(estimated_distribution), length(true_distribution))
  probs <- seq(0, 1, length.out = num_quantiles)
  
  quantiles_estim <- quantile(estimated_distribution, probs)
  quantiles_true <- quantile(true_distribution, probs)
  
  qq_df <- data.frame(estimated_quantiles = quantiles_estim, true_quantiles = quantiles_true)
  
  ggplot2::ggplot(qq_df, ggplot2::aes(x = .data$estimated_quantiles, y = .data$true_quantiles)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = color) + # Adds a reference line y = x
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Estimated Distribution", y = "True Distribution", title = "QQ Plot of Estimated vs True Distribution")
}

#' Plot the histogram of the approximation of the RR distribution 
#' @param estimated The ScoreDistribution element in the list 
#' returned by the DistributionApproximation function
#' @param binwidth The width of the histogram bins
#' @param hist_color The fill color for the histogram bars
#' @param density_color The color for the density curve
#' @param sqrt A Boolean to specify whether we normalize the estimated or not, it is recommended on large random walk.
#' @param xlab Label of X axis
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_minimal after_stat
#' @importFrom dplyr data_frame
#' @importFrom rlang .data
#' @examples
#' \donttest{
#' data("ATC_Tree_UpperBound_2024")
#' data("FAERS_myopathy")
#' 
#' estimation = DistributionApproximation(epochs = 10, ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy)
#' 
#' plot_frequency(estimated = estimation$ScoreDistribution)
#'}
#' @return no returned value, should plot the histogram of the estimated distribution
#' (estimated).
#' @export
plot_frequency <- function(estimated, sqrt = FALSE, binwidth = 0.1, hist_color = "#69b3a2", density_color = "#FF5733",
                           xlab = "Score") {
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  # Create a data frame from the returned value array
  if (sqrt) {
    df <- data.frame(x = histogramToDitribution(sqrt(estimated)))
    y_lab <- "sqrt(frequency)"
  } else {
    df <- data.frame(x = histogramToDitribution(estimated))
    y_lab <- "frequency"
  }
  
  # Create histogram plot
  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(binwidth = binwidth, fill = hist_color, color = "#e9ecef") +
    ggplot2::labs(title = "The approximation of the distribution histogram", x = xlab, y = y_lab) +
    ggplot2::theme_minimal()
}
