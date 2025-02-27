## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(emcAdr)
library(gridExtra)
library(ggplot2)

## -----------------------------------------------------------------------------
data("ATC_Tree_UpperBound_2024") ## The ATC tree containing the upper bound
head(ATC_Tree_UpperBound_2024)

## -----------------------------------------------------------------------------
data("FAERS_myopathy") ## The Individual Case Safety Reports in the following format
head(FAERS_myopathy) ## First column is "patientATC" containing vector of index of drug intake for patient at row i in the ATC tree.
## Second column is "patientADR", value is true if the patient at row i experience the considered AE and false otherwise.

## -----------------------------------------------------------------------------
estimated_distribution_size1_300ksteps_t500 <- DistributionApproximation(10,
                        ATC_Tree_UpperBound_2024, FAERS_myopathy,
                        temperature = 500, nbResults = 200,
                        Smax = 1,num_thread = 8)

## -----------------------------------------------------------------------------
true_distribution_size1 <- trueDistributionDrugs(ATC_Tree_UpperBound_2024, FAERS_myopathy,
                                beta = 4, num_thread = 8)

## -----------------------------------------------------------------------------
qq_plot_output(estimated_distribution_size1_300ksteps_t500,
               true_distribution_size1, filtered = T)

## -----------------------------------------------------------------------------
plot_estimated_distribution <- plot_frequency(
    estimated_distribution_size1_300ksteps_t500$Filtered_score_distribution[2:length(estimated_distribution_size1_300ksteps_t500$Filtered_score_distribution)], binwidth = .3, sqrt = F, xlab = "H(C)") + labs(title = "Estimated distribution of risks among size 1 cocktails") + theme(plot.title = element_text(size=20)) + ylim(c(0,0.35))

plot_true_distribution <- plot_frequency(
    true_distribution_size1$Filtered_score_distribution[2:length(true_distribution_size1$Filtered_score_distribution)], binwidth = .3, xlab = "H(C)") + labs(title = "True Distribution of Risks among size 1 cocktails") + theme(plot.title = element_text(size=20)) + ylim(c(0,0.35))

grid.arrange(plot_estimated_distribution, plot_true_distribution ,nrow = 2)

## -----------------------------------------------------------------------------
genetic_algorithm_results <- GeneticAlgorithm(epochs = 20, nbIndividuals = 100,
                                              ATC_Tree_UpperBound_2024, FAERS_myopathy,
                                              num_thread = 2, diversity = T,
                                              p_mutation = 0.2, alpha = 1.5)

## -----------------------------------------------------------------------------
## We put the estimation of risk distribution of each cocktails size in a list
distribution_list <- list(estimated_distribution_size1_300ksteps_t500)

p_values <- p_value_genetic_results(distribution_list, genetic_algorithm_results)
p_values

