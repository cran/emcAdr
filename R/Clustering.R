#' Clustering of the solutions of the genetic algorithm using the hclust algorithm
#' 
#' @param genetic_results The return value of the genetic algorithm
#' @param ATCtree ATC tree with upper bound of the DFS
#' @param method (from hclust function) the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param dist.normalize Do we normalize the distance (so it bellongs to [0;1])
#' @return the hierarchical clustering of the results of the genetic algorithm
#' @examples
#'\donttest{
#' data("ATC_Tree_UpperBound_2024")
#' data("FAERS_myopathy")
#' 
#' results = GeneticAlgorithm(epochs = 10, nbIndividuals = 10, 
#'             ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy)
#' 
#' hclust_genetic_solution(genetic_results = results,
#'                  ATCtree = ATC_Tree_UpperBound_2024)
#'}
#' @export
hclust_genetic_solution <- function(genetic_results,ATCtree, dist.normalize = TRUE,
                                    method = "complete"){
  requireNamespace("stats")
  
  if(dist.normalize){
    divergence <- get_dissimilarity_from_genetic_results(genetic_results, ATCtree, T)
  }else{
    divergence <- get_dissimilarity_from_genetic_results(genetic_results, ATCtree, F)
  }
  divergence <- do.call(rbind,divergence)
  divergence <- stats::as.dist(divergence)
  hc <- stats::hclust(divergence, method = method)
  
  return (hc)
}

#' Clustering of the solutions of the genetic algorithm using the hclust algorithm
#' 
#' @param genetic_results A list of cocktails in the form of integer vector
#' @param ATCtree ATC tree with upper bound of the DFS
#' @param dist.normalize Do we normalize the distance (so it belongs to [0;1])
#' @param umap_config The configuration to use in order to project the cocktails in a smaller space (umap::umap.defaults by default)
#' @return A dataframe containing UMAP 1/2 the two coordinates of each cocktails in the plane as well as the cluster number of each cocktails
#' @examples
#' \donttest{
#'  data("ATC_Tree_UpperBound_2024")
#' 
#'  results = GeneticAlgorithm(epochs = 10, nbIndividuals = 10, 
#'             ATCtree = ATC_Tree_UpperBound_2024,
#'             observations = FAERS_myopathy)
#' 
#'  hclust_genetic_solution(genetic_results = results,
#'                  ATCtree = ATC_Tree_UpperBound_2024)
#'}
#' @export
clustering_genetic_algorithm <- function(genetic_results,ATCtree,dist.normalize = TRUE,
                                         umap_config=NULL){
  requireNamespace("umap")
  requireNamespace("dbscan")

  divergence <- get_dissimilarity_from_genetic_results(genetic_results, ATCtree, dist.normalize)
  divergence <- do.call(rbind,divergence)
  divergence <- as.matrix(divergence)
  
  if(is.null(umap_config)){
    umap_config = umap::umap.defaults
  }
  
  umap_results <- umap::umap(divergence, config = umap_config)
  layout <- umap_results$layout
  dbscan_results <- dbscan::dbscan(layout)
  return (data.frame(cocktails = genetic_results,
                     UMAP1 = umap_results$layout[,1],
                     UMAP2 = umap_results$layout[,2],
                     cluster = dbscan_results$cluster))
}