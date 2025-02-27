#ifndef POPULATION_H
#define POPULATION_H

#include "MCMC.h"
#include <queue>

using RealMatrix = std::vector<std::vector<double>>;
using IntMatrix = std::vector<std::vector<int>>;

class Population{
public:
  Population() = default;
  Population(int nbIndividuals);
  Population(const std::vector<std::vector<int>>&);
  Population(int treeSize, int nbIndividuals,double meanMedic,
             const std::vector<int>& upperBounds);
  Population(const Population& pop);
  
  inline std::vector<std::pair<double,Individual>> getIndividuals() const{
    return individuals_;
  }
  
  inline void setIndividuals(const std::vector<std::pair<double,Individual>>& newPopulation){
    individuals_ = newPopulation;
  }
  
  void addAnIndividualToPopulation(const std::pair<double,Individual>& ind);
  
  void evaluate(const std::vector<std::vector<int>>& medications,
                const Rcpp::LogicalVector& ADR,
                const std::vector<int>& upperBound, int ADR_proportion,
                int not_ADR_proportion, int geom_max, int num_thread);
  
  void keepElite(int nbElite, Population& matingPool) const;
  
  void tournamentSelection(int tournamentSize, Population& matingPool, int nbDrawing) const;
  
  void crossover(int nbElite, const std::vector<int>& ATClength, const std::vector<int>& upperBounds,
                 const Rcpp::DataFrame& ATCtree, double p_crossover);
  
  void mutate(int nbElite, double p_mutation, const Rcpp::DataFrame& ATCtree,
              const std::vector<int>& upperBounds, double alpha); 
  
  RealMatrix initSimilarityMatrix() const;
  
  std::pair<IntMatrix, std::vector<int>> pretraitement(const std::vector<int>& depth,
                                                       const std::vector<int>& father) const;
  
  double dist_norm(int i, int j, const IntMatrix& M, const std::vector<int>& idx) const;
  
  double dist(int i, int j, const IntMatrix& M, const std::vector<int>& idx) const;
  
  RealMatrix similarity(const IntMatrix& M, const std::vector<int>& idx) const;
  
  RealMatrix dissimilarity(const IntMatrix& M, const std::vector<int>& idx,
                           bool normalize = true) const;
  
  void penalize(const std::vector<int>& depth, const std::vector<int>& father);
  
  void clear();
  
  double getMean() const;
  
  int bestIndividual() const;
  
  std::vector<std::vector<int>> getMedications() const;
  
  std::vector<double> getRR() const;
  
  void printPopulation(std::ostream& ost) const;
  
  void printSummary(int epoch, double populationMean, int populationBestIndex) const;
  
  Population& operator=(const Population& pop);
  
  friend std::ostream& operator<<(std::ostream& ost, const Population& pop);
  
private:
  std::vector<std::pair<double,Individual>> individuals_;
};

#endif
