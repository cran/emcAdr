#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include <iostream>
#include <algorithm>
#ifdef _OPENMP
  #include<omp.h>
#endif
#include "RcppArmadillo.h"

class Individual{
public:
  Individual() = default;
  Individual(const std::vector<int>& medic, double temperature=1);
  
  void printMedications() const;
  void printTemperature() const; 
  
  inline std::vector<int> getMedications() const{
    return medications_;
  }
  inline double getTemperature() const{
    return temperature_;
  }
  
  inline void setMedications(const std::vector<int>& newMed){
    medications_ = newMed;
  }
  inline void setTemperature(const double newTemp){
    temperature_ = newTemp;
  }
  
  bool matches(const std::vector<int>& observation, const std::vector<int>& upperBound) const;
  
  std::pair<double, std::pair<int,int>> computeRR(const std::vector<std::vector<int>>& medications, 
                                                  const Rcpp::LogicalVector& ADR,
                                                  const std::vector<int>& upperBound,
                                                  int RRmax, int num_thread=1) const;
  
  // compute the -log(phyper) given the number of people having ADR and taking 
  // this cocktail
  std::pair<double, std::pair<int,int>> computePHypergeom(const std::vector<std::vector<int>>& medications,
                                                          const Rcpp::LogicalVector& ADR,
                                                          const std::vector<int>& upperBound,
                                                          int ADRProportion, int notADRProportion,
                                                          int geomMax,
                                                          int num_thread) const;
  

  std::vector<std::pair<int,int>> getVertexList(const Rcpp::DataFrame& ATCtree) const;
  // we use this function to determine whether a cocktail is true or not
  // a true cocktail is a cocktail that does not contain a vertice and its son
  // or a vertice and its father
  bool isTrueCocktail(const std::vector<int>& upperBound) const;
  
  bool operator==(const Individual& ind) const;
  //the result is not important, we redefined this operator because of his utilization in "keepElite",
  // result does not matter since it is used in an std::pair<> and the comparaison rely on the first
  // element of the pair
  bool operator<(const Individual& ind) const;
  
private:
  std::vector<int> medications_;
  double temperature_;
};

#endif
