#include "Individual.h"
// [[Rcpp::depends(RcppArmadillo)]]

Individual::Individual(const std::vector<int>& medic, double temperature) :
  medications_{medic}, temperature_{temperature}
  {}

void Individual::printMedications() const{
  for(const int& med : medications_){
    Rcpp::Rcout << med<< ' ';
  }
  Rcpp::Rcout << '\n';
}
void Individual::printTemperature() const{
  Rcpp::Rcout << temperature_ << '\n';
} 

bool Individual::matches(const std::vector<int>& observation, const std::vector<int>& upperBound) const{
  int idx;
  bool inIt;
  
  for(const int& medIdx : medications_){
    inIt = false;
    idx = 0;
    
    while (idx < observation.size() && !inIt){
      if(observation[idx] >= medIdx && observation[idx] < upperBound[medIdx]){
        inIt = true;
      }else{
        ++idx;
      }
      
    }

    // if this condition is true, then no observation matched a medication of the individual, we return false
    if(idx == observation.size())
      return false;
  }

  return true;
  
}

std::pair<double, std::pair<int,int>> Individual::computeRR(const std::vector<std::vector<int>>& medications,const Rcpp::LogicalVector& ADR,
                                 const std::vector<int>& upperBound, int RRmax, int num_thread) const{
  int yesInDelt = 0, noInDelt = 0;
  int yesNotInDelt = 0, noNotInDelt = 0;
  double sumInDelt, sumNotInDelt;
  
  //if the cocktail is empty, the related risk is zero.
  if(this->medications_.size() == 0)
    return std::make_pair(0.0, std::make_pair(0,0));
#ifdef _OPENMP
  omp_set_num_threads(num_thread);
#pragma omp parallel for reduction(+:yesInDelt) reduction(+:noInDelt) reduction(+:yesNotInDelt) reduction(+:noNotInDelt)
#endif
    for(int i = 0; i < medications.size() ; ++i){
    //does the i-th observations is included in the individual ?
    bool isInDelta = this->matches(medications[i], upperBound);
    
    if(isInDelta){
      if(ADR[i]){
        ++yesInDelt;
      }else{
        ++noInDelt;
      }
    }
    else{
      if(ADR[i]){
        ++yesNotInDelt;
      }else{
        ++noNotInDelt;
      }
    }
    
  }
  sumInDelt = yesInDelt + noInDelt;
  sumNotInDelt = yesNotInDelt + noNotInDelt;
  int returnedYesInDelt = yesInDelt;
  int returnedSumInDelt = sumInDelt;
  //during the test the denominator was frequently 0 so I make it really small if it is 0 
  //because if the denominator is 0 the numerator has to be 0 so the result would be 0 
  //no matter the denominator, we have 
  sumInDelt = (sumInDelt == 0) ? 1 : sumInDelt;
  
  double P_ADR_SEQ = static_cast<double>(yesInDelt) / static_cast<double>(sumInDelt);
  double P_ADR_NotSEQ = static_cast<double>(yesNotInDelt) / static_cast<double>(sumNotInDelt);
  
  //same as the sumInDelt
  P_ADR_NotSEQ = P_ADR_NotSEQ == 0 ? 0.00001 : P_ADR_NotSEQ;
  double RR = std::min((P_ADR_SEQ / P_ADR_NotSEQ), static_cast<double>(RRmax));
  return std::make_pair(RR, std::make_pair(returnedYesInDelt, returnedSumInDelt));
  
}

std::pair<double, std::pair<int,int>> Individual::computePHypergeom(const std::vector<std::vector<int>>& medications,
                                                        const Rcpp::LogicalVector& ADR,
                                                        const std::vector<int>& upperBound,
                                                        int ADRProportion, int notADRProportion,
                                                        int geomMax, int num_thread) const{
  int takingCocktail = 0;
  int takingCocktailHavingADR = 0;
  
  //if the cocktail is empty, the phyper is zero.
  if(this->medications_.size() == 0)
    return std::make_pair(0.0, std::make_pair(0,0));
#ifdef _OPENMP
  omp_set_num_threads(num_thread);
#pragma omp parallel for reduction(+:takingCocktail) reduction(+:takingCocktailHavingADR)
#endif
  for(int i = 0; i < medications.size() ; ++i){
    //does the i-th observations is included in the individual ?
    bool isInDelta = this->matches(medications[i], upperBound);
    if(isInDelta){
      if(ADR[i]){
        ++takingCocktailHavingADR;
      }
      ++takingCocktail;
    }
  }
  
  Rcpp::IntegerVector tmp{takingCocktailHavingADR-1};
  
  double phyper = -(Rcpp::phyper(tmp, ADRProportion, notADRProportion, takingCocktail, false, true)[0]);
  phyper = std::min(phyper, static_cast<double>(geomMax));
  
  return std::make_pair(phyper, std::make_pair(takingCocktailHavingADR, takingCocktail));
}


std::vector<std::pair<int,int>> Individual::getVertexList(const Rcpp::DataFrame& ATCtree) const{
  std::vector<std::pair<int,int>> returnedVec{0};
  
  std::vector<int> upperBound = ATCtree["upperBound"];
  std::vector<int> depth = ATCtree["ATC_length"];
  
  int idx,depthMed,nextDepth,upperBMed;
  for(const auto& med : medications_){
    idx = med;
    depthMed = depth[med];
    upperBMed = upperBound[med]-1; //shift index to cpp zero indexing
    //get the next depth
    if(depthMed == 1 || depthMed == 5)
      nextDepth = depthMed+2;
    else if(depthMed < 7)
      nextDepth = depthMed+1;
    
    //find the lower depth medications if we are not on a leaf
    if(depthMed != 7){
      
      while(idx <= upperBMed){
        //if we are on the lower depth and the medication is not on the current medications vector
        if(depth[idx] == nextDepth && (std::find(medications_.begin(),medications_.end(),idx) == medications_.end())){
          returnedVec.emplace_back(med,idx);
        }
        ++idx;
      }
    
    }
    //we come back above the medication 
    idx = med-1;
    //now find the upper depth medication (it should only have one because it is a tree), we do it if we are not on the first depth
    if(depthMed != 1){
      //the test >=0 may be deleted because we always should land on a lower depth if our current depth is not 1
      while(idx >= 0 && depthMed <= depth[idx]){
        idx--;
      }
      //we add it if it is not already on the medications vector
      if(std::find(medications_.begin(),medications_.end(),idx) == medications_.end())
        returnedVec.emplace_back(med,idx);
      
    }
    
  }

  return returnedVec;
}


bool Individual::isTrueCocktail(const std::vector<int>& upperBound) const{
  if(medications_.size() <= 1)
    return true;
  
  auto med = medications_;
  std::sort(med.begin(), med.end());
  for(int i = 0; i < med.size() - 1 ; ++i){
    //if the cocktail contain a vertice and its son, we return false 
    if(upperBound[med[i]] > med[i+1]){
      return false;
    }
  }
  // if we are here no father-son medications has been found
  return true;
}


bool Individual::operator==(const Individual& ind) const{
  if(ind.medications_.size() != medications_.size())
    return false;
  
  for(int i : medications_){
    int j =0;
    while(j < ind.medications_.size() && ind.medications_[j] != i){
      ++j;
    }
    if(j == ind.medications_.size())
      return false;
  }
  
  return true;
}

bool Individual::operator<(const Individual& ind) const{
  return true;
}


