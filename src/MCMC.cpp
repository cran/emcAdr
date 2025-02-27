#include "MCMC.h"

// [[Rcpp::depends(RcppArmadillo)]]


double meanMedications(const std::vector<std::vector<int>>& observations){
  double final_sum = std::accumulate(
    observations.begin(), observations.end(), 0.0,
    [](double sum, const std::vector<int>& innerVector) {
      return sum + innerVector.size();
    }
  );
  
  return observations.size() > 0 ? final_sum / static_cast<double>(observations.size()):
    0;
}


std::pair<Individual,double> largerScore(const std::pair<Individual,double>& firstScore, const std::pair<Individual,double>& secScore,
                                      const std::pair<Individual,double>& thirdScore, const std::pair<Individual,double>& fourthScore){
  if(firstScore.second >= secScore.second && firstScore.second >= thirdScore.second && firstScore.second >= fourthScore.second)
    return firstScore;
  else if(secScore.second >= firstScore.second && secScore.second >= thirdScore.second && secScore.second >= fourthScore.second)
    return secScore;
  else if(thirdScore.second >= firstScore.second && thirdScore.second >= secScore.second &&thirdScore.second >= fourthScore.second)
    return thirdScore;
  else
    return fourthScore;
}

void addScoretoDistribution(const double score,std::vector<unsigned int>& vec){
  double dIndex = (score) * 10;
  int index = static_cast<int>(dIndex);
  vec[index]++;
}

double addToBestCocktails(std::vector<std::pair<Individual,double>>& bestResults,
                        const std::pair<Individual,double>& currentResult,
                        int nbResults, double minScore,const std::vector<int>& upperBound){
  double newMinScore = minScore;
  if(isNotInResultList(bestResults,currentResult) && currentResult.first.isTrueCocktail(upperBound)){
    if(bestResults.size() < nbResults){
      bestResults.emplace_back(currentResult);
      newMinScore = minScore < currentResult.second ? minScore : currentResult.second;
    }
    else if(minScore < currentResult.second){
      auto it = std::find_if(bestResults.begin(),bestResults.end(),
                             [minScore](const std::pair<Individual,double>& p){return p.second == minScore;});
      if(it != bestResults.end()){
        bestResults.erase(it);
      }
      bestResults.emplace_back(currentResult);
      auto tmpMin = *std::min_element(bestResults.begin(),bestResults.end(),
                                      [](const std::pair<Individual,double>& lp,const std::pair<Individual,double>& rp){
                                        return lp.second < rp.second; 
                                      });
      newMinScore = tmpMin.second;
    }
  }
  return newMinScore;
}

void addPairToSet(const Individual& i, std::set<std::pair<int,int>>& p){
  int minPair,maxPair;
  
  std::vector<int> mutIndivTmp = i.getMedications();
  if(mutIndivTmp.size() == 2){
    if(mutIndivTmp[0] > mutIndivTmp[1]){
      maxPair = mutIndivTmp[0];
      minPair = mutIndivTmp[1];
    }else{
      maxPair = mutIndivTmp[1];
      minPair = mutIndivTmp[0];
    }
    p.insert(std::make_pair(minPair, maxPair));
  }
}

bool isNotInResultList(const std::vector<std::pair<Individual,double>>& bestResults,
                       const std::pair<Individual,double>& bestResult){
  return std::find(bestResults.begin(),bestResults.end(),bestResult) == bestResults.end();
}


bool mutatedByType2(std::vector<int> X, std::vector<int> Y, const std::vector<std::pair<int,int>>& vertex){
  std::sort(X.begin(),X.end());
  std::sort(Y.begin(),Y.end());
  std::vector<int> diff{};
  diff.reserve(std::max(X.size(),Y.size()));
  std::set_symmetric_difference(X.begin(),X.end(),Y.begin(),Y.end(),
                                std::back_inserter(diff));
  diff.shrink_to_fit();
  
  std::pair<int,int> find{diff[0],diff[1]};
  for(const auto& i : vertex){
    if(PermutEqual(find,i)){
      return true;
    }
  }
  return false;
}


std::vector<Individual> DFtoCPP_WOtemp(const Rcpp::List& startingInd){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(startingInd.size());
  Rcpp::NumericVector temp;
  int itemp = 1;
  
  for(const std::vector<int> ind : startingInd){
    temp = Rcpp::runif(1,itemp-1,itemp);
    returnedVec.push_back(Individual{ind,temp[0]});
    //because the temperatures have to be increasing
    ++itemp;
  }
  return returnedVec;
}

std::vector<Individual> DFtoCPP_WOIndividual(int treeSize, int nbIndividuals,double meanMedic, const Rcpp::NumericVector& temperatures){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(nbIndividuals);
  int nbMedic,num;
  std::vector<int> medicVec(0);
  
  for(int i = 0; i < nbIndividuals; ++i){
    nbMedic = 1 + Rcpp::rpois(1,meanMedic)[0];
    medicVec.reserve(nbMedic);
    
    for(int j = 0; j < nbMedic ; ++j){
      //Draw every medications for a patient -> consider using sample ?
      num = trunc(Rcpp::runif(1,0,treeSize)[0]);
      num = num == treeSize ? treeSize-1 : num;
      medicVec.push_back(num);
    }
    returnedVec.push_back(Individual{medicVec,temperatures[i]});
    medicVec.clear();
  }
  return returnedVec;
}

std::vector<Individual> DFtoCPP_WOIndividual(int treeSize, int nbIndividuals,double meanMedic){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(nbIndividuals);
  int nbMedic,num;
  std::vector<int> medicVec(0);
  
  for(int i = 0; i < nbIndividuals; ++i){
    nbMedic = 1 + Rcpp::rpois(1,meanMedic)[0];
    medicVec.reserve(nbMedic);
    
    for(int j = 0; j < nbMedic ; ++j){
      //Draw every medications for a patient -> consider using sample ?
      num = trunc(Rcpp::runif(1,0,treeSize)[0]);
      num = num == treeSize ? treeSize-1 : num;
      medicVec.push_back(num);
    }
    returnedVec.push_back(Individual{medicVec, 1});
    medicVec.clear();
  }
  return returnedVec;
}

std::vector<Individual> DFtoCPP_Wtemp(const Rcpp::List& startingInd,const Rcpp::NumericVector& startingTemp){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(startingInd.size());

  for(int i = 0 ; i < startingInd.length(); ++i){
    returnedVec.push_back(Individual{startingInd[i],startingTemp[i]});
  }
  return returnedVec;
  
}

std::vector<Individual> DFtoCPP_WOtempAndIndividual(int treeSize, int nbIndividuals,double meanMedic){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(nbIndividuals);
  int nbMedic,num, itemp = 1;
  double temp;
  std::vector<int> medicVec(0);
  
  for(int i = 0; i < nbIndividuals; ++i){
    nbMedic = 1 + Rcpp::rpois(1,meanMedic)[0];
    medicVec.reserve(nbMedic);
    
    for(int j = 0; j < nbMedic ; ++j){
      //Draw every medications for a patient -> consider using sample ?
      num = trunc(Rcpp::runif(1,0,treeSize)[0]);
      num = num == treeSize ? treeSize-1 : num;
      medicVec.push_back(num);
    }
    temp = Rcpp::runif(1,itemp-1,itemp)[0];
    
    returnedVec.push_back(Individual{medicVec,temp});
    ++itemp;  
    medicVec.clear();
  }
  return returnedVec;
}

std::vector<Individual> newIndividualWithCocktailSize(int treeSize, int cocktailSize, int nbIndividuals, double temperature){
  std::vector<Individual> returnedVec;
  returnedVec.reserve(nbIndividuals);
  int num;
  
  std::vector<int> medicVec;
  for(int i = 0 ; i < nbIndividuals; ++i){
    medicVec.reserve(cocktailSize);
    for(int j = 0; j < cocktailSize ; ++j){
      //Draw every medications for a patient 
      do
      {
        num = trunc(Rcpp::runif(1,0,treeSize)[0]);
      } while(std::find(medicVec.begin(), medicVec.end(), num) != std::end(medicVec));
      medicVec.push_back(num);

    }
    returnedVec.push_back(Individual{medicVec,temperature});
    medicVec.clear();
  }
  
  return returnedVec;
}


std::set<std::pair<int,int>> getADRPairs(const Rcpp::List& observationsMed, const Rcpp::LogicalVector& ADR){
  std::set<std::pair<int,int>> retSet;
  int i = 0;
  int max, min;
  for(const std::vector<int> tab : observationsMed){
    
    if(tab.size() == 2 && ADR[i]){
      if(tab[0] > tab[1]){
        max = tab[0];
        min = tab[1];
      }
      else{
        max = tab[1];
        min = tab[0];
      }
      retSet.insert(std::make_pair(min,max));
    }
    ++i;
  }

  return retSet;
}

Individual type1Mutation(const Individual& indiv, int treeSize, double alpha, bool emptyCocktail){
  //peut optimiser les appels à getMedication()
  int mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
  mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
  std::vector<int> newMed = indiv.getMedications();
  
  if(!emptyCocktail){
    
    double addAcceptation = (alpha/indiv.getMedications().size());
    double draw = Rcpp::runif(1,0,1)[0];
    
    if(addAcceptation >= draw){
      while (std::find(newMed.begin(), newMed.end(), mutateMed) != std::end(newMed)){
        // we draw a medication that is not inside the cocktail for the moment
        mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
        mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
      }
      
      newMed.push_back(mutateMed);
    }
    else{ // here we remove an element from the medications of the individuals
      int removeIndex = trunc(Rcpp::runif(1,0,indiv.getMedications().size())[0]);
      removeIndex = removeIndex == indiv.getMedications().size() ? removeIndex-1 : removeIndex;
      newMed.erase(newMed.begin() + removeIndex);
    }
    
  }
  else{ // here the cocktail is empty so we add a drug with probability 1 
    newMed.push_back(mutateMed);
  }
  
  return {newMed,indiv.getTemperature()};
}

Individual adjustedType1Mutation(const Individual& indiv, int treeSize, double alpha, bool emptyCocktail){
  int mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
  mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
  std::vector<int> newMed = indiv.getMedications();
  
  if(emptyCocktail){
    newMed.push_back(mutateMed);
  }
  else{
    
    double draw = Rcpp::runif(1,0,1)[0];
    //here alpha is a probability 
    
    if(draw <= alpha){
      //here we add since alpha is the probability to add a drug to the cocktail
      while (std::find(newMed.begin(), newMed.end(), mutateMed) != std::end(newMed)){
        // we draw a medication that is not inside the cocktail for the moment
        mutateMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
        mutateMed = mutateMed == treeSize ? treeSize-1 : mutateMed;
      }
      
      newMed.push_back(mutateMed);
    }
    else{
      //the modified medication
      int chosenMedic = Rcpp::sample(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(newMed)),1)[0];
      auto end = std::remove(newMed.begin(),newMed.end(), chosenMedic);
      newMed.erase(end, newMed.end());
      int replacementMed;
      do{
        replacementMed = trunc(Rcpp::runif(1,0,treeSize)[0]);
        replacementMed = replacementMed == treeSize ? treeSize-1 : replacementMed;
      }while(std::find(newMed.begin(), newMed.end(), replacementMed) != std::end(newMed));
      // here replacementMed is either a new medication or the chosen medic, we add it only if it is not the chosen medic
      if( replacementMed != chosenMedic){
        newMed.push_back(replacementMed);
      }
      
    }
  }
  
  return {newMed, indiv.getTemperature()};
}

Individual type2Mutation(const Individual& indiv, int treeSize, const std::pair<int,int>& p){
  // if a patient got every medication
  if(indiv.getMedications().size() == treeSize)
    return indiv;
  
  std::vector<int> prevMedic = indiv.getMedications();
  auto rmEnd = std::remove(prevMedic.begin(),prevMedic.end(),p.first);
  prevMedic.erase(rmEnd,prevMedic.end());
  prevMedic.push_back(p.second);
  
  return {prevMedic,indiv.getTemperature()};
}

Individual crossoverMutation(const Individual& indiv1, const Individual& indiv2,const Rcpp::DataFrame& ATCtree,
                             int selectedNode, int upperBound){
  std::vector<int> newMedi{};
  newMedi.reserve(indiv1.getMedications().size() + indiv2.getMedications().size());
  
  for(int med : indiv1.getMedications()){
    if(med < selectedNode || med >= upperBound){
      newMedi.push_back(med);
    }
  }
  for(int med : indiv2.getMedications()){
    if(med >= selectedNode && med < upperBound){
      newMedi.push_back(med);
    }
  }
  newMedi.shrink_to_fit();
  
  return {newMedi, indiv1.getTemperature()};
}

std::pair<std::vector<int>,std::vector<int>> treeDepthFather(const std::vector<int>& ATC_length){
 std::vector<int> depth, pere, pere_crt;
 depth.reserve(ATC_length.size());
 pere.reserve(ATC_length.size());
 pere_crt.resize(5); // we need the father for each depth
 pere_crt[0] = -1; // father of depth 1 node is -1
 
 for(int i = 0; i < ATC_length.size(); ++i){
   if(ATC_length[i] > 1 && ATC_length[i] < 7){
     depth.push_back(ATC_length[i] -1);
   }
   else if(ATC_length[i] == 7){
     depth.push_back(5);
   }
   else{
     depth.push_back(1);
   }
   
   pere.push_back(pere_crt[depth[i]-1]);
   if(depth[i] != 5){
     pere_crt[depth[i]] = i;
   }
 }
 return {depth, pere};
}


