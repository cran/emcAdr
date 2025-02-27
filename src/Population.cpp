#include "Population.h"

// [[Rcpp::depends(RcppArmadillo)]]

Population::Population(int nbIndividuals) : individuals_{}
{
  individuals_.reserve(nbIndividuals);
}

Population::Population(const std::vector<std::vector<int>>& drugs_set) : individuals_{}
{
  individuals_.reserve(drugs_set.size());
  
  for(const auto& drugs : drugs_set){
    individuals_.emplace_back(std::make_pair(0.0, Individual(drugs)));
  }
}

Population::Population(int treeSize, int nbIndividuals,double meanMedic,
                       const std::vector<int>& upperBounds) : individuals_{}{
  std::vector<Individual> individuals = DFtoCPP_WOIndividual(treeSize, nbIndividuals, meanMedic);
  std::vector<double> fitnessVector(nbIndividuals);
  
  //we need to ensure that there is no false cocktail before we start
  std::vector<int> newMed;
  for(auto& ind : individuals){
    while(!ind.isTrueCocktail(upperBounds)){
      ind = newIndividualWithCocktailSize(treeSize, ind.getMedications().size(), 1,
                                         ind.getTemperature())[0];
    }
  }
  
  individuals_.reserve(individuals.size());
  for(int i = 0 ; i < individuals.size(); ++i){
    individuals_.push_back({fitnessVector[i], individuals[i]});
  }
}


Population::Population(const Population& pop) : individuals_{pop.individuals_}
{}

void Population::addAnIndividualToPopulation(const std::pair<double,Individual>& ind){
  individuals_.push_back(ind);
}

void Population::evaluate(const std::vector<std::vector<int>>& medications,
                          const Rcpp::LogicalVector& ADR,
                          const std::vector<int>& upperBound, int ADR_proportion,
                          int not_ADR_proportion, int geom_max, int num_thread){
  for(auto& indiv : individuals_){
    indiv.first = indiv.second.computePHypergeom(medications, ADR, upperBound,
                                                 ADR_proportion, not_ADR_proportion,
                                                 geom_max, num_thread).first;
  }
}

void Population::keepElite(int nbElite, Population& matingPool) const{
  if (nbElite > 0){
    std::priority_queue<std::pair<double,Individual>, std::vector<std::pair<double,Individual>>,
                        std::greater<std::pair<double,Individual>>> pq;
    
    for (const auto & i : individuals_) {
      if (pq.size() < nbElite) {
        pq.push(i);
      } else if (i > pq.top()) {
        pq.pop();
        pq.push(i);
      }
    }
    
    while (!pq.empty()) {
      matingPool.addAnIndividualToPopulation(pq.top());
      pq.pop();
    }
  }
}

//we may want to return the mating pool instead of passing it by reference in parameter
void Population::tournamentSelection(int tournamentSize, Population& matingPool, int nbDrawing) const{
  int popSize = individuals_.size();
  int indexDrawing;
  double indexMax;
  std::vector<int> selectedIndividuals;
  selectedIndividuals.reserve(tournamentSize);
  //we have to add nbDrawing in the mating pool
  for(int i = 0 ; i < nbDrawing; ++i){
    //select the index of the players of the tournament
    selectedIndividuals.clear();
    for(int j = 0; j < tournamentSize ; ++j){
      do{
        indexDrawing = Rcpp::runif(1, 0, popSize)[0];
      } while (std::find(selectedIndividuals.begin(), selectedIndividuals.end(), indexDrawing) != selectedIndividuals.end());
      selectedIndividuals.push_back(indexDrawing);
    }
    //now make them compete again each other
    indexMax = selectedIndividuals[0];
    for(int j = 1; j < tournamentSize ; ++j){
      indexMax = individuals_[indexMax] > individuals_[selectedIndividuals[j]] ? indexMax : selectedIndividuals[j];
    }
    matingPool.addAnIndividualToPopulation(individuals_[indexMax]);
  }
}

void Population::crossover(int nbElite, const std::vector<int>& ATClength, const std::vector<int>& upperBounds,
                           const Rcpp::DataFrame& ATCtree, double p_crossover){
  int remainingIndividuals = individuals_.size() - nbElite;
  int selectedNode;
  int upperBound;
  Individual tmp1, tmp2;
  double draw;
  
  for(int i = nbElite ; i < individuals_.size()-2; i+=2){
    draw = Rcpp::runif(1,0,1)[0];
    // With probability p_crossover we opperate a crossover, otherwise we let the individuals as it is
    if(draw <= p_crossover){
      do{
        
        do{
          selectedNode = trunc(Rcpp::runif(1,0,ATCtree.nrow())[0]);
          } while (ATClength[selectedNode] == 7);
        
        upperBound = upperBounds[selectedNode];
        
        tmp1 = crossoverMutation(individuals_[i].second, individuals_[i+1].second,
                                 ATCtree, selectedNode, upperBound);

        tmp2 = crossoverMutation(individuals_[i+1].second, individuals_[i].second,
                                 ATCtree, selectedNode, upperBound);

        }while (!tmp1.isTrueCocktail(upperBounds) || !tmp2.isTrueCocktail(upperBounds));
      individuals_[i].second = tmp1;
      individuals_[i+1].second = tmp2;
    }
  }
  
  //if the number of remaining individuals is even, we have to operate one more crossover on the last 2 individuals
  //otherwise we just have one individual so we just copy it (equivalent of doing nothing)
  if(remainingIndividuals %2 == 0){
    draw = Rcpp::runif(1,0,1)[0];
    if(draw <= p_crossover){
      do{
        do{
          selectedNode = trunc(Rcpp::runif(1,0,ATCtree.nrow())[0]);
        } while (ATClength[selectedNode] == 7);
        upperBound = upperBounds[selectedNode];
        
        tmp1 = crossoverMutation(individuals_[individuals_.size()-2].second, 
                                 individuals_[individuals_.size()-1].second, ATCtree,
                                 selectedNode, upperBound);

        tmp2 = crossoverMutation(individuals_[individuals_.size()-1].second, 
                                 individuals_[individuals_.size()-2].second,
                                 ATCtree, selectedNode, upperBound);

      }while (!tmp1.isTrueCocktail(upperBounds) || !tmp2.isTrueCocktail(upperBounds));
      individuals_[individuals_.size()-2].second = tmp1;
      individuals_[individuals_.size()-1].second = tmp2;
    }
  }
  
}

void Population::mutate(int nbElite, double p_mutation, const Rcpp::DataFrame& ATCtree,
                        const std::vector<int>& upperBounds, double alpha){
  double draw, drawMutation, chosenVertexIdx;
  bool emptyCocktail;
  std::vector<std::pair<int,int>> vertex;
  Individual tmp;
  
  for(int i = nbElite ;  i < individuals_.size() ; ++i){
      draw = Rcpp::runif(1,0,1)[0];
      //for each individual in the population we draw a number uniformly in [0;1]
      // with probability p_mutation we mutate the individual i
    if(draw <= p_mutation){
      //we mutate with type 1 mutation or type 2 equiprobably

        emptyCocktail = (individuals_[i].second.getMedications().size() == 0);
        drawMutation = Rcpp::runif(1,0,1)[0];
        if(drawMutation <= 0.5){
          //type1 mutation
          tmp = type1Mutation(individuals_[i].second, ATCtree.nrow(), alpha, emptyCocktail);
        }
        else{
          //type2 mutation
          if(!emptyCocktail){
            vertex = individuals_[i].second.getVertexList(ATCtree);
            chosenVertexIdx = trunc(Rcpp::runif(1,0,vertex.size())[0]);
            tmp = type2Mutation(individuals_[i].second, ATCtree.nrow(), vertex[chosenVertexIdx]);
            
          }
        }
      if(tmp.isTrueCocktail(upperBounds))
        individuals_[i].second = tmp;
    }
  }
}

std::pair<IntMatrix, std::vector<int>> Population::pretraitement(const std::vector<int>& depth,
                                                                 const std::vector<int>& father) const{
  IntMatrix M;
  int number_of_medication = std::accumulate(individuals_.begin(), individuals_.end(), 0,
                  [](int a, const std::pair<double,Individual>& p){
                    return a + p.second.getMedications().size();
                  });
  M.resize(number_of_medication);
  for(auto& line : M){
    line.reserve(6); // depth of the tree
  }
  
  std::vector<int> idx;
  idx.reserve(individuals_.size());
  int index = 0;
  int medline = 0;
  for(int i = 0; i < individuals_.size(); ++i){
    for(const auto& med : individuals_[i].second.getMedications()){
      int pred = med;
      //if we are not on a leaf 
      for(int j = 0; j < 5-depth[med]; ++j){
        M[medline].push_back(med);
      }
      //then we use the regular treatment
      while(pred != -1){
        M[medline].push_back(pred);
        pred = father[pred];
      }
      M[medline].push_back(-1);
      ++medline;
    }
    idx.push_back(index);
    index+= individuals_[i].second.getMedications().size();
  }
  
  return {M, idx};
}

RealMatrix Population::initSimilarityMatrix() const{
  RealMatrix sim;
  sim.resize(individuals_.size());
  for(auto& vec : sim){
    vec.resize(individuals_.size(),-1);
  }
  
  return sim;
}
//refactor with dist
double Population::dist_norm(int i, int j, const IntMatrix& M, const std::vector<int>& idx) const{
  //if j == idx.size -1 the last index of a medication of C_j is M.size()-1 (last row of M)
  int max_idx_j = j == idx.size()-1 ? M.size() : idx[j+1];
  int max_idx_i = idx[i+1];
  int min_idx_i = idx[i], min_idx_j = idx[j];
  double ATC_height = 5;

  std::vector<int> indexC1(max_idx_i - min_idx_i);
  std::iota(indexC1.begin(), indexC1.end(), min_idx_i);
  std::vector<int> indexC2(max_idx_j - min_idx_j);
  std::iota(indexC2.begin(), indexC2.end(), min_idx_j);
  
  std::vector<int> deleteC1;
  deleteC1.reserve(indexC1.size());
  int depth = 0;
  double cost = 0;
  int initial_length = indexC1.size() + indexC2.size();
  
  while(!indexC1.empty() && !indexC2.empty()){
    for(const auto& idx1 : indexC1){
      for(const auto& idx2 : indexC2){
        if(M[idx1][depth] == M[idx2][depth]){
          //if exactly the same node
          if(M[idx1] == M[idx2]){
            cost+=0;
          }else{
            cost+= depth - 
              std::min(std::count(M[idx1].begin(),M[idx1].end(), M[idx1][0]),
                       std::count(M[idx2].begin(),M[idx2].end(), M[idx2][0])) +
                         1;
          }
          deleteC1.push_back(idx1);
          indexC2.erase(std::remove(indexC2.begin(), indexC2.end(), idx2),
                        indexC2.end());
          break;
        }
      }
    }
    indexC1.erase(std::remove_if(indexC1.begin(), indexC1.end(), 
                                 [&deleteC1](int value){
                                   return std::find(deleteC1.begin(), deleteC1.end(), value) != deleteC1.end();
                                 }), indexC1.end());
    deleteC1.clear();
    ++depth;
  }
  
  // we add to the cost the cost of adding remaining drugs, the cost is :
  //number of drugs to add times cost of adding which is (ATC_tree height / 2)
  // and the ATC height is 5
  double insertion_cost = (ATC_height/2.0);
  cost += (indexC1.size() + indexC2.size()) * insertion_cost; 
  
  return cost / (static_cast<double>(initial_length) * insertion_cost);
}

double Population::dist(int i, int j, const IntMatrix& M, const std::vector<int>& idx) const{
  //if j == idx.size -1 the last index of a medication of C_j is M.size()-1 (last row of M)
  int max_idx_j = j == idx.size()-1 ? M.size() : idx[j+1];
  int max_idx_i = idx[i+1];
  int min_idx_i = idx[i], min_idx_j = idx[j];
  double ATC_height = 5;
  
  std::vector<int> indexC1(max_idx_i - min_idx_i);
  std::iota(indexC1.begin(), indexC1.end(), min_idx_i);
  std::vector<int> indexC2(max_idx_j - min_idx_j);
  std::iota(indexC2.begin(), indexC2.end(), min_idx_j);
  
  std::vector<int> deleteC1;
  deleteC1.reserve(indexC1.size());
  int depth = 0;
  double cost = 0;
  
  while(!indexC1.empty() && !indexC2.empty()){
    for(const auto& idx1 : indexC1){
      for(const auto& idx2 : indexC2){
        if(M[idx1][depth] == M[idx2][depth]){
          //allow to compare cocktail that are not drug cocktails
          if(M[idx1] == M[idx2]){
            cost+=0;
          }else{
            cost+= depth - 
              std::min(std::count(M[idx1].begin(),M[idx1].end(), M[idx1][0]),
                       std::count(M[idx2].begin(),M[idx2].end(), M[idx2][0])) +
                         1;
          }
          deleteC1.push_back(idx1);
          //cost+= depth;
          indexC2.erase(std::remove(indexC2.begin(), indexC2.end(), idx2),
                        indexC2.end());
          break;
        } 
      }
    }
    indexC1.erase(std::remove_if(indexC1.begin(), indexC1.end(), 
                                 [&deleteC1](int value){
                                   return std::find(deleteC1.begin(), deleteC1.end(), value) != deleteC1.end();
                                 }), indexC1.end());
    deleteC1.clear();
    ++depth;
  }
  
  // we add to the cost the cost of adding remaining drugs, the cost is :
  //number of drugs to add times cost of adding which is (ATC_tree height / 2)
  // and the ATC height is 5
  double insertion_cost = (ATC_height/2.0);
  cost += (indexC1.size() + indexC2.size()) * insertion_cost; 
  
  return cost;
}

RealMatrix Population::similarity(const IntMatrix& M, const std::vector<int>& idx) const{
  RealMatrix S = initSimilarityMatrix();

  for(int i = 0 ; i < idx.size()-1; ++i){
    S[i][i] = 1; // a cocktail have a perfect similarity with himself
    for(int j = i+1; j < idx.size(); ++j){
      if(S[i][j] == -1){
        //distance from cocktail i to cocktail j is the same as the distance of
        //cocktail j to cocktail i
        
        double sim = 1-dist_norm(i,j,M,idx);
        
        S[i][j] = sim;
        S[j][i] = sim;
      }
    }
  }
  S[idx.size()-1][idx.size()-1] = 1;
  return S;
}

RealMatrix Population::dissimilarity(const IntMatrix& M,
                                     const std::vector<int>& idx,
                                     bool normalize) const{
  RealMatrix D = initSimilarityMatrix();
  double dissim;
  for(int i = 0; i < idx.size()-1; ++i){
    D[i][i] = 0;
    for(int j = i+1; j < idx.size(); ++j){
      if(D[i][j] == -1){
        
        if(normalize){
          dissim = dist_norm(i,j,M,idx);
        }else{
          dissim = dist(i,j,M,idx);
        }
        D[i][j] = dissim;
        D[j][i] = dissim;
      }
    }
  }
  
  D[idx.size()-1][idx.size()-1] = 0;
  return D;
}


void Population::penalize(const std::vector<int>& depth, const std::vector<int>& father){
  IntMatrix M;
  std::vector<int> indexM;

  std::tie(M, indexM) = pretraitement(depth,father);

  RealMatrix S = similarity(M, indexM);
  for(int i = 0 ; i < individuals_.size(); ++i){
    individuals_[i].first = (individuals_[i].first / std::accumulate(S[i].begin(), S[i].end(), 0.0));
  }
}

void Population::clear(){
  individuals_.clear();
}

double Population::getMean() const{
  const auto size = static_cast<double>(individuals_.size());
  return std::accumulate(individuals_.begin(), individuals_.end(),
                         0.0, [](double& a, const std::pair<double,Individual>& b){
                           return a + b.first;
                           }) / size;
}

int Population::bestIndividual() const{
  int i_max = 0;
  for(int i = 1; i < individuals_.size() ; ++i){
    i_max = individuals_[i].first > individuals_[i_max].first ? i : i_max;
  }
  return i_max;
}

std::vector<std::vector<int>> Population::getMedications() const{
  std::vector<std::vector<int>> ret;
  ret.reserve(individuals_.size());
  
  for(const auto& ind : individuals_){
    ret.push_back(ind.second.getMedications());
  }
  
  return ret;
}

std::vector<double> Population::getRR() const{
  std::vector<double> ret;
  ret.reserve(individuals_.size());
  
  for(const auto& ind : individuals_)
    ret.push_back(ind.first);
  
  return ret;
}

void Population::printPopulation(std::ostream& ost) const{
  for(const auto& indiv : individuals_){
    ost << "RR : " << indiv.first << "\n medication : ";
    indiv.second.printMedications();
  }
}

void Population::printSummary(int epoch, double populationMean, int populationBestIndex) const{
  Rcpp::Rcout << "epoch : " << epoch << " | mean : " << populationMean << " | best score : ";
  Rcpp::Rcout << individuals_[populationBestIndex].first << " | best cocktail : ";
  individuals_[populationBestIndex].second.printMedications();
}

Population& Population::operator=(const Population& pop){
  individuals_ = pop.individuals_;
  return *this;
}

std::ostream& operator<<(std::ostream& ost, const Population& pop){
  pop.printPopulation(ost);
  return ost;
}

