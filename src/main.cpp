// we only include MCMC.h which pulls RcppArmadillo.h and Rcpp.h in for us
#include "MCMC.h"
#include "Population.h"
#include <iostream>
#include <string_view>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <set>

#ifdef _OPENMP
  #include <omp.h>
#endif

using Rcpp::DataFrame;
// [[Rcpp::plugins(openmp)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]


//'The MCMC method that runs the random walk on a single cocktail in order to estimate the distribution of score among cocktails of size Smax.
//'
//'@param epochs : number of steps for the MCMC algorithm
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root, also see on the github repo for an example)
//'@param observations : real observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'
//'@param temperature : starting temperature, default = 1 (denoted T in the article)
//'@param nbResults : Number of returned solution (Cocktail of size Smax with the best oberved score during the run), 5 by default
//'@param Smax : Size of the cocktail we approximate the distribution from
//'@param p_type1 : probability to operate type1 mutation. Note :
//'the probability to operate the type 2 mutation is then 1 - P_type1. P_type1 must be in [0;1]. Default is .01
//'@param beta : filter the minimum number of patients that must have taken the 
//'cocktail for his risk to be taken into account in the DistributionScoreBeta default is 4
//'@param max_score : maximum number the score can take. Score greater than this 
//'one would be added to the distribution as the value max_score. Default is 500
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'@param verbose : Output summary (default is false)
//'
//'@return I no problem, return a List containing :
//' - ScoreDistribution : the distribution of the score as an array with each cells
//' representing the number of risks =  (index-1)/ 10
//' - Outstanding_score : An array of the score greater than max_score,
//' - Best_cocktails : the nbResults bests cocktails encountered during the run.
//' - Best_scores : Score corresponding to the bestCocktails.
//' - FilteredDistribution : Distribution containing score for cocktails taken by at
//' least beta patients.
//' - Best_cocktails_beta : the nbResults bests cocktails taken by at least beta patients
//' encountered during the run.
//' - Best_scores_beta : Score corresponding to the bestCocktailsBeta.
//' - cocktailSize : Smax parameter used during the run.
//'; Otherwise the list is empty
//'
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' estimation = DistributionApproximation(epochs = 10, ATCtree = ATC_Tree_UpperBound_2024,
//'             observations = FAERS_myopathy)
//'}
//'@export
//[[Rcpp::export]]
Rcpp::List DistributionApproximation(int epochs, const DataFrame& ATCtree, const DataFrame& observations,
                                              int temperature = 1, int nbResults = 5, int Smax = 2,
                                              double p_type1 = .01, int beta = 4, int max_score = 500,
                                              int num_thread = 1, bool verbose = false){
  //arguments verification
  if(p_type1 > 1 || p_type1 < 0 || epochs < 1){
    Rcpp::Rcerr << "problem in the values of the parameter in the call of this function \n";
    return Rcpp::List();
  }
  
  // OMP SET NUM THREAD = k, s.t. 1 <= k <= omp_get_num_procs()
#ifdef _OPENMP
  if(num_thread < 1 || num_thread > omp_get_num_procs()){
    Rcpp::Rcerr << "Wrong thread number, it should be between 1 and " 
              << omp_get_num_procs() << " \n";
    return Rcpp::List();
  }
#endif
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  Individual cocktail, mutatedIndividual;
  
  //for the moment the distribution is bounded by 0 and RRmax
  //const int distribSize = 300;
  std::vector<unsigned int> score_distribution((max_score*10) +1); // +1 stand for every RR over the RRmax value
  unsigned int nbCocktailNotInPopulation = 0;
  std::vector<double> outstanding_score{};
  outstanding_score.reserve(10);
  
  std::vector<unsigned int> score_distribution_beta((max_score*10) +1);
  
  //used in the phyper function
  int ADRCount = 0;
  for(const auto& adr : observationsADR){
    if(adr)
      ++ADRCount;
  }
  int notADRCount = observationsMedication.size() - ADRCount;
    
  std::pair<double, std::pair<int, int>> currentGeom = std::make_pair(0.0,std::make_pair(0,0));
  std::pair<double, std::pair<int, int>> computeGeomOutput; // pair< phypergeometric, <N° of people taking the cocktail and having the ADR, N° of people taking the cocktail>>
  double minGeom = 0;
  double minGeomBeta = 0;

  //acceptance rate
  int acceptedMove = 0;
  int accepted_type1 =0;
  int accepted_type2 =0;
  int type1_move=0, type2_move=0;
  int type1_move_inF = 0, type2_move_inF = 0;
  int falseAcceptedCocktailCount = 0, falseSampledCocktailCount= 0;
  
  double RRx_k, RRy_k, pMutation, pAcceptation, pDraw;
  
  std::vector<std::pair<int,int>> vertexX;
  std::vector<std::pair<int,int>> vertexY;
  int chosenVertexidx;
  
  std::vector<std::pair<Individual,double>> bestResults;
  bestResults.reserve(nbResults);
  std::vector<std::pair<Individual,double>> bestResultsBeta;
  bestResultsBeta.reserve(nbResults);
  
  std::pair<Individual, double> currentResult;
  
  //if p.second is greater than 0, it means that the cocktail correspond to at least one person
  auto belongToF = [](const std::pair<double,std::pair<int,int>>& p){
    return (p.second.second > 0);
  };

  //initialization (every cocktail has to contain Smax medication)
  do{
    cocktail = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature)[0];
    computeGeomOutput = cocktail.computePHypergeom(observationsMedication, observationsADR,
                                                   upperBounds, ADRCount, notADRCount,
                                                   max_score, num_thread);
  } while (!belongToF(computeGeomOutput));
  currentGeom = computeGeomOutput;
  
#ifdef _OPENMP
  Rcpp::Rcout<< "openMP available \n";
#endif
  
  minGeom = currentGeom.first;
  bestResults.emplace_back(cocktail, currentGeom.first);
  if(currentGeom.second.second > beta){
    minGeomBeta = currentGeom.first;
    bestResultsBeta.emplace_back(cocktail, currentGeom.first);
  }
  
  for(int i = 0; i < epochs; ++i){
      pMutation = Rcpp::runif(1,0,1)[0];
      
      if(pMutation < p_type1){
        //type 1 mutation
        RRx_k = currentGeom.first;
        
        //here the current type 1 mutation consist in drawing a new cocktail of the same size
        mutatedIndividual = newIndividualWithCocktailSize(ATCtree.nrow(), Smax, 1, temperature)[0];
        
        computeGeomOutput = mutatedIndividual.computePHypergeom(observationsMedication, observationsADR,
                                                                upperBounds, ADRCount,
                                                                notADRCount,
                                                                max_score,
                                                                num_thread);
        
        RRy_k = computeGeomOutput.first;

        //to have an overview of the explored space (not in this method for the moment)
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        
        if(belongToF(computeGeomOutput)){
          // with this mutation, our ration q(X|Y) / q(Y|X) = 1
          pAcceptation = exp(((RRy_k - RRx_k)/static_cast<double>(cocktail.getTemperature()))); 
          pDraw = Rcpp::runif(1,0,1)[0];
          ++type1_move_inF;
          if(pAcceptation > pDraw){
            cocktail = mutatedIndividual;
            currentGeom = computeGeomOutput;
            ++acceptedMove;
            ++accepted_type1;
          }
        }else{
          ++nbCocktailNotInPopulation;
        }
        
        ++type1_move;
      }
      else{
        //type 2 mutation
        //if the selected indivudual is empty, the type 2 mutation cannot occur

        RRx_k = currentGeom.first;
        //get every vertex 0/1 for this patient
        vertexX = cocktail.getVertexList(ATCtree);
        
        chosenVertexidx = trunc(Rcpp::runif(1,0,vertexX.size())[0]);
        chosenVertexidx = chosenVertexidx == vertexX.size() ? vertexX.size()-1 : chosenVertexidx;
        
        std::pair<int,int> chosenVertex = vertexX[chosenVertexidx];
        
        mutatedIndividual = type2Mutation(cocktail, ATCtree.nrow(), chosenVertex);
        
        computeGeomOutput = mutatedIndividual.computePHypergeom(observationsMedication, observationsADR,
                                                                upperBounds, ADRCount,
                                                                notADRCount,
                                                                max_score, 
                                                                num_thread);
        RRy_k = computeGeomOutput.first;
        vertexY = mutatedIndividual.getVertexList(ATCtree);
        
        //to have an overview of the explored space
        //addPairToSet(mutatedIndividual_k, exploredPairs);
        
        if(belongToF(computeGeomOutput)){
          pAcceptation = (exp(((RRy_k - RRx_k) / cocktail.getTemperature()))) * 
            (static_cast<double>(vertexX.size())/static_cast<double>(vertexY.size()));
          
          pDraw = Rcpp::runif(1,0,1)[0];
          ++type2_move_inF;
          if(pAcceptation > pDraw){
            cocktail = mutatedIndividual;
            currentGeom = computeGeomOutput;
            ++acceptedMove;
            ++accepted_type2;
          }
        }else{
          ++nbCocktailNotInPopulation;
        }
        ++type2_move;
      
    }
    if(currentGeom.first < max_score){
      int index = 10 * currentGeom.first;
      ++score_distribution[index];
      // in every cases we add the RR to the "normal returned" distribution, and if
      // more than beta persons take it, we add the RR to the other ditribution named
      // score_distribution_beta
      if(currentGeom.second.second > beta){ // second.first = N° of people taking the cocktail and having the ADR
        ++score_distribution_beta[index];
      }
    }
    else{ // since we have an RR max, we just increment the last elements of the distribution
      // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
      outstanding_score.push_back(currentGeom.first); // could remove this line ?
      ++score_distribution[score_distribution.size()-1];
    }
    
    if(!cocktail.isTrueCocktail(upperBounds))
      falseAcceptedCocktailCount++;
    if(!mutatedIndividual.isTrueCocktail(upperBounds))
      falseSampledCocktailCount++;
        
    currentResult = std::make_pair(cocktail, currentGeom.first);
    //adding the result to the best result if we need to 
    minGeom = addToBestCocktails(bestResults, currentResult, nbResults, minGeom,
                               upperBounds);
    
    if(currentGeom.second.second > beta){
      minGeomBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                       minGeomBeta, upperBounds);
    }
  }
  
  //create the returned vector
  std::vector<std::vector<int>> returnedMed{};
  returnedMed.reserve(bestResults.size());
  std::vector<double>returned_score{};
  returned_score.reserve(bestResults.size());
  
  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returned_score.push_back(pair.second);
  }
  
  //create the returned vector with the cocktail taken by more than beta person
  std::vector<std::vector<int>> returnedMedBeta{};
  returnedMedBeta.reserve(bestResultsBeta.size());
  std::vector<double>returned_scoreBeta{};
  returned_scoreBeta.reserve(bestResultsBeta.size());
  
  for(const auto &pair : bestResultsBeta){
    returnedMedBeta.push_back(pair.first.getMedications());
    returned_scoreBeta.push_back(pair.second);
  }
  
  if(verbose){
    Rcpp::Rcout << "acceptance rate : " << static_cast<double>(acceptedMove) / static_cast<double>(epochs)<< "\n";
    Rcpp::Rcout << "acceptance rate type1 mutation : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move)<< "\n";
    Rcpp::Rcout << "acceptance rate type2 mutation : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move)<< "\n";
    Rcpp::Rcout << "acceptance rate type1 mutation when the proposal is in F : " << static_cast<double>(accepted_type1) / static_cast<double>(type1_move_inF)<< "\n";
    Rcpp::Rcout << "acceptance rate type2 mutation when the proposal is in F : " << static_cast<double>(accepted_type2) / static_cast<double>(type2_move_inF)<< "\n";
    Rcpp::Rcout << "number of proposed cocktail that was taken by nobody in the population : " << nbCocktailNotInPopulation << '\n';
    Rcpp::Rcout << "number of false cocktail sampled : " << falseSampledCocktailCount << '\n';
    Rcpp::Rcout << "number of false cocktail concidered in the distribution during the run : " << falseAcceptedCocktailCount << '\n';
  }
  
  return Rcpp::List::create(Rcpp::Named("ScoreDistribution") = score_distribution, Rcpp::Named("Outstanding_score") = outstanding_score, 
                            Rcpp::Named("Best_cocktails") = returnedMed, Rcpp::Named("Best_scores") = returned_score,
                            Rcpp::Named("Filtered_score_distribution") = score_distribution_beta,
                            Rcpp::Named("Best_cocktails_beta") = returnedMedBeta,
                            Rcpp::Named("Best_scores_beta") = returned_scoreBeta,
                            Rcpp::Named("cocktailSize") = Smax);
}

//'Genetic algorithm, trying to reach riskiest cocktails (the ones which maximize
//'the fitness function, Hypergeometric score in our case)
//'
//'@param epochs : number of step or the algorithm 
//'@param nbIndividuals : size of the population
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : real observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'@param diversity : enable the diversity mechanism of the algorithm
//' (favor the diversity of cocktail in the population),  default is false
//'@param p_crossover : probability to operate a crossover on the crossover phase. Default is 80\%
//'@param p_mutation : probability to operate a mutation after the crossover phase. Default is 1\%
//'@param nbElite : number of best individual we keep from generation to generation. Default is 0
//'@param tournamentSize : size of the tournament (select the best individual 
//'between tournamentSize sampled individuals) 
//'@param alpha : when making a type 1 mutation you have (alpha / size of cocktail) chance to add a drug. 
//'@param summary : print the summary of population at each steps ? 
//'
//'@return If no problem, return a List :
//' - meanFitnesses : The mean score of the population at each epochs of the algorithm.
//' - BestFitnesses : The best score of the population at each epochs of the algorithm.
//' - FinalPopulation : The final population of the algorithm when finished (medications
//' and corresponding scores)
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' results = GeneticAlgorithm(epochs = 10, nbIndividuals = 10, 
//'             ATCtree = ATC_Tree_UpperBound_2024,
//'             observations = FAERS_myopathy)
//'}
//'@export
//[[Rcpp::export]]
Rcpp::List GeneticAlgorithm(int epochs, int nbIndividuals, const DataFrame& ATCtree, 
                            const DataFrame& observations, int num_thread = 1, 
                            bool diversity = false, double p_crossover = .80,
                            double p_mutation = .01, int nbElite = 0, 
                            int tournamentSize = 2, double alpha = 1,
                            bool summary = true){ 
  //arguments verification
  if(p_crossover > 1 || p_crossover < 0 || nbIndividuals < 1 || p_mutation > 1 || p_mutation < 0 || epochs < 1){
    Rcpp::Rcerr << "problem in the values of the parameter in the call of this function \n";
    return Rcpp::List();
  }
  if(tournamentSize > nbIndividuals || nbElite > nbIndividuals){
    Rcpp::Rcerr << "the tournament size and the nbElite parameters must be less equal than the number of individuals \n";
    return Rcpp::List();
  }
  // OMP SET NUM THREAD = k, s.t. 1 <= k <= omp_get_num_procs()
#ifdef _OPENMP
  if(num_thread < 1 || num_thread > omp_get_num_procs()){
    Rcpp::Rcerr << "Wrong thread number, it should be between 1 and " 
              << omp_get_num_procs() << " \n";
    return Rcpp::List();
  }
#endif
  
  
  
  //since there is a diversity mechanism, we may not want to set an upper bound
  // to the metric (here Phypergeometric), so we put it equal to INT_MAX
  // we may want to put this as a parameter
  int max_metric = INT_MAX;
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
    
  int ADRCount = 0;
  for(const auto& adr : observationsADR){
    if(adr)
      ++ADRCount;
  }
  int notADRCount = observationsMedication.size() - ADRCount;
  
  //generate the initial population randomly (do we consider an Smax ?)
  double meanMedicationPerObs = meanMedications(observationsMedication) - 1;
  Population population(ATCtree.nrow(), nbIndividuals, meanMedicationPerObs,
                        upperBounds);
  Population matingPool(nbIndividuals);
  
  std::vector<double> meanScore;
  meanScore.reserve(epochs);
  std::vector<double> bestScore;
  bestScore.reserve(epochs);
  int bestIndividualIndex;
  
  int remainingDraw = 0;
  
  std::vector<double> score_before_penalization;
  score_before_penalization.reserve(population.getIndividuals().size());
  
  //here we may want to have a more sophisticated stopping condition (like, if the RR is 
  //significantly high given the previous calculated distribution)
  for(int i =0; i < epochs; ++i){
    
    //do we apply the diversity mechanism ?
    if(diversity){
      population.penalize(depth,father);
    }
    
    //1st : fit every individual
    population.evaluate(observationsMedication, observationsADR, upperBounds,
                        ADRCount, notADRCount, max_metric, num_thread);
    
    for(const auto& ind : population.getIndividuals()){
      score_before_penalization.push_back(ind.first);
    }
    
    
    //do we apply the diversity mechanism ?
    if(diversity){
      population.penalize(depth,father);
    }
    
    //2nd : make a selection and apply the crossover on the selected individual
    //keep the elite first
    population.keepElite(nbElite, matingPool);

    //select the individual according to the fitness
    remainingDraw = nbIndividuals - nbElite;
    population.tournamentSelection(tournamentSize, matingPool, remainingDraw);

    //operate a crossover over the mating pool 
    matingPool.crossover(nbElite, ATClength, upperBounds, ATCtree, p_crossover);
    
    //operate a mutation over the mating pool
    matingPool.mutate(nbElite, p_mutation, ATCtree, upperBounds, alpha);
    

    //3rd : replace the population
    population = matingPool;
    matingPool.clear();
    
    meanScore.push_back(population.getMean());
    bestIndividualIndex = population.bestIndividual();
    bestScore.push_back(population.getIndividuals()[bestIndividualIndex].first);
    
    if(summary)
      population.printSummary(i, meanScore[i], bestIndividualIndex);
    
    score_before_penalization.clear();
    
  }
  population.evaluate(observationsMedication, observationsADR, upperBounds,
                      ADRCount, notADRCount, max_metric, num_thread);

  //output the population 
  std::vector<std::vector<int>> medications;
  medications.reserve(nbIndividuals);
  std::vector<double> populationScore;
  populationScore.reserve(nbIndividuals);
  
  medications = population.getMedications();
  populationScore = population.getRR();

  
  Rcpp::List returnedPop = Rcpp::List::create(Rcpp::Named("cocktails") = medications,
                                              Rcpp::Named("score") = populationScore);
  
  return Rcpp::List::create(Rcpp::Named("meanFitnesses") = meanScore,
                            Rcpp::Named("BestFitnesses") = bestScore,
                            Rcpp::Named("FinalPopulation") = returnedPop);
}


//'The true distribution of the score among every single nodes of the ATC
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//'@param beta : minimum number of person taking the cocktails in order to consider it
//'in the beta score distribution 
//'@param max_score : maximum number the score can take. Score greater than this 
//'one would be added to the distribution as the value max_score. Default is 1000
//'@param nbResults : Number of returned solution (Cocktail with the
//' best oberved score during the run), 100 by default
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'
//'@return Return a List containing :
//' - ScoreDistribution : the distribution of the score as an array with each cells
//' representing the number of risks =  (index-1)/ 10
//' - Filtered_score_distribution : Distribution containing score for cocktails taken by at
//' least beta patients.
//' - Outstanding_score : An array of the score greater than max_score,
//' - Best_cocktails : the nbResults bests cocktails encountered during the run.
//' - Best_cocktails_beta : the nbResults bests cocktails taken by at least beta patients
//' encountered during the run.
//' - Best_scores : Score corresponding to the Best_cocktails.
//' - Best_scores_beta : Score corresponding to the Best_cocktails_beta.
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' size_1_score_distribution = trueDistributionDrugs(ATCtree = ATC_Tree_UpperBound_2024,
//'             observations = FAERS_myopathy[1:100,], beta = 4)
//'}
//'@export
//[[Rcpp::export]]
Rcpp::List trueDistributionDrugs(const DataFrame& ATCtree, const DataFrame& observations,
                                 int beta, int max_score = 1000, int nbResults = 100,
                                 int num_thread = 1){
 
#ifdef _OPENMP
 if(num_thread < 1 || num_thread > omp_get_num_procs()){
   Rcpp::Rcerr << "Wrong thread number, it should be between 1 and " 
             << omp_get_num_procs() << " \n";
   return Rcpp::List();
 }
#endif
 
 Rcpp::List observationsMedicationTmp = observations["patientATC"];
 std::vector<std::vector<int>> observationsMedication;
 observationsMedication.reserve(observationsMedicationTmp.size());
 Rcpp::LogicalVector observationsADR = observations["patientADR"];
 std::vector<int> upperBounds = ATCtree["upperBound"];
 int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
 
 for(int i =0; i < observationsMedicationTmp.size(); ++i){
   observationsMedication.push_back(observationsMedicationTmp[i]);
 }
 
 int notADRCount = observationsMedication.size() - ADRCount;
 
 
 //for the moment the distribution is bounded by 0 and 1000
 const int distribSize = max_score *10;
 std::vector<unsigned int> score_distribution(distribSize);
 std::vector<unsigned int> score_distributionGreaterBeta(distribSize);
 
 std::vector<double> outstanding_score{};
 outstanding_score.reserve(10);
 std::vector<double> outstanding_scoreBeta{};
 outstanding_scoreBeta.reserve(10);
 
 std::pair<double, std::pair<int,int>> computePhyperOutput;
 
 double minPhyper = 0, minPhyperBeta = 0;
 std::vector<std::pair<Individual,double>> bestResults;
 bestResults.reserve(nbResults);
 std::vector<std::pair<Individual,double>> bestResultsBeta;
 bestResultsBeta.reserve(nbResults);
 std::pair<Individual, double> currentResult;
 
 Individual indiv{};
 indiv.setTemperature(1);
 int med;
 for(int i = 0 ; i < ATCtree.nrow() - 1 ; ++i){
   med = i;

   indiv.setMedications({med});
   computePhyperOutput = indiv.computePHypergeom(observationsMedication, observationsADR,
                                                 upperBounds, ADRCount, notADRCount,
                                                 8000, num_thread);
   if(computePhyperOutput.second.second > 0){ // if the cocktail belongs to F
     if(computePhyperOutput.first < max_score){
       int index = 10 * computePhyperOutput.first;
       ++score_distribution[index];
       if(computePhyperOutput.second.second > beta){
         ++score_distributionGreaterBeta[index];
       }
     }
     else{
       if(computePhyperOutput.second.second > beta){
         outstanding_scoreBeta.push_back(computePhyperOutput.first);
       }
       outstanding_score.push_back(computePhyperOutput.first);
     }
     
     //keep the best results
     currentResult = std::make_pair(indiv, computePhyperOutput.first);
     minPhyper = addToBestCocktails(bestResults, currentResult, nbResults, minPhyper,
                                    upperBounds);
     
     if(computePhyperOutput.second.second > beta){
       minPhyperBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                          minPhyperBeta, upperBounds);
     }
     
   
   }
 }
 
 //create the vector of the best cocktail
 std::vector<std::vector<int>> returnedMed{};
 returnedMed.reserve(bestResults.size());
 std::vector<double>returned_score{};
 returned_score.reserve(bestResults.size());
 
 for(const auto &pair : bestResults){
   returnedMed.push_back(pair.first.getMedications());
   returned_score.push_back(pair.second);
 }
 
 //create the returned vector of the best cocktail taken by more than beta person
 std::vector<std::vector<int>> returnedMedBeta{};
 returnedMedBeta.reserve(bestResultsBeta.size());
 std::vector<double>returned_scoreBeta{};
 returned_scoreBeta.reserve(bestResultsBeta.size());
 
 for(const auto &pair : bestResultsBeta){
   returnedMedBeta.push_back(pair.first.getMedications());
   returned_scoreBeta.push_back(pair.second);
 }
 
 return Rcpp::List::create(Rcpp::Named("ScoreDistribution") = score_distribution,
                           Rcpp::Named("Filtered_score_distribution") = score_distributionGreaterBeta,
                           Rcpp::Named("Outstanding_score") = outstanding_score,
                           Rcpp::Named("Best_cocktails") = returnedMed,
                           Rcpp::Named("Best_cocktails_beta") = returnedMedBeta,
                           Rcpp::Named("Best_scores") = returned_score,
                           Rcpp::Named("Best_scores_beta") = returned_scoreBeta);
}

//'The true distribution of the score among every size-two cocktails
//'
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//'@param beta : minimum number of person taking the cocktails in order to consider it
//'in the beta score distribution 
//'@param max_score : maximum number the score can take. Score greater than this 
//'one would be added to the distribution as the value max_score. Default is 1000
//'@param nbResults : Number of returned solution (Cocktail with the
//' best oberved score during the run), 100 by default
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'
//'@return Return a List containing :
//' - ScoreDistribution : the distribution of the score as an array with each cells
//' representing the number of risks =  (index-1)/ 10
//' - Filtered_score_distribution : Distribution containing score for cocktails taken by at
//' least beta patients.
//' - Outstanding_score : An array of the score greater than max_score,
//' - Best_cocktails : the nbResults bests cocktails encountered during the run.
//' - Best_cocktails_beta : the nbResults bests cocktails taken by at least beta patients
//' encountered during the run.
//' - Best_scores : Score corresponding to the Best_cocktails.
//' - Best_scores_beta : Score corresponding to the Best_cocktails_beta.
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' size_2_score_distribution = trueDistributionSizeTwoCocktail(ATCtree = ATC_Tree_UpperBound_2024,
//'             observations = FAERS_myopathy[1:100,], beta = 4)
//'}
//'@export
//[[Rcpp::export]]
Rcpp::List trueDistributionSizeTwoCocktail(const DataFrame& ATCtree, const DataFrame& observations,
                                           int beta, int max_score = 100, int nbResults = 100,
                                           int num_thread = 1){
  
#ifdef _OPENMP
  if(num_thread < 1 || num_thread > omp_get_num_procs()){
    Rcpp::Rcerr << "Wrong thread number, it should be between 1 and " 
              << omp_get_num_procs() << " \n";
    return Rcpp::List();
  }
#endif
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
  
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  int notADRCount = observationsMedication.size() - ADRCount;

  
  //for the moment the distribution is bounded by 0 and 30
  const int distribSize = max_score *10;
  std::vector<unsigned int> score_distribution(distribSize);
  std::vector<unsigned int> score_distributionGreaterBeta(distribSize);
  
  std::vector<double> outstanding_score{};
  outstanding_score.reserve(10);
  std::vector<double> outstanding_scoreBeta{};
  outstanding_scoreBeta.reserve(10);
  
  std::pair<double, std::pair<int,int>> computePhyperOutput;
  
  double minPhyper = 0, minPhyperBeta = 0;
  std::vector<std::pair<Individual,double>> bestResults;
  bestResults.reserve(nbResults);
  std::vector<std::pair<Individual,double>> bestResultsBeta;
  bestResultsBeta.reserve(nbResults);
  std::pair<Individual, double> currentResult;

  Individual indiv{};
  indiv.setTemperature(1);
  std::vector<int> med;
  med.resize(2);
  for(int i = 0 ; i < ATCtree.nrow() - 1 ; ++i){
    med[0] = i;
    for(int j = i+1 ; j < ATCtree.nrow(); ++j){
      med[1] = j;
      indiv.setMedications(med);
      //computeRROutput = indiv.computeRR(observationsMedication, observationsADR, upperBounds, 8000, num_thread);
      computePhyperOutput = indiv.computePHypergeom(observationsMedication, observationsADR,
                                                    upperBounds, ADRCount, notADRCount,
                                                    8000, num_thread);
      if(computePhyperOutput.second.second > 0){ // if the cocktail belongs to F
        if(computePhyperOutput.first < max_score){
          int index = 10 * computePhyperOutput.first;
          ++score_distribution[index];
          if(computePhyperOutput.second.second > beta){
            ++score_distributionGreaterBeta[index];
          }
        }
        else{
          if(computePhyperOutput.second.second > beta){
            outstanding_scoreBeta.push_back(computePhyperOutput.first);
          }
          // if we are on a good RR we have a huge probability to stay on it -> maybe add the current RR only if it does not belong to the vector already ?
          outstanding_score.push_back(computePhyperOutput.first);
        }
        
        //keep the best results
        currentResult = std::make_pair(indiv, computePhyperOutput.first);
        minPhyper = addToBestCocktails(bestResults, currentResult, nbResults, minPhyper,
                                   upperBounds);
        
        if(computePhyperOutput.second.second > beta){
          minPhyperBeta = addToBestCocktails(bestResultsBeta, currentResult, nbResults,
                                             minPhyperBeta, upperBounds);
        }
        
      }
    }
  }
  
  //create the vector of the best cocktail
  std::vector<std::vector<int>> returnedMed{};
  returnedMed.reserve(bestResults.size());
  std::vector<double>returned_score{};
  returned_score.reserve(bestResults.size());
  
  for(const auto &pair : bestResults){
    returnedMed.push_back(pair.first.getMedications());
    returned_score.push_back(pair.second);
  }
  
  //create the returned vector of the best cocktail taken by more than beta person
  std::vector<std::vector<int>> returnedMedBeta{};
  returnedMedBeta.reserve(bestResultsBeta.size());
  std::vector<double>returned_scoreBeta{};
  returned_scoreBeta.reserve(bestResultsBeta.size());
  
  for(const auto &pair : bestResultsBeta){
    returnedMedBeta.push_back(pair.first.getMedications());
    returned_scoreBeta.push_back(pair.second);
  }
  
  return Rcpp::List::create(Rcpp::Named("ScoreDistribution") = score_distribution,
                            Rcpp::Named("Filtered_score_distribution") = score_distributionGreaterBeta,
                            Rcpp::Named("Outstanding_score") = outstanding_score,
                            Rcpp::Named("Best_cocktails") = returnedMed,
                            Rcpp::Named("Best_cocktails_beta") = returnedMedBeta,
                            Rcpp::Named("Best_scores") = returned_score,
                            Rcpp::Named("Best_scores_beta") = returned_scoreBeta);
}



std::vector<double> MetricCalc_2(const std::vector<int> &cocktail, 
                                  const std::vector<int>& ATClength,
                                  const std::vector<int>& upperBounds,
                                  const std::vector<std::vector<int>>& observationsMedication,
                                  const Rcpp::LogicalVector& observationsADR,
                                  int ADRCount, int num_thread = 1){
   
   std::vector<double> solution;
   solution.reserve(5);
   
   Individual ind{cocktail};
   Individual ind_D1{{cocktail[0]}};
   Individual ind_D2{{cocktail[1]}};
   
   
   int n000 = 0,n001 = 0;
   int n100 = 0, n101 = 0;
   int n011 = 0, n010 = 0;
   int n111 = 0, n110 = 0;
   bool have_D1 = false, have_D2 = false;
   for(int i = 0 ; i < observationsMedication.size(); ++i){
     have_D1 = ind_D1.matches(observationsMedication[i], upperBounds);
     have_D2 = ind_D2.matches(observationsMedication[i], upperBounds);
     if(!have_D1 && !have_D2){
       if(observationsADR[i]){
         ++n001;
       }
       else{
         ++n000;
       }
     } 
     else if(have_D1){
       if(have_D2){
         if(observationsADR[i])
           ++n111;
         else
           ++n110;
       }else{
         if(observationsADR[i])
           ++n101;
         else
           ++n100;
       }
     }
     else{
       if(observationsADR[i]){
         ++n011;
       } 
       else{
         ++n010;
       } 
     }
   } 
   
   
   double D1_rate = static_cast<double>(n101) / (static_cast<double>(n100 + n101));
   double D2_rate = static_cast<double>(n011) / (static_cast<double>(n010 + n011));
   double background_rate = static_cast<double>(n001) / (static_cast<double>(n000 + n001));
   double D1_D2_rate = static_cast<double>(n111) / (static_cast<double>(n110 + n111));
   int n00 = n001 + n000;
   int n10 = n101 + n100;
   int n01 = n011 + n010;
   int n11 = n111 + n110;
   
   // RR
   double RR_cocktail = 0.0;
   if(n111 > 0){
     RR_cocktail = (D1_D2_rate) / 
       ((n001 + n011 + n101) / (static_cast<double>(n00 + n10 + n01)));
   }
   
   
   // PRR
   double signal_PRR = 0;
   
   double RR1 = ((n101 + n111) / (static_cast<double>(n11 + n10))) /
     ((n011 + n001) / (static_cast<double>(n01 + n00)));
   
   double RR2 = ((n011 + n111) / (static_cast<double>(n11 + n01))) /
     ((n101 + n001) / (static_cast<double>(n10 + n00)));
   
   double SD_D1;
   if((n101 != 0 || n111 !=0) && (n001 !=0 || n011 !=0)){
     SD_D1 = std::sqrt((1.0/(static_cast<double>(n101+n111))) -
       (1.0 / (static_cast<double>(n10+n11))) + (1.0/(static_cast<double>(n001+n011))) -
       (1.0/(static_cast<double>(n01+n00))));
   }else{
     SD_D1 = 0;
   }
   
   double SD_D2;
   if((n101 != 0 || n001 !=0) && (n111 !=0 || n011 !=0)){
     SD_D2 = std::sqrt((1.0/(static_cast<double>(n011+n111))) -
       (1.0 / (static_cast<double>(n01+n11))) + (1.0/(static_cast<double>(n001+n101))) -
       (1.0/(static_cast<double>(n10+n00))));
   }else{
     SD_D2 = 0;
   }
   
   double SD_D1D2;
   if((n101 + n001 + n011 !=0) && (n111 !=0)){
     SD_D1D2 = std::sqrt((1.0/(static_cast<double>(n111))) -
       (1.0 / (static_cast<double>(n11))) + (1.0/(static_cast<double>(n001+n101+n011))) -
       (1.0/(static_cast<double>(n10+n00+n01))));
   }else{
     SD_D1D2 = 0;
   }
   
   double PRR_D1_025 = std::exp((std::log(RR1) - 1.96*SD_D1));
   double PRR_D2_025 = std::exp((std::log(RR2) - 1.96*SD_D2));
   double PRR_D1D2_025 = std::exp((std::log(RR_cocktail) - 1.96*SD_D1D2));
   if(n111 > 0)
     signal_PRR = PRR_D1D2_025 > std::max(PRR_D1_025, PRR_D2_025) ? 1 : 0;
   else
     signal_PRR = 0;
   
   
   //CSS
   double PRR_D1_975 = std::exp((std::log(RR1) + 1.96*SD_D1));
   double PRR_D2_975 = std::exp((std::log(RR2) + 1.96*SD_D2));
   double CSS = 0.0;
   if(n111 > 0)
     CSS = PRR_D1D2_025 / std::max(PRR_D1_975,PRR_D2_975);
   
   
   
   //omega
   double odds_ratio_r00 = background_rate / (1 - background_rate);
   double odds_ratio_r10 = D1_rate / (1 - D1_rate);
   double odds_ratio_r01 = D2_rate / (1 - D2_rate);
   
   double s11 = 1 - (1.0 / (std::max(odds_ratio_r00, odds_ratio_r10) + 
                     std::max(odds_ratio_r00, odds_ratio_r01) - 
                     odds_ratio_r00 + 1));

   double omega = std::log2((n111 + 0.5) / (s11*n11 + 0.5));
   Rcpp::NumericVector tmp{0.975};
   double lower_bound_omegaIC = INT_MIN;
   if(n111 > 0)
     lower_bound_omegaIC = omega - (Rcpp::qnorm(tmp)[0] / (std::log(2) * std::sqrt(n111)));
   
   
   //phyper
   double phyper = ind.computePHypergeom(observationsMedication, observationsADR,
                                         upperBounds, ADRCount, 
                                         observationsMedication.size() - ADRCount,
                                         10000, num_thread).first;
   
   solution.push_back(lower_bound_omegaIC);
   solution.push_back(n110);
   solution.push_back(signal_PRR);
   solution.push_back(RR_cocktail);
   solution.push_back(n111);
   solution.push_back(phyper);
   solution.push_back(CSS);
   
   return solution;
 }

//'Function used in the reference article to compare diverse Disproportionality Analysis metrics 
//'
//'@param CocktailList : A list of cocktails on which the Disproportionality analysis metrics should be computed
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//'
//'@return Multiple DA metrics computed on CocktailList cocktails
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' cocktails = list(c(561, 904),
//'                c(1902, 4585)) # only size 2 cocktails allowed for this function
//' 
//' scores_of_cocktails = computeMetrics_size2(CocktailList = cocktails,
//'                               ATCtree = ATC_Tree_UpperBound_2024, 
//'                               observations = FAERS_myopathy[1:100,])
//'}
//'@export
//[[Rcpp::export]]
Rcpp::DataFrame computeMetrics_size2(const std::vector<std::vector<int>>& CocktailList,
                                     const DataFrame& ATCtree,
                                     const DataFrame& observations,
                                     int num_thread = 1 ){
  
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> upperBounds = ATCtree["upperBound"];
  
  Rcpp::LogicalVector observationsADR = observations["patientADR"];
  int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
  
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  
  
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }

  std::unordered_map<std::string, std::vector<double>> metrics{
    {"n110", {}},
    {"n111", {}},
    {"RR" , {}}, 
    {"PRR" , {}},
    {"CSS" , {}},
    {"omega_025" , {}},
    {"phyper" , {}}};
  
  std::vector<double> metricsResults;
  metricsResults.reserve(metrics.size());
  
  for(auto& cocktail : CocktailList){
    metricsResults = MetricCalc_2(cocktail, ATClength, upperBounds, 
                                observationsMedication, observationsADR, ADRCount,
                                num_thread);
   
    int i = 0;
    for(auto& pair : metrics){
      pair.second.push_back(metricsResults[i++]);
    }
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("RR") = metrics["RR"],
                                 Rcpp::Named("phyper") = metrics["phyper"],
                                 Rcpp::Named("PRR") = metrics["PRR"],
                                 Rcpp::Named("CSS") = metrics["CSS"],
                                 Rcpp::Named("omega_025") = metrics["omega_025"],
                                 Rcpp::Named("n110") = metrics["n110"],
                                 Rcpp::Named("n111") = metrics["n111"]);            
}  


void print_list_in_file(const Rcpp::List& resultsGeneticAlgorithm,
                        const std::string& filename){
  std::ofstream ost(filename, std::ios::app);
  if(!ost.is_open()){
    Rcpp::Rcerr << "erreur ouverture fichier \n";
  }
  
  Rcpp::List final_population = resultsGeneticAlgorithm["FinalPopulation"];
  std::vector<std::vector<int>> cocktails = final_population["cocktails"];
  std::vector<double> score = final_population["score"];
  
  std::vector<int> mean_score_per_epoch = resultsGeneticAlgorithm["meanFitnesses"];
  std::vector<int> best_score_per_epoch = resultsGeneticAlgorithm["BestFitnesses"];
  
  //print cocktail + score 
  int i = 0;
  for(const auto& vec : cocktails){
    for(int med : vec){
      ost << med << ' ';
    }
    ost << score[i++] << '\n';
  }
  
  //for each epoch we print "mean_score best_score"
  for(int j = 0 ; j < mean_score_per_epoch.size() ; ++j){
    ost << mean_score_per_epoch[j] << ' ' << best_score_per_epoch[j] << '\n';
  }
  ost.close();
}


//' This function can be used in order to try different set of parameters for the genetic
//' algorithm in a convenient way. This will run each combination of mutation_rate,
//' nb_elite and alphas possible nb_test_desired times. For each sets of parameters,
//' results will be saved in a file named according to the set of parameter. One
//' can regroup the results of each run in a csv file by using the print_csv function
//' specifying the names of each file that needs to be treated and the number of 
//' performed runs on each parameter set
//' 
//' @param epochs : the number of epochs for the genetic algorithm
//' @param nb_individuals : the size of the population in the genetic algorithm
//' @param ATCtree : ATC tree with upper bound of the DFS (without the root)
//' @param observations : observation of the AE based on the medications of each patients
//' (a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//' @param nb_test_desired : number of genetic algorithm runs on each sets of parameters
//' @param mutation_rate : a vector with each mutation_rate to be tested
//' @param nb_elite : a vector with each nb_elite to be tested
//' @param alphas : a vector with each alphas to be tested
//' @param path : the path where the resulting files should be written
//' @param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//' @return No return value, this function should output results of the runs of the 
//' genetic algorithm in a specific format supported by function print_csv
//' and p_value_csv_file. The files are outputed in path which is current 
//' directory by default.
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' # different parameter to test for
//' mutation_rate = c(.1,.2,.3)
//' nb_elite = c(0,1,2)
//' alphas = c(0.5,1,2)
//' hyperparam_test_genetic_algorithm(epochs = 2, nb_individuals = 2,
//'                               ATCtree = ATC_Tree_UpperBound_2024, 
//'                               observations = FAERS_myopathy,
//'                               nb_test_desired = 5, mutation_rate = mutation_rate,
//'                               nb_elite = nb_elite, alphas = alphas)
//'}
//[[Rcpp::export]]
void hyperparam_test_genetic_algorithm(int epochs, int nb_individuals, 
                                       const DataFrame& ATCtree, 
                                       const DataFrame& observations,
                                       int nb_test_desired, 
                                       const std::vector<double>& mutation_rate,
                                       const std::vector<int>& nb_elite,
                                       const std::vector<double>& alphas,
                                       const std::string& path = "./",
                                       int num_thread=1){
  for(const auto& mr : mutation_rate){
    for(const auto& ne : nb_elite){
      for(const auto& alpha : alphas){
        std::ostringstream filename;
        filename << path << epochs << "e_" << nb_individuals << "ind_" << mr << "mr_" <<
          ne << "ne_" << alpha << "alpha.txt";
        for(int i = 0; i < nb_test_desired; ++i){
          auto out = GeneticAlgorithm(epochs, nb_individuals, ATCtree, observations,
                                      num_thread, true, 0.8, mr, ne, 2, alpha,
                                      false);
          print_list_in_file(out, filename.str());
        }
      }
    }
  }
}

using answer_set = std::vector<std::vector<int>>;

std::vector<int> recup_cocktail(const std::string& line){
  std::istringstream stream(line.data());
  int medoc;
  std::vector<int> returned_vec;
  
  while(stream >> medoc){
    returned_vec.push_back(medoc);
  }
  returned_vec.pop_back();
  
  return returned_vec;
}

std::pair<double, std::vector<int>> recup_solution(const std::string& line){
  std::istringstream stream(line.data());
  
  std::vector<int> vec;
  double number;
  
  while(stream >> number){
    vec.push_back(number);
  }
  
  vec.pop_back();
  std::sort(vec.begin(),vec.end());
  return {number, vec};
}


//'Print every cocktails found during the genetic algorithm when used with the 
//'hyperparam_test_genetic_algorithm function. This enables to condense the solutions 
//'found in each files by collapsing similar cocktail in a single row by cocktail.
//'
//'
//' @param input_filenames : A List containing filename of hyperparam_test_genetic_algorithm output file
//' @param observations : observation of the AE based on the medications of each patients
//' (a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//' @param repetition : The parameter nb_test_desired used in the hyperparam test function
//' @param ATCtree : ATC tree with upper bound of the DFS (without the root)
//' @param csv_filename : Name of the output file, "solutions.csv" by default
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  data("FAERS_myopathy")
//'  files = c('250e_700ind_0.2mr_0ne_2alpha.txt') # results of hyperparam_test_genetic_algorithm
//' 
//'  print_csv(input_filenames = files, observations = FAERS_myopathy,
//'           repetition = 5, ATCtree = ATC_Tree_UpperBound_2024)
//' }
//' @return No return value, should process the output of the genetic algorithm in 
//' files produced by hyperparam_test_genetic_algorithm and output a summary csv file.
//' The csv file is outputed in current directory and named after the csv_filename
//' variable (solutions.csv by default).
//' @export
//[[Rcpp::export]]
void print_csv(const std::vector<std::string>& input_filenames,
               const DataFrame& observations,
               int repetition, const DataFrame& ATCtree,
               const std::string& csv_filename = "solutions.csv" ){
  Rcpp::List observationsMedicationTmp = observations["patientATC"];
  std::vector<std::vector<int>> observationsMedication;
  observationsMedication.reserve(observationsMedicationTmp.size());
  for(int i =0; i < observationsMedicationTmp.size(); ++i){
    observationsMedication.push_back(observationsMedicationTmp[i]);
  }
  
  Rcpp::LogicalVector observationsADR = observations["patientADR"];

  
  std::vector<std::pair<double, std::vector<int>>> solutions;
  
  for(const auto& filename : input_filenames){
    std::ifstream input(filename);
    if(!input.is_open()){
      Rcpp::Rcerr << "erreur ouverture du fichier " << filename << "\n";
      return;
    }
    
    int epochs = std::stoi(filename.substr(filename.find('/')+1,
                                           filename.find('e')).data());
    int nb_individuals = std::stoi(filename.substr(filename.find('_')+1,
                                                   filename.find('i')).data());
    
    std::string line;
    solutions.reserve(solutions.capacity() + repetition*nb_individuals);
    
    for(int i = 0; i < repetition; ++i){
      for(int j = 0; j < nb_individuals;++j){
        std::getline(input, line);
        auto tmp = recup_solution(line);
        if(tmp.first > 0)
          solutions.push_back(tmp);
      } 
      
      for(int j = 0; j < epochs; ++j){
        std::getline(input, line);
      } 
    }
    
    input.close();
  }
  
  std::set<
    std::pair<double, std::vector<int>>,
    std::greater<std::pair<double, std::vector<int>>>
  > set_sol(solutions.begin(), solutions.end());
  
  std::ofstream output(csv_filename);
  if(!output.is_open()){
    Rcpp::Rcerr << "erreur ouverture du fichier " << csv_filename << "\n";
    return;
  }
  std::vector<std::string> ATCName = ATCtree["Name"];
  
  output << "score ; Cocktail ; n patient taking C ; n patient taking C and having AE ; RR \n";
  
  for(const auto& sol : set_sol){
    Individual c{sol.second};
    auto pair = c.computePHypergeom(observationsMedication, observationsADR,
                                    ATCtree["upperBound"], 1,1,1,1).second;
    
    auto RR = c.computeRR(observationsMedication, observationsADR, 
                            ATCtree["upperBound"], 100000);
    output << sol.first << ";";
    for(auto ite = sol.second.begin(); ite != sol.second.end()-1; ++ite){
      output << ATCName[*ite] << ":"; 
    }
    output << ATCName[*(sol.second.end()-1)] << ";";
    output << pair.second << ";" << pair.first << ";" << RR.first << "\n";
  }
  
  output.close();
}

std::vector<std::vector<double>> dissim(const Population& pop,
                                        const std::vector<int>& depth,
                                        const std::vector<int>& father,
                                        bool normalization){
  IntMatrix M;
  std::vector<int> indexM;
  
  std::tie(M, indexM) = pop.pretraitement(depth,father);
  RealMatrix D = pop.dissimilarity(M, indexM, normalization);
  
  return D;
}

//' Recover the square matrix of distance between cocktails where the index (i,j)
//' of the matrix is the distance between cocktails i and j in the genetic_results
//' list. 
//' @param genetic_results the List returned by the genetic algorithm.
//' @param ATCtree : ATC tree with upper bound of the DFS (without the root)
//' @param normalization : Do we keep the distance between cocktail in the range [0;1] ? 
//' 
//' @return The square matrix of distances between cocktails
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  data("FAERS_myopathy")
//'  
//'  genetic_results = GeneticAlgorithm(epochs = 10, nbIndividuals = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024,
//'             observations = FAERS_myopathy)
//'  distance_matrix = get_dissimilarity_from_genetic_results(genetic_results = genetic_results,
//'                         ATCtree = ATC_Tree_UpperBound_2024, normalization = TRUE)
//' }
//' @export
//[[Rcpp::export]]
std::vector<std::vector<double>> get_dissimilarity_from_genetic_results(const Rcpp::List& genetic_results,
                                                   const DataFrame& ATCtree,
                                                   bool normalization){
  Rcpp::List population_list = genetic_results["FinalPopulation"];
  std::vector<std::vector<int>> cocktails = population_list["cocktails"];
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  Population population(cocktails);

  return dissim(population, depth, father, normalization);
}

//' Recover the square matrix of distance between cocktails where the index (i,j)
//' of the matrix is the distance between cocktails i and j in the csv file containing
//' results of genetic algorithm
//' 
//' @param filename : the name of the file returned by the print_csv function.
//' @param ATCtree : ATC tree with upper bound of the DFS (without the root)
//' @param normalization : Do we keep the distance between cocktail in the range [0;1] ? 
//' 
//' @return The square matrix of distances between cocktails
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  
//'  distance_matrix = get_dissimilarity_from_txt_file(filename = '250e_700ind_0.2mr_0ne_2alpha.txt',
//'                         ATCtree = ATC_Tree_UpperBound_2024, normalization = TRUE)
//' }
//' @export
//[[Rcpp::export]]
 std::vector<std::vector<double>> get_dissimilarity_from_txt_file(
                                                   const std::string& filename,
                                                   const DataFrame& ATCtree,
                                                   bool normalization = true){
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<std::vector<int>> cocktails;
  std::vector<int> current_cocktail;
  
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  std::string line;
  std::ifstream input(filename);
  if(!input.is_open()){
    Rcpp::Rcerr << "erreur ouverture du fichier " << filename << "\n";
    return {};
  }
  
  std::getline(input,line); // we don't want the first line as it is a title
  
  while(std::getline(input, line)){
    int beg_cocktail = line.find('|');
    int end_cocktail = line.substr(beg_cocktail+1).find('|');
    std::string cocktail = line.substr(beg_cocktail + 1, end_cocktail);
    // 3.8 is the expected number of char in the id of a randomly picked 
    // drug in the ATC tree (2014tree)
    int approximated_cocktail_size = cocktail.size() / 3.8;
    current_cocktail.reserve(approximated_cocktail_size);
    
    int med;
    std::istringstream iss(cocktail);
    while(iss >> med){
      current_cocktail.push_back(med);
    }
    cocktails.push_back(current_cocktail);
    current_cocktail.clear();
    
  }
  
  Population population(cocktails);
  
  return dissim(population, depth, father, normalization);
}


//' Recover the square matrix of distance between cocktails where the index (i,j)
//' of the matrix is the distance between cocktails i and j in an arbitrary
//' cocktail list
//' 
//' @param cocktails : A list of cocktails in the form of a vector of integer
//' @param ATCtree : ATC tree with upper bound of the DFS (without the root)
//' @param normalization : Do we keep the distance between cocktail in the range [0;1] ? 
//' 
//' @return The square matrix of distances between cocktails
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' 
//' cocktails = list(c(561, 904),
//'                c(1902, 4585)) # only size 2 cocktails allowed for this function
//' 
//' distance_matrix = get_dissimilarity_from_cocktail_list(cocktails = cocktails,
//'                               ATCtree = ATC_Tree_UpperBound_2024, 
//'                               normalization = TRUE)
//'}
//[[Rcpp::export]]
std::vector<std::vector<double>> get_dissimilarity_from_cocktail_list(const std::vector<std::vector<int>>& cocktails,
                                                                const Rcpp::DataFrame& ATCtree,
                                                                bool normalization = true){
  std::vector<int> ATClength = ATCtree["ATC_length"];
  std::vector<int> depth, father;
  std::tie(depth, father) = treeDepthFather(ATClength);
  
  Population population(cocktails);
  
  return dissim(population, depth, father, normalization);
}

