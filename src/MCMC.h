#ifndef MCMC_H
#define MCMC_H


#include "Individual.h"
#include <memory>
#include <cmath>
#include <algorithm>
#include <set>
#include <numeric>

//' Get the average number of medications taken by a patient (made for the observations set)
//'
//'@param observations : a C++ vector of medications for each patient 
double meanMedications(const std::vector<std::vector<int>>& observations);

//' return the larger score and the corresponding individual
//' 
//' @param a list of 4 score to compare
std::pair<Individual,double> largerScore(const std::pair<Individual,double>& firstScore,const std::pair<Individual,double>& secScore,
                                      const std::pair<Individual,double>& thirdScore,const std::pair<Individual,double>& fourthScore);
 
//' add the score in the right box of the score distributions array
//' 
//' @param score : a score
//' @param the list of score distribution
void addScoretoDistribution(const double,std::vector<unsigned int>&);

//' add a pair cocktail/score in the bestResult set if it is needed
//' 
//' @param bestResults : The best current results of the random walk
//' @param currentResult : The current state (cocktail) and the associated score
//' @param nbResults : the number of results the user want to have at the end of the random walk
//' @param minScore : the minimum score contained in the bestResults vector
double addToBestCocktails(std::vector<std::pair<Individual,double>>& bestResults,
                        const std::pair<Individual,double>& currentResult,
                        int nbResults, double minScore, const std::vector<int>& upperBound);

//' Add the pair to the set of explored pair if needed
//' 
//' @param i : the individual that may be added
//' @param the current list of explored pairs
void addPairToSet(const Individual& i, std::set<std::pair<int,int>>&);
 
 //' Check if a Result is already in the result vector or not
 //' 
 //' @param bestResults : the current list of returned result (on the emc algorithm)
 //' @param bestResult : the result which needs to be tested
 //' 
 //' @return true if the result is already in the results vector, false otherwise
bool isNotInResultList(const std::vector<std::pair<Individual,double>>& bestResults,
                       const std::pair<Individual,double>& bestResult);
 
template<class T1>
bool PermutEqual(const std::pair<T1,T1>& P1, const std::pair<T1,T1>& P2){
 return ((P1.first == P2.first) && (P1.second == P2.second)) || ((P1.first == P2.second) && (P1.second == P2.first));
}

 //' Check if the Y cocktail could come from a type 2 mutation of the X cocktail
 //' 
 //' @param X : original cocktail
 //' @param Y : mmutated cocktail
 //' @param vertex : all the (1,0) vertex that are available for the cocktail X.
 //' 
 //' @return true if the Y cocktail could come from X using the type 2 mutation
bool mutatedByType2(std::vector<int> X, std::vector<int> Y, const std::vector<std::pair<int,int>>& vertex);
 
 
//'Create the individuals vector with starting individuals only
//'
//'@param startingInd : starting individuals given by the user
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_WOtemp(const Rcpp::List& startingInd);

//'Create the individuals vector with starting temperatures only
//'
//'@param treeSize : size of the ATC tree (to get the DFS index interval)
//'@param nbIndividuals : number of individuals in the population
//'@param meanMedic : the average number of medications took by the observations patients
//'@param temperatures : starting temperatures given by the user
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_WOIndividual(int treeSize, int nbIndividuals,double meanMedic, const Rcpp::NumericVector& temperatures);;

//' Overload of the DFtoCPP_WOtemp function here we don't have to specify the temperatures
//' every individual has the same temperature which is 1.
//' 
//'@param treeSize : size of the ATC tree (to get the DFS index interval)
//'@param nbIndividuals : number of individuals in the population
//'@param meanMedic : the average number of medications took by the observations patients
//'@return Individual vector which would be use by the emc algorithm and the GeneticAlgorithm
std::vector<Individual> DFtoCPP_WOIndividual(int treeSize, int nbIndividuals,double meanMedic);

//' Made for the distributionApproximation function, return an individual with a specified cocktail size
//' 
//'@param treeSize : size of the ATC tree (to get the DFS index interval)
//'@param cocktailSize : the number of drugs the individual must take
//'@param nbIndividuals : number of individuals in the population
//'@param temperature : the temperature of the individual
//'
//'@return Individual vector randomly initialized with cocktailSize drugs
std::vector<Individual> newIndividualWithCocktailSize(int treeSize, int cocktailSize, int nbIndividuals = 1, double temperature = 1);

//'Create the individuals vector with starting individuals and starting temperatures
//'
//'@param startingInd : starting individuals given by the user
//'@param startingTemp : starting temperatures given by the user
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_Wtemp(const Rcpp::List& startingInd,const Rcpp::NumericVector& startingTemp);

//'Create the individuals vector randomly
//'
//'@param treeSize : size of the ATC tree (to get the DFS index interval)
//'@param nbIndividuals : number of individuals in the population
//'@param meanMedic : the average number of medications took by the observations patients
//'
//'@return Individual vector which would be use by the emc algorithm
std::vector<Individual> DFtoCPP_WOtempAndIndividual(int treeSize, int nbIndividuals,double meanMedic);

//'Return the Pairs cocktails causing the ADR (used to check which part of the space is explored)
//'
//'@param observationsMed : The List containing the medications taken by the real person
//'@param ADR : A logical vector containing the ADR associated for each patients of the first array (observationsMed)
//'
//'@return the Pairs cocktails causing the ADR 
std::set<std::pair<int,int>> getADRPairs(const Rcpp::List& observationsMed, const Rcpp::LogicalVector& ADR);

//' Return a Mutated version of the individual in parameter (using the 1st mutation)
//' 
//' @param indiv : the individual chosen to be mutated
//' @param treeSize : size of the ATC tree
//' @param alpha : a hyperparameter allowing us to manage to probability of adding a drug to the cocktail. The probability
//' to add a drug to the cocktail is the following : $$ \alpha / n$$ Where n is the original size of the cocktail. 
//' 
//' @return the mutated individual
Individual type1Mutation(const Individual& indiv, int treeSize, double alpha, bool emptyCocktail);

//' Test adjusted type 1 mutation
Individual adjustedType1Mutation(const Individual& indiv, int treeSize, double alpha, bool emptyCocktail);

//' Return a mutated version of the individual in parameter (using the 2nd mutation)
//'
//'@param indiv : the individual chosen to be mutated
//'@param treeSize : size of the ATC tree
//'@param p : a pair of <int,int> representing the vertex being modified
//'
//'@return the mutated individual
Individual type2Mutation(const Individual& indiv, int treeSize, const std::pair<int,int>& p);

//' Return a Mutated version of the individual in parameter (using the crossover mutation)
//' 
//' @param indiv1 : the individual chosen to be mutated
//' @param indiv2 : the individual with which the indiv1 will be mutated
//' @param ATCtree : tree with every medication
//' @param selectedNode : represent the internal node of the tree on which we will perform the crossover
//' @param upperBound : the upper bound of the set to consider when performing a crossover. Note the interval 
//' to swap between indiv1 and indiv2 is [selectedNode ; upperBound[
//' 
//' @return the mutated individual
Individual crossoverMutation(const Individual& indiv,const Individual& indiv2,const Rcpp::DataFrame& ATCtree,
                             int selectedNode, int upperBound);

//'Return the depth and the father of each node of the ATC tree given in parameter
//'
//'@param ATC_length : the atc_length vector of the ATC tree
//'
//'@return : first = the depth of each node ; second = father of each node
std::pair<std::vector<int>,std::vector<int>> treeDepthFather(const std::vector<int>& ATC_length);

#endif
