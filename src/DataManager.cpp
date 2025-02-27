#include "RcppArmadillo.h"
#include "MCMC.h"
#include <vector>
#include <numeric>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <string>
using Rcpp::DataFrame;


// [[Rcpp::depends(RcppArmadillo)]]
//'Convert ATC Code for each patients to the corresponding DFS number of the ATC tree 
//'
//' @param tree : ATC tree (we assume that there is a column 'ATCCode' )
//' @param patientATC : patients observations, for each patient we got a string 
//'containing taken medications (ATC code)
//' @examples
//'  ATC_code <- c('A01AA30 A01AB03', 'A10AC30')
//'  ATCtoNumeric(ATC_code, ATC_Tree_UpperBound_2024)
//'
//' @return a matrix of the same size as patientATC but containing integer 
//' that are the index of the corresponding ATC code.
//' @export
// [[Rcpp::export]]
std::vector<std::vector<int>> ATCtoNumeric(const std::vector<std::string>&
  patientATC,const DataFrame& tree) {
  std::vector<std::string> cppTree= tree["ATCCode"];
  std::vector<std::vector<int>> newPatientATC;
  
  newPatientATC.reserve(patientATC.size());
  
  std::string codeI;
  std::string delimiter = " ";
  int posOcc = 0,posTree = 0;
  std::vector<int> patientI(0);
  for(const std::string& code : patientATC){
    patientI.clear();
    patientI.reserve(3);
    
    codeI = code;
    posOcc = 0;
    posTree = 0;
    
    while(posOcc > -1){
      posOcc = codeI.find(delimiter);
      
      //token is the ATC code being converted to numeric 
      std::string token = codeI.substr(0,posOcc);
      //new code without the token (because it is already converted)
      codeI = codeI.substr(posOcc+1);
      //first condition should be useless but could detect error -> when the ATC code is not found in the ATC tree
      while(posTree < cppTree.size() && token != cppTree[posTree]){
        ++posTree;
      }
      
      if(posTree == cppTree.size()){
        Rcpp::Rcerr<<"error : a patient take a medication that is not in the tree" << '\n';
        return {};
      }
      //+1 because of the cpp indexes (starting at 0)
      patientI.push_back(posTree);
    }
    patientI.shrink_to_fit();
    newPatientATC.push_back(patientI);
    
  }
  return newPatientATC;
  
}


//'Convert the histogram returned by the DistributionApproximation function, to a real number distribution
//'(that can be used in a test for example) 
//'
//'@param vec : distribution returned by the DistributionAproximationFunction
//'
//'@return A vector containing sampled risk during the MCMC algorithm 
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//'  DistributionApproximationResults = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy)
//'   histogramToDitribution(DistributionApproximationResults$ScoreDistribution)
//' }
//'@export
// [[Rcpp::export]]
Rcpp::NumericVector histogramToDitribution(const std::vector<int>& vec){
  std::vector<double> returnedVec;
  returnedVec.reserve(std::accumulate(vec.begin(),vec.end(),0));
  int count;
  for(int i = 0 ; i < vec.size(); ++i){
    count = vec[i];
    for(int j = 0 ; j < count ; ++j){
      returnedVec.push_back(static_cast<double>(i)/10.0);
    }
  }
  return Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(returnedVec));
}

//' Output the outstanding score (Outstanding_score) outputed by the MCMC algorithm
//' in a special format
//' 
//' @param outstanding_score : Outstanding_score outputed by MCMC algorithm to be converted
//' to the ScoreDistribution format
//' @param max_score : max_score parameter used during the MCMC algorithm
//' 
//' @return outstanding_score in a format compatible with MCMC algorithm output
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  data("FAERS_myopathy")
//' 
//'   DistributionApproximationResults = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy)
//'   OutsandingScoreToDistribution(DistributionApproximationResults$Outstanding_score, max_score = 100)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector OutsandingScoreToDistribution(const std::vector<double>& outstanding_score, int max_score){
  std::vector<double> returnedVec;
  returnedVec.resize((max_score*10)+1);
  
  for(const auto& score: outstanding_score){
    int index;
    if(score < max_score){
      index= score*10;
    }
    else{
      index = returnedVec.size()-1;
    }
    ++returnedVec[index];
  }
  
  return Rcpp::wrap(returnedVec);
}


bool hasExtension(const std::string& filename, const std::string& extension){
  size_t dotPos = filename.rfind('.');
  
  std::string fileExtension = dotPos == std::string::npos ? "" : filename.substr(dotPos);
  
  return fileExtension == extension;
}


std::vector<std::vector<std::string>> read_csv_genetic(std::ifstream& ifstr, char sep = ';'){
  std::vector<std::vector<std::string>> file;
  
  file.reserve(500);

  std::string line;
  //read the header of the csv file
  //std::getline(ifstr,line);
  
  while(std::getline(ifstr,line)){
    std::vector<std::string> row;
    row.reserve(6);
    
    std::stringstream ss(line);
    std::string cell;
    
    while(std::getline(ss, cell, sep)){
      row.push_back(cell);
    }
    file.push_back(row);
  }
  
  file.shrink_to_fit();
  return file;
}


std::vector<std::vector<std::string>> check_extension_and_read_csv(
  const std::string& filename, char sep){
  
  std::vector<std::vector<std::string>> file;
  if(!hasExtension(filename,".csv")){
    Rcpp::Rcerr << "file extension not supported for now \n";
    file.resize(0);
    return file;
  }
  
  std::ifstream ifstr(filename);
  if(!ifstr.is_open()){
    Rcpp::Rcerr << "the file " << filename << " has failed to open\n";
    file.resize(0);
    return file;
  }
  
  file = read_csv_genetic(ifstr, sep);
  
  ifstr.close();
  
  return file;
}

//'Function used to compute the Relative Risk on a list of cocktails
//'
//'@param cocktails : A list containing cocktails in the form of vector of integers (ATC index)
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//' 
//'@return RR score among "cocktails" parameters
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' cocktails = list(c(561, 904),
//'                c(1902, 4585))
//' 
//' RR_of_cocktails = compute_RR_on_list(cocktails = cocktails,
//'                               ATCtree = ATC_Tree_UpperBound_2024, 
//'                               observations = FAERS_myopathy)
//'}
//'@export
//[[Rcpp::export]]
std::vector<double> compute_RR_on_list(const std::vector<std::vector<int>> &cocktails, 
                                      const DataFrame& ATCtree, 
                                      const DataFrame& observations,
                                      int num_thread = 1)
{
 std::vector<double> RR;
 RR.reserve(cocktails.size());
 Rcpp::LogicalVector observationsADR = observations["patientADR"];
 std::vector<std::vector<int>> observationMed = observations["patientATC"];
 
 std::vector<int> upperBounds = ATCtree["upperBound"];
 
 for(const auto& cocktail : cocktails){
   RR.push_back(Individual(cocktail).computeRR(observationMed,
                           observationsADR,
                           upperBounds, 100000,
                           num_thread).first);
 }
 return RR;
}


//'Function used to compute the Hypergeometric score on a list of cocktails
//'
//'@param cocktails : A list containing cocktails in the form of vector of integers (ATC index)
//'@param ATCtree : ATC tree with upper bound of the DFS (without the root)
//'@param observations : observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//'@param num_thread : Number of thread to run in parallel if openMP is available, 1 by default
//' 
//'@return Hypergeometric score among "cocktails" parameters
//'@examples
//'\donttest{
//' data("ATC_Tree_UpperBound_2024")
//' data("FAERS_myopathy")
//' 
//' cocktails = list(c(561, 904),
//'                c(1902, 4585))
//' 
//' Hypergeom_of_cocktails = compute_hypergeom_on_list(cocktails = cocktails,
//'                               ATCtree = ATC_Tree_UpperBound_2024, 
//'                               observations = FAERS_myopathy)
//'}
//'@export
//[[Rcpp::export]]
std::vector<double> compute_hypergeom_on_list(const std::vector<std::vector<int>> &cocktails, 
                                             const DataFrame& ATCtree, 
                                             const DataFrame& observations,
                                             int num_thread = 1)
{
 std::vector<double> Phyper;
 Phyper.reserve(cocktails.size());
 Rcpp::LogicalVector observationsADR = observations["patientADR"];
 std::vector<std::vector<int>> observationsMedication = observations["patientATC"];
 
 std::vector<int> upperBounds = ATCtree["upperBound"];
 int ADRCount = std::count(observationsADR.begin(), observationsADR.end(), true);
 int patient_number = observationsMedication.size();
 
 for(const auto& cocktail : cocktails){
   Phyper.push_back(Individual(cocktail).
                      computePHypergeom(observationsMedication, observationsADR,
                                        upperBounds, ADRCount, 
                                        patient_number - ADRCount,
                                        10000, num_thread).first);
 }
 
 return Phyper;
}


//' Used to add the p_value to each cocktail of a csv_file that is an
//' output of the genetic algorithm
//' @param distribution_outputs A list of distribution of cocktails of different sizes
//' in order to compute the p_value for multiple cocktail sizes
//' @param filename The file name of the .csv file containing the output
//' @param filtred_distribution Does the p-values have to be computed using filtered distribution
//' or normal distribution (filtered distribution by default)
//' @param sep The separator used in the csv file (';' by default)
//' 
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  data("FAERS_myopathy")
//' 
//'   DistributionApproximationResults_size2 = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy, Smax = 2)
//'             
//'   DistributionApproximationResults_size3 = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy, Smax = 3)
//'             
//'   score_distribution_list = list(DistributionApproximationResults_size2,
//'                               DistributionApproximationResults_size3)
//'   p_value_csv_file(score_distribution_list, "path/to/output.csv")
//' }
//' @return A real valued number vector representing the p-value of the inputed
//' csv file filename, computed on the distribution_outputs List.
//' @export
//[[Rcpp::export]]
void p_value_csv_file(const std::vector<Rcpp::List>& distribution_outputs, const std::string& filename,
                                            bool filtred_distribution = false,
                                            const std::string & sep = ";"){
  
  std::vector<std::vector<std::string>> file = check_extension_and_read_csv(filename, sep[0]);
  if(file.size() == 0 ){
    Rcpp::Rcerr << "No cocktail to recover\n";
    return;
  }
  
  Rcpp::Function compute_p_value = Rf_findFun(Rf_install("p_value_on_sampled"),
                                              R_GlobalEnv);
  
  file[0].push_back(" p_value");
  
  for(const Rcpp::List& list : distribution_outputs){
    int k = list["cocktailSize"];
    
    
    
    for(auto it = file.begin()+1 ; it != file.end(); ++it){
      if(std::count((*it)[1].begin(), (*it)[1].end(), ':') == k-1){
        (*it).push_back(std::to_string(
                          Rcpp::as<double>(
                            compute_p_value(list, std::stod((*it)[0]),filtred_distribution)
                          )
                        ));
      }
    }
    
  }

  std::ofstream ofstr(filename);
  if(!ofstr.is_open()){
    Rcpp::Rcerr << "the file " << filename << " has failed to open (to output results)\n";
    return;
  }

  
  for(const auto& line : file){
    int line_size = line.size();
    for(int i = 0; i < line_size -1; ++i){
      ofstr << line[i] << ';' ;
    }
    ofstr << line[line_size-1] << "\n";
  }
  
}

//' Used to add the p_value to each cocktail of an output of the genetic algorithm
//' @param distribution_outputs A list of distribution of cocktails of different sizes
//' in order to compute the p_value for multiple cocktail sizes
//' @param genetic_results outputs of the genetic algorithm
//' @param filtred_distribution Does the p-values have to be computed using filtered distribution
//' or normal distribution (filtered distribution by default)
//' 
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  data("FAERS_myopathy")
//'   DistributionApproximationResults_size2 = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy, Smax = 2)
//'             
//'   DistributionApproximationResults_size3 = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy, Smax = 3)
//'             
//'   score_distribution_list = list(DistributionApproximationResults_size2,
//'                               DistributionApproximationResults_size3)
//'   genetic_results = GeneticAlgorithm(epochs = 10, nbIndividuals = 20, 
//'             ATCtree = ATC_Tree_UpperBound_2024,
//'             observations = FAERS_myopathy)
//'   p_value_genetic_results(score_distribution_list, genetic_results)
//' }
//' @return A real valued number vector representing the p-value of the inputed
//' genetic algorithm results (genetic_results) computed on the 
//' distribution_outputs List.
//' @export
//[[Rcpp::export]]
std::vector<double> p_value_genetic_results(const std::vector<Rcpp::List>& distribution_outputs, 
                      const Rcpp::List& genetic_results,
                      bool filtred_distribution = false){
 
 Rcpp::List population_list = genetic_results["FinalPopulation"];
 std::vector<std::vector<int>> cocktails = population_list["cocktails"];
 std::vector<double> score = population_list["score"];
 std::vector<double> p_value;
 p_value.resize(cocktails.size(), std::numeric_limits<double>::infinity());
 
 Rcpp::Function compute_p_value = Rf_findFun(Rf_install("p_value_on_sampled"),
                                             R_GlobalEnv);
 
 for(const Rcpp::List& list : distribution_outputs){
    int k = list["cocktailSize"];
    
    for(int i = 0 ; i < cocktails.size(); ++i){
      if(cocktails[i].size() == k){
        p_value[i] = Rcpp::as<double>(compute_p_value(list, score[i], filtred_distribution));
      }
    }
 }
 
 return p_value;
}


//' Used to add the p_value to each cocktail of cocktail list
//' @param distribution_outputs A list of distribution of cocktails of different sizes
//' in order to compute the p_value for multiple cocktail sizes
//' @param cocktails A list containing cocktails in the form of vector of integers (ATC index)
//' @param ATCtree ATC tree with upper bound of the DFS (without the root)
//' @param observations observation of the AE based on the medications of each patients
//'(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
//' on which we want to compute the risk distribution
//' @param filtred_distribution Does the p-values have to be computed using filtered distribution
//' or normal distribution (filtered distribution by default)
//' @param num_thread Number of thread to run in parallel if openMP is available, 1 by default
//' @examples
//' \donttest{
//'  data("ATC_Tree_UpperBound_2024")
//'  data("FAERS_myopathy")
//'  
//'   DistributionApproximationResults_size2 = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy, Smax = 2)
//'             
//'   DistributionApproximationResults_size3 = DistributionApproximation(epochs = 10,
//'             ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy, Smax = 3)
//'             
//'   score_distribution_list = list(DistributionApproximationResults_size2,
//'                               DistributionApproximationResults_size3)
//' 
//'   cocktails = list(c(561, 904),
//'                c(1902, 4585))
//'  
//'   p_value_cocktails(score_distribution_list, cocktails, ATC_Tree_UpperBound_2024,
//'                     FAERS_myopathy)
//' }
//' @return A real valued number vector representing the p-value of the inputed
//' cocktails computed on the distribution_outputs List.
//' @export
//[[Rcpp::export]]
std::vector<double> p_value_cocktails(const std::vector<Rcpp::List>& distribution_outputs, 
                                           const std::vector<std::vector<int>>& cocktails,
                                           const DataFrame& ATCtree,
                                           const DataFrame& observations,
                                           int num_thread = 1,
                                           bool filtred_distribution = false){
 
 std::vector<double> score = compute_hypergeom_on_list(cocktails,
                                                       ATCtree,
                                                       observations,
                                                       num_thread);
 std::vector<double> p_value;
 p_value.resize(cocktails.size(), std::numeric_limits<double>::infinity());
 
 Rcpp::Function compute_p_value = Rf_findFun(Rf_install("p_value_on_sampled"),
                                             R_GlobalEnv);
 
 for(const Rcpp::List& list : distribution_outputs){
   int k = list["cocktailSize"];
   
   for(int i = 0 ; i < cocktails.size(); ++i){
     if(cocktails[i].size() == k){
       p_value[i] = Rcpp::as<double>(compute_p_value(list, score[i], filtred_distribution));
     }
   }
 }
 
 return p_value;
}

//' Function used to convert your genetic algorithm results that are stored into 
//' a .csv file to a Data structure that can be used by the clustering algorithm
//' @param ATC_name the ATC_name column of the ATC tree
//' @param filename Name of the file where the results are located
//' @param sep the separator to use when opening the csv file (';' by default)
//' @return An R List that can be used by other algorithms (e.g. clustering algorithm)
//' @examples
//' \donttest{
//'   data("ATC_Tree_UpperBound_2024")
//'   genetic_results = csv_to_population(ATC_Tree_UpperBound_2024$Name,
//'                     "path/to/output.csv")
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List csv_to_population(const std::vector<std::string>& ATC_name,
                                  const std::string& filename,
                                  const std::string & sep = ";"){
  
  std::vector<std::vector<std::string>> file = check_extension_and_read_csv(filename, sep[0]);
  if(file.size() == 0 ){
    Rcpp::Rcerr << "No cocktail to recover\n";
    return Rcpp::List();
  }
  std::vector<std::vector<int>> cocktails;
  cocktails.reserve(file.size());
  
  for(const auto& row : file){
    std::string drug;
    std::stringstream ss(row[1]);
    
    std::vector<int> rowth_cocktail;
    rowth_cocktail.reserve(7);
    
    while(std::getline(ss, drug, ':')){
      auto it = std::find(ATC_name.begin(), ATC_name.end(), drug);
      if(it != ATC_name.end()){
        rowth_cocktail.push_back(std::distance(ATC_name.begin(), it));
      }
    }
    rowth_cocktail.shrink_to_fit();
    cocktails.push_back(rowth_cocktail);
  }
  
  return Rcpp::wrap(cocktails);
}

//' Function used to convert a string vector of drugs in form "drug1:drug2" to 
//' a vector of index of the ATC tree ex: c(ATC_index(drug1), ATC_index(drugs2))
//' @param ATC_name the ATC_name column of the ATC tree
//' @param lines A string vector of drugs cocktail in the form "drug1:drug2:...:drug_n"
//' @return An R List that can be used by other algorithms (e.g. clustering algorithm)
//' @examples
//' \donttest{
//'   data("ATC_Tree_UpperBound_2024")
//'   string_list = c('hmg coa reductase inhibitors:nervous system',
//'                   'metformin:prasugrel')
//'   string_list_to_int_cocktails(ATC_Tree_UpperBound_2024$Name,
//'                               string_list)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List string_list_to_int_cocktails(const std::vector<std::string>& ATC_name,
                                        const std::vector<std::string>& lines){
  std::vector<std::vector<int>> cocktails;
  cocktails.reserve(lines.size());
  
  for(const auto& line : lines){
    std::string drug;
    std::stringstream ss(line);
    
    std::vector<int> rowth_cocktail;
    rowth_cocktail.reserve(7);
    
    while(std::getline(ss, drug, ':')){
      auto it = std::find(ATC_name.begin(), ATC_name.end(), drug);
      if(it != ATC_name.end()){
        rowth_cocktail.push_back(std::distance(ATC_name.begin(), it));
      } 
    }
    rowth_cocktail.shrink_to_fit();
    cocktails.push_back(rowth_cocktail);
  }
   return Rcpp::wrap(cocktails);
}

//' Function used to convert integer cocktails (like the one outputed by the distributionApproximation function)
//' to string cocktail in order to make them more readable
//' 
//' @param cocktails cocktails vector to be converted (index in the ATC tree)
//' @param ATC_name The ATC_name column of the ATC tree
//' 
//' @return The name of integer cocktails in cocktails
//' @examples
//' \donttest{
//'   data("ATC_Tree_UpperBound_2024")
//'   int_list = list(c(561, 904),
//'                c(1902, 4585))
//'   int_cocktail_to_string_cocktail(int_list, ATC_Tree_UpperBound_2024$Name)
//' }
//' @export
// [[Rcpp::export]]
std::vector<std::vector<std::string>> int_cocktail_to_string_cocktail(
    const std::vector<std::vector<int>>& cocktails, const std::vector<std::string>& ATC_name){
  std::vector<std::vector<std::string>> string_cocktails;
  string_cocktails.reserve(cocktails.size());
  
  for(const auto& cocktail : cocktails){
    std::vector<std::string> current_string_cocktail;
    current_string_cocktail.reserve(cocktail.size());
    for(const auto& med : cocktail)
      current_string_cocktail.push_back(ATC_name[med]);
    
    string_cocktails.push_back(current_string_cocktail);
  }
  
  return string_cocktails;
}



