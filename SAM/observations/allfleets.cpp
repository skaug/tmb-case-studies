#include <TMB.hpp>


// Function that checks if data is missing 
// Returns logical true or false
template<class Type> 
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type> 
Type objective_function<Type>::operator()(){
  
 //--------------DATA import -------------------------
  DATA_MATRIX(N); // Recruitment and survival
  DATA_MATRIX(F); // Fishing
  DATA_MATRIX(M); // Natural mortality 
  DATA_IMATRIX(aux); // row i corresponds to observation i in obs
  DATA_IMATRIX(idx1); // Not used
  DATA_IMATRIX(idx2); // Not used
  DATA_IVECTOR(fleetTypes); // 0 means catch, 2 means survey 
  DATA_VECTOR(sampleTimes); // Survey times
  DATA_INTEGER(minAge); // Lowest age 
  DATA_INTEGER(minYear); // First observed year
  DATA_VECTOR(obs); // observations accross age and year
  DATA_IMATRIX(keyQ); // Index for proportional constant 
  DATA_IMATRIX(keySd); // Index for observed age group within survey 
  
  
  // ------- PARAMETER import  and transformations----------------------
  PARAMETER_VECTOR(logQ); // Proportional constant for all observed age groups and all surveys 
  PARAMETER_VECTOR(log_sigma_s); // Standard error for a survey
  PARAMETER_VECTOR(missing); // Estimate missing data as a random effect
  
  // Transform log_sigma_s
  vector<Type> sigma_s = exp(log_sigma_s); 
  ADREPORT(sigma_s); 

  
  //-----Define variables used later 
  int n = obs.size(); // Number of observations 
  vector<Type> logPred(n);   // Estimated predictions
  logPred.setZero();
  
  // Get year and age for each observation
  int year, fleet, age = 0; 
  
  // Negative log likelihood
  Type nll = 0;
  
  // Total mortality rate
  Type Z = 0; 

  
 //--------------- Build likelihood function--------------------
 
 
 // Impute missing values 
 
 // count missing values 
 int index_is_missing = 0; 
 for(int i = 0; i < n; i++){
   // Estimate missing data
   if(isNA(obs(i))){
     obs(i) = exp(missing(index_is_missing++));
   }
 }
  
  // Loop over all observations
  for(int i = 0; i < n; i++){
    
    year = aux(i, 0) - minYear;
    fleet = aux(i, 1) - 1; 
    age = aux(i, 2) - minAge; 
    Z = F(year, age) + M(year, age); 
    
    switch(fleetTypes(fleet)){
    case 0: // catch fleet 
      logPred(i) = log(F(year, age)) - log(Z) + log(N(year, age)) + log(Type(1) - exp(-Z));  
    break; 
    
    case 2: // survey fleet
      logPred(i) = logQ(keyQ(fleet, age)) + log(N(year, age)) - sampleTimes(fleet) * Z; 
    break; 
      
    default:
      std::cout << "This fleet is not implemented yet. Input must be 0 or 2." << std::endl; 
      exit(EXIT_FAILURE);
    break; 
    }
    
    nll -= dnorm(log(obs(i)), logPred(i), sigma_s(keySd(fleet, age)), true);
  }
  
  REPORT(logPred);
  return nll; 
  
}
