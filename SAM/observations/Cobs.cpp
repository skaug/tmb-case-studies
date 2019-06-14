#include <TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_VECTOR(Cobs); // Catch-at-age
  DATA_MATRIX(N); // Stocksize
  DATA_MATRIX(F); // Fishing mortality
  DATA_MATRIX(M); // Natural mortality
  DATA_IMATRIX(aux); // Integer matrix with year and age for each observation
  DATA_INTEGER(minAge); // Lowest age 
  DATA_INTEGER(minYear); // First observed year
  
  // Note: row i in aux corresponds to element i in Cobs
  
  PARAMETER(log_sigma); 
  
  Type sigma = exp(log_sigma); 
  ADREPORT(sigma); 
  
  int n = Cobs.size(); // Number of observations 
  vector<Type> logPred(n); // Predictions
  logPred.setZero();
  
  // Get year and and for each observation
  int year, age = 0; 
  
  // Negative log likelihood
  Type nll = 0;
  
  // Total mortality rate
  Type Z = 0; 
  
  // Loop over all observations
  for(int i = 0; i < n; i++){
  
    year = aux(i, 0) - minYear; 
    age = aux(i, 2) - minAge; 
    Z = F(year, age) + M(year, age); 
    logPred(i) = log(F(year, age)) - log(Z) + log(N(year, age)) + log(Type(1) - exp(-Z)); 
    
    nll -= dnorm(log(Cobs(i)), logPred(i), sigma, true);
  }
  
  // Report uncertainty about prediction 
  ADREPORT(logPred);
  
  return nll; 
  
}
