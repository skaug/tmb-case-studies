#include <TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_VECTOR(Cobs); 
  DATA_MATRIX(N);
  DATA_MATRIX(F);
  DATA_MATRIX(M);
  DATA_IMATRIX(aux); // Integer matrix
  DATA_INTEGER(minAge); // Lowest age 
  DATA_INTEGER(minYear); // First observed year
  
  // row i in aux corresponds to element i in Cobs
  
  PARAMETER(log_sigma); 
  
  Type sigma = exp(log_sigma); 
  ADREPORT(sigma); 
  
  // Report uncertainty about prediction 
  int n = Cobs.size(); // Number of observations 
  vector<Type> logPred(n);
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
  
  ADREPORT(logPred);
  
  return nll; 
  
}