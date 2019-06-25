#include<TMB.hpp>


template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_VECTOR(x); // Observations
  DATA_INTEGER(n); // number of observations
  
  PARAMETER(mu); // parameter to estimate
  
  Type nll = 0; // negative log likelihood 
  
  // loop over observations
  for(int i = 0; i < n; i++){
    nll -= dnorm(x(i), mu, Type(1), true);
  }
  
  // return negative log likelihood
  return nll; 
}