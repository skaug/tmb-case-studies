#include <TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator()(){
  
  // Data 
  DATA_VECTOR(Robs);
  
  // Parameter
  PARAMETER(log_sigma_r);

  Type sigma_r = exp(log_sigma_r);
  ADREPORT(sigma_r); // Report standard deviation of sigma_r back to R
  
  Type nll = 0; // negative log likelihood
  vector<Type> log_Robs = log(Robs); // log of recruitment
  
  for(int i = 1; i < Robs.size(); i++){
    nll -= dnorm(log_Robs(i), log_Robs(i - 1), sigma_r, true);
  }
  
  return(nll); 
}
