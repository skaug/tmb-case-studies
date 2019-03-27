#include <TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_VECTOR(Robs); 
  
  PARAMETER_VECTOR(logR); // Latent process 
  PARAMETER(log_sigma_Robs); 
  PARAMETER(log_sigma_logR); 
  
  // Transform data
  vector<Type> log_Robs = log(Robs);
  
 // Transform standard deviations
 Type sigma_Robs = exp(log_sigma_Robs); 
 Type sigma_logR = exp(log_sigma_logR);
 
 ADREPORT(sigma_Robs); 
 ADREPORT(sigma_logR); 
 
 // Negative log likelihood
 Type nll = 0; 
 // Number of observations 
 Type n = Robs.size(); 
 
 // Contribution to likelihood from latent process logR
 for(int i = 1; i < n; i++){
   nll -= dnorm(logR(i), logR(i - 1), sigma_logR, true);
 }
 
 // Contribution to likelihood from observations Robs
 for(int i = 0; i < n; i++){
   nll -= dnorm(log_Robs(i), logR(i), sigma_Robs, true);
 }
 
 return(nll);
  
}