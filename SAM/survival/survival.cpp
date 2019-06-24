#include <TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_MATRIX(Nobs); 
  DATA_MATRIX(M); 
  DATA_MATRIX(F); 
  
  PARAMETER_MATRIX(logN); // Survival matrix for all age groups and years 
  PARAMETER(log_sigma_Nobs); // standard deviation for observations 
  PARAMETER(log_sigma_logN); // standard devitation for state equation
  PARAMETER(log_sigma_logR); // stadard deviation for recruitment 
  
  // Transform parameters
  Type sigma_Nobs = exp(log_sigma_Nobs); 
  Type sigma_logN = exp(log_sigma_logN); 
  Type sigma_logR = exp(log_sigma_logR);
  
  // Report to R
  ADREPORT(sigma_logN);
  ADREPORT(sigma_Nobs);
  ADREPORT(sigma_logR);
  

  
  Type nll = 0; 
  int n_year = Nobs.rows();
  int n_age = Nobs.cols(); 
  
  // Recruitment for age zero
  for(int y = 1; y < n_year; y++){
    nll -= dnorm(logN(y, 0), logN(y - 1, 0), sigma_logR, true);
  }
  
  // survival for age greater than zero
  
  // Contribution form latent process
  for(int y = 1; y < n_year; y++){
    for(int a = 1; a < n_age; a++){
      
      // If we are at max age we get contribution from same group last year
      if(a == n_age - Type(1)){
        Type pred = log(exp(logN(y - 1, a - 1) - F(y - 1, a - 1) - M(y - 1, a - 1))) +
                     log(exp(logN(y - 1, a) - F(y - 1, a) - M(y - 1, a)));
        nll -= dnorm(logN(y, a), pred, sigma_logN, true);
      } else{
        nll -= dnorm(logN(y ,a), logN(y - 1, a - 1) - F(y - 1, a - 1) - M(y - 1, a - 1), sigma_logN, true);
      }
    }
  }
  
  // Contribution from observations given latent process
  
  for(int y = 0; y < n_year; y++){
    for(int a = 0; a < n_age; a++){
     nll -= dnorm(log(Nobs(y, a)), logN(y, a), sigma_Nobs, true);
    }
  }
  
  return nll;
  
}
