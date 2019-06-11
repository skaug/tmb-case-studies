#include <TMB.hpp>


// Function that calculates spawning stock biomass (ssb). 
// Paramereters: 
// logN: Survival process 
// F: Fishing 
// M: Mortality 
// SW: Stock mean weight 
// MO: Fraction of stock that are mature
// PF: Fraction applied before spawning. Applied to fishing process
// PM: Fraction applied before spawning. Applied to mortality process
template<class Type> 
vector<Type> ssbF(matrix<Type> logN, matrix<Type> F, matrix<Type> M, 
                  matrix<Type> SW, matrix<Type> MO, matrix<Type> PF,
                  matrix<Type> PM){
  
  
  // Number of years and age groups
  int n_year = logN.rows();
  int n_age = logN.cols(); 
  
  vector<Type> ssb(n_year);
  ssb.setZero();
  
  // sum over all age groups for each year
  for(int y  = 0; y < n_year; y++){
    for(int a = 0; a < n_age; a++){
      ssb(y) += MO(y, a) * SW(y, a) * exp(logN(y, a)) * exp(- PF(y, a) * F(y, a) - PM(y, a) * M(y, a));
    }
  }
  return ssb;
}

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_MATRIX(Nobs); 
  DATA_MATRIX(M); 
  DATA_MATRIX(F); 
  DATA_MATRIX(SW);
  DATA_MATRIX(MO); 
  DATA_MATRIX(PF);
  DATA_MATRIX(PM);
  DATA_INTEGER(method);
  

  PARAMETER_MATRIX(logN);
  PARAMETER(log_sigma_Nobs); 
  PARAMETER(log_sigma_logN);
  PARAMETER(log_sigma_logR); 
  PARAMETER_VECTOR(ricker);
  PARAMETER_VECTOR(bh);
  
  Type sigma_Nobs = exp(log_sigma_Nobs); 
  Type sigma_logN = exp(log_sigma_logN); 
  Type sigma_logR = exp(log_sigma_logR);
  
  ADREPORT(sigma_logN);
  ADREPORT(sigma_Nobs);
  ADREPORT(sigma_logR);
  

  // Negative log likelihood
  Type nll = 0; 
  
  //Number of years and age groups
  int n_year = Nobs.rows();
  int n_age = Nobs.cols(); 
  
  // Prediction for logN
  Type pred = 0; 
  
  // Calculate ssb 
  vector<Type> ssb = ssbF(logN, F, M, SW, MO, PF, PM);
  
  
  // Recruitment for age zero
  for(int y = 1; y < n_year; y++){
    
    switch(method){
    case 0:
      pred = logN(y - 1, 0);
    break;
    
    // Ricker   
    case 1: 
      pred = ricker(0) + log(ssb(y - 1)) - exp(ricker(1)) * ssb(y - 1);
    break;
      
    // Berverton-Holt  
    case 2:
      pred = bh(0) + log(ssb(y - 1)) - log(Type(1.0) + exp(bh(1)) * ssb(y - 1));
    break;
      
    default:
      std::cout << "This method is not implementet. Method has to be a number 0, 1, 2" << std::endl;
    break;
    }
    
    nll -= dnorm(logN(y, 0), pred, sigma_logR, true);
    
  }
  
  // survival for age greater than zero
  
  // Contribution form latent process
  for(int y = 1; y < n_year; y++){
    for(int a = 1; a < n_age; a++){

      // If we are at max age we get contribution from same group last year
      if(a == (n_age - Type(1))){
        Type pred =  log(exp(logN(y - 1, a - 1) - F(y - 1, a - 1) - M(y - 1, a - 1))) +
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
  
  // Calculate standard error for ssb
  ADREPORT(ssb);
  
  return nll;
  
}
