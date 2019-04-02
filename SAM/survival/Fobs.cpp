#include<TMB.hpp>

// Helper function to ensure valid correlation 

template<class Type> 
Type trans(Type trans_rho){
  Type rho = (exp(trans_rho) - Type(1)) / (exp(trans_rho) + Type(1));
  return rho;
}

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_MATRIX(Fobs);
  DATA_INTEGER(cormode); 
  
  PARAMETER_MATRIX(logF);
  PARAMETER_VECTOR(log_sigma_logF);
  PARAMETER_VECTOR(trans_rho); 
  PARAMETER(log_sigma_Fobs); 
  
  
  // Number of years and age groups
  int n_year = Fobs.rows();
  int n_age = Fobs.cols();
  
  // Transform parameters 
  Type sigma_Fobs = exp(log_sigma_Fobs); 
  vector<Type> sigma_logF = exp(log_sigma_logF); 
  Type rho = 0;
  
  
  // Report standard deviation of parameters
  ADREPORT(sigma_Fobs);
  ADREPORT(sigma_logF);
  //ADREPORT(rho);
  
  // Covariance matrix for age groups 
  matrix<Type> Sigma(n_age, n_age);
  Sigma.setZero();
  // Negative log likelihood
  Type nll = 0; 
  
  
  switch(cormode){
  
    // Independent
    case 0:
      Sigma.diagonal() = sigma_logF * sigma_logF;
    break;
      
    // Compound symmertry
    case 1:
      Sigma.diagonal() = sigma_logF * sigma_logF;
      rho = trans(trans_rho(0)); 
      for(int i = 0; i < n_age; i++){
        for(int j = 0; j < i; j++){
          Sigma(i, j) = rho * sigma_logF(i) * sigma_logF(j);
          Sigma(j, i) = Sigma(i, j);
        }
      }
    break;
      
    // AR(1) 
    case 2: 
      Sigma.diagonal() = sigma_logF * sigma_logF;
      rho = trans(trans_rho(0)); 
      for(int i = 0; i < n_age; i++){
        for(int j = 0; j < i; j++){
          Sigma(i, j) = pow(rho, Type(i - j)) * sigma_logF(i) * sigma_logF(j);
          Sigma(j, i) = Sigma(i, j);
        }
      }
    break;
      
    default:
      std::cout<<"Error: This cormode not implemented yet."<<std::endl;
      exit(EXIT_FAILURE);
    break;
  }
  
  // Contribution from latent process
  
  // Make zero mean multivariate object with covariance Sigma
  density::MVNORM_t<Type> negative_mvn(Sigma);
  ADREPORT(Sigma);
  
  
  for(int y = 1; y < n_year; y++){
    nll += negative_mvn(logF.row(y) - logF.row(y - 1)); // Returns negative log likelihood  
  }
  
  // Contribution from observations
  for(int y = 0; y < n_year; y++){
    for(int a = 0; a < n_age; a++){
      nll -= dnorm(log(Fobs(y, a)), logF(y, a), sigma_Fobs, true);
    }
  }
  
  return nll;
  
  
  
}