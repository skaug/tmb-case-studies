#include <TMB.hpp>

template<class Type> 
Type objective_function<Type>::operator()(){
  
  DATA_VECTOR(Robs);
  DATA_VECTOR(SSB);
  DATA_INTEGER(method);
  
  PARAMETER_VECTOR(logR); // Latent process 
  PARAMETER(log_sigma_Robs); 
  PARAMETER(log_sigma_logR); 
  PARAMETER_VECTOR(ricker); 
  PARAMETER_VECTOR(bh); 
  
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
 Type pred = 0; 
 // Contribution to likelihood from latent process logR
 
 
 for(int i = 1; i < n; i++){
   
   switch(method){
   case 0:
     pred = logR(i - 1);
     break;
     
   case 1:
     pred = ricker(0) + log(SSB(i - 1)) - exp(ricker(1)) * SSB(i - 1);
     break;

   case 2:
     pred = bh(0) + log(SSB(i - 1)) - log(Type(1.0) + exp(bh(1)) * SSB(i - 1));
     break;

   default:
     std::cout << "This method is not implementet. Method has to be a number 0, 1, 2" << std::endl;
   break;
   }
   
   nll -= dnorm(logR(i), pred, sigma_logR, true);

 }
 
 // Contribution to likelihood from observations Robs
 for(int i = 0; i < n; i++){
   nll -= dnorm(log_Robs(i), logR(i), sigma_Robs, true);
 }
 
 return(nll);
  
}