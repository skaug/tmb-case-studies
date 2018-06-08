#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_IVECTOR(i);
  DATA_IVECTOR(j);
    
  PARAMETER_VECTOR(beta);      
  PARAMETER_VECTOR(u);      
  PARAMETER(sigma);      
  PARAMETER(tau);      
  
  Type nll = 0.0;

  nll -= sum(dnorm(u,0,tau,true));
  for(int k=0; k<y.size(); k++){
    Type eta = beta(i(k)) + u(j(k));
    nll -= dnorm(y(k),eta,sigma,true);
  }

  return nll;
}
