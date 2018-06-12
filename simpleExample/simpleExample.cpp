#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(NL);
  DATA_VECTOR(x);
  
  PARAMETER(beta0);
  PARAMETER(beta1);
  PARAMETER(logSigma);
  
  Type sigma= exp(logSigma);
  Type nll = 0;
  
  vector<Type> mu = beta0 + x*beta1;
  nll += -sum(dnorm(NL,mu,sigma,TRUE));

  return(nll);
}
