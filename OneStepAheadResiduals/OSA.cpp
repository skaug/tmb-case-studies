// Code from Thygesen et al., 2017 is used.
#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);                  //Observations
  DATA_VECTOR_INDICATOR(keep, y);  //Indicator variable for one-step-ahead predictions
  DATA_SCALAR(huge);    //Used to integreat out X_1
  PARAMETER_VECTOR(x);  //The latent AR1 random variable
  PARAMETER(mu);        //The drift
  PARAMETER(logsigma);  //Log sd of epsilon
  PARAMETER(logs);      //Log sd of psi
  
  // Initial condition
  Type nll = -dnorm(x(0), Type(0), huge, true);
  
  // Increments
  for (int i = 1; i < x.size(); ++i)
    nll -= dnorm(x(i), x(i - 1) + mu, exp(logsigma), true);
  
  // Observations
  for (int i = 0; i < y.size(); ++i)
    nll -= keep(i) * dnorm(y(i), x(i), exp(logs), true);
  
  return nll;
}
