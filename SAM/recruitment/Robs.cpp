#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(year)
  DATA_VECTOR(ssb)
  DATA_VECTOR(Robs)
  DATA_INTEGER(minAge)    
  DATA_INTEGER(mode)
  
  PARAMETER(logsdo);
  PARAMETER(logsdp);
  PARAMETER_VECTOR(logR);
  PARAMETER_VECTOR(rickerpar);
  PARAMETER_VECTOR(bhpar);
  Type sdo = exp(logsdo);
  Type sdp = exp(logsdp);
  
  Type jnll=0;
  
  Type pred;
  for(int i=1; i<logR.size(); ++i){
    switch(mode){
    case 0:
      pred = logR(i-1);
      break;
      
    case 1:
      pred = rickerpar(0)+log(ssb(i-minAge))-exp(rickerpar(1))*ssb(i-minAge);
      break;
      
    case 2:
      pred = bhpar(0)+log(ssb(i-minAge))-log(1.0+exp(bhpar(1))*ssb(i-minAge));
      break;
      
    default:
      std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
    break;
    }
    jnll += -dnorm(logR(i),pred,sdp,true);
  }
  for(int i=0; i<Robs.size(); ++i){
    jnll += -dnorm(log(Robs(i)),logR(i),sdo,true);
  }
  return jnll;
}