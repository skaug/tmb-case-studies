#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density; 
  using namespace Eigen; //Needed for utilisation of sparse structures
  
  //Load data and parameters----------------
  DATA_VECTOR(time);      //The response
  DATA_IVECTOR(notcens);  //indicator vector stating which persons that were censored
  DATA_MATRIX(X);         //Design matrix for fixed effects
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles 
  //DATA_INTEGER(flag); // if flag=0 the prior for x is calculated
  
  PARAMETER_VECTOR(beta);  
  PARAMETER(log_tau);
  PARAMETER(log_kappa);
  PARAMETER(log_omega);  
  PARAMETER_VECTOR(x);  
  //---------------------------------------
  
  //Transform parameters-------------------
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type omega = exp(log_omega);  // Parameter of Weibull distribution
  //------------------------------------------

  // Spatial interpolation
  vector<Type> delta = (A*x)/tau;
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_spde(spdeMatrices,kappa);
  //---------------------------------------------
  
  //Calculates nll-------------------------------
  Type nll = 0.0;
  nll = GMRF(Q)(x);                              
  
  //nll = GMRF(Q, false)(x);  
  // Return un-normalized density on request
  //if (flag == 0) return nll;
  
  vector<Type> eta = X*beta + delta;
  for(int i=0; i<time.size(); i++){    
    Type lambda = exp(eta(i));
    Type t_omega = pow(time(i),omega);
    Type S = exp(-lambda*t_omega);        // Survival function
    Type f = lambda*omega*t_omega/time(i)*S;// Weibull density
    if(notcens(i)){
      nll -= log(f);
    }else{
      nll -= log(S); //The pasient survived until cencoring
    }
  }
  //---------------------------------------------
  
  //Report what we want to report----------------
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  ADREPORT(range);
  //---------------------------------------------
  
  return nll;
}
