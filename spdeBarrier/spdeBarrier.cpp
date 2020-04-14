#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
  using namespace density;
  using namespace Eigen; //Needed for utilisation of sparse structures

  //Load data and parameters----------------
  DATA_VECTOR(y);      //The response
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_STRUCT(spdeMatricesBarrier,spde_b); //Structure needed for the barrier procedure
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles
  DATA_INTEGER(barrier); // if barrier==1, barrier procedure is used
  DATA_VECTOR(c);      // Scaling of range
  DATA_MATRIX(X);         //Design matrix for fixed effects
  
  

  PARAMETER_VECTOR(beta); //Regresion coefficients for intercept and covariates
  PARAMETER(log_tau); //log precision for spatial effect
  PARAMETER(log_kappa); //log spatial scale parameter in spatial effect
  PARAMETER(log_sigmaIID); //log standard deviation in iid effect
  PARAMETER_VECTOR(x); //Spatial effect
  PARAMETER_VECTOR(xiid);  //Iid effect
  
  //---------------------------------------

  //Transform parameters-------------------
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type sigmaIID = exp(log_sigmaIID);
  //------------------------------------------

  // Spatial interpolation
  vector<Type> delta = (A*x)/tau;
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q;
  if(barrier==0){
    Q = Q_spde(spdeMatrices,kappa);
  }else if(barrier ==1){
    Q = Q_spde(spdeMatricesBarrier,kappa,c);
  }
  //---------------------------------------------

  //Calculates nll-------------------------------
  Type nll = 0.0;
  nll += GMRF(Q)(x);

  for(int i =0; i< y.size(); ++i){
    nll-= dnorm(xiid(i), Type(0), sigmaIID,true);
  }

  vector<Type> eta = exp(X*beta + delta + xiid);
  for(int i=0; i<y.size(); i++){
    nll -= dpois(y(i), eta(i),true);
  }
  //---------------------------------------------

  //Report what we want to report----------------
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  ADREPORT(range);
  //---------------------------------------------

  return nll;
}
