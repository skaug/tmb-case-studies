#include <TMB.hpp>
#include <math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;

  //Load data--------------
  DATA_VECTOR(logPM10); //The response
  DATA_MATRIX(X); //Design matrix for fixed effects
  DATA_INTEGER(n_data); //Number of data points
  DATA_INTEGER(maxDt);//Number of intervalls in the AR1 structure
  DATA_INTEGER(lengthDt); //Length of time intervalls in the AR1 strucutre
  DATA_STRUCT(spdeMatrices,spde_t);//The three matrices from R-INLA which defines the sparce spatial precision structure.
  DATA_SPARSE_MATRIX(A);//Matrix for interpolating points within triangles
  DATA_VECTOR(aLoc); //Help variable for interpoalation of points within triangles
  DATA_INTEGER(flag); // flag=0 => only prior retuned, used when we normalize outside of C++
  //-----------------------

  //lOAD PARAMETERS--------
  PARAMETER_VECTOR(beta); //Regression coefficients
  PARAMETER(log_tau); //Parameter in the Matern covariance function
  PARAMETER(log_kappa);//Parameter in the Matern covariance function
  PARAMETER(rhoTan);//Parameter in the AR1 covariance function
  PARAMETER_ARRAY(x);//The spatio-temporal latent field
  PARAMETER(logSigmaE);//Parameter for unexplained variation
  //------------------------

  //Transform variables-----
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type sigmaE = exp(logSigmaE);
  Type rho = tanh(rhoTan);
  //------------------------

  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q_s = Q_spde(spdeMatrices,kappa);
  //------------------------

  //Calculates nll-------------------------------
  Type nll = 0;

  //Add prior
  nll = GMRF(Q_s,false)(x.col(0)*sqrt(1-rho*rho)); //By some reason the optimization fails if we do not remove the noramlizing constant
  for(int i=1;i<maxDt;i++){
    nll += GMRF(Q_s,false)(x.col(i)-rho*x.col(i-1));
  }
//  nll = SEPARABLE(AR1(rho), GMRF(Q_s))(x); //Alternative approach with use of seperability, includes normalizing constant here because TMB seems to be both faster and more stable when including it.

  if (flag == 0) return nll; // Return un-normalized density on request

  vector<Type> eta = X*beta;

  int counter = 0;
  for(int i = 0; i<maxDt; ++i){
  vector<Type> gammaDt = (A*vector<Type> (x.col(i)));
    for(int j =0; (j<lengthDt & counter <n_data); ++j){
      Type mu;

      mu = eta(counter) + gammaDt(CppAD::Integer(aLoc(j)))/tau;
      if(logPM10(counter)>(-99)){
         nll -= dnorm(logPM10(counter), mu, sigmaE,true);
      }
      counter++;
    }
  }
  //---------------------------------------------

  //Report what we want to report----------------
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  ADREPORT(range);
  //---------------------------------------------

  return nll;
}
