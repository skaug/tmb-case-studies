#include <TMB.hpp>                                // Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;

  DATA_VECTOR(Y);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(S);
  DATA_VECTOR(Sdims);
  DATA_SPARSE_MATRIX(designMatrixForReport);

  PARAMETER(beta0);
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(log_lambda);
  PARAMETER(log_sigma);

  Type sigma = exp(log_sigma);
  vector<Type> lambda = exp(log_lambda);

  Type nll=0;

  vector<Type> S_beta = S*beta;
  nll -= 0.5*(log_lambda*Sdims).sum();
  int counter = 0;
  for(int i=0;i<Sdims.size(); i++){
    for(int j=0;j<Sdims(i); j++){
      nll -= -0.5*lambda(i)*beta(counter)*S_beta(counter);
      counter++;
    }
  }

  vector<Type> mu(Y.size());
  mu = beta0 + X*beta;
  for(int i=0; i<Y.size(); i++){
    nll -= dnorm(Y(i), mu(i), sigma, true);
  }

  vector<Type> splineForReport = designMatrixForReport*beta;
  ADREPORT(splineForReport);
  ADREPORT(beta);

  return nll;
}
