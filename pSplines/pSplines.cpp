#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;

  //Read data from R------------------
  DATA_VECTOR(Y);       //The response
  DATA_MATRIX(X);       //Design matrix for splines
  DATA_SPARSE_MATRIX(S);//Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR(Sdims);   //Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX(designMatrixForReport);//Design matrix for report of splines
  DATA_INTEGER(flag);     //Used for calculating the normalizing constant of the prior for beta on R side
  //----------------------------------
  
  //Read parameters from R------------
  PARAMETER(beta0);       //Intercept
  PARAMETER_VECTOR(beta); //Spline regression parameters
  PARAMETER_VECTOR(log_lambda);//Penalization parameters
  PARAMETER(log_sigma);   
  //----------------------------------
  
  //Transform some of the parameters--
  Type sigma = exp(log_sigma);
  vector<Type> lambda = exp(log_lambda);
  //----------------------------------
  
  //Calculate the objective function--
  Type nll=0;


  int k=0;  // Counter
  for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    SparseMatrix<Type> S_i = S.block(k,k,k+m_i,k+m_i);  // Recover Si
    vector<Type> beta_i = beta.segment(k,k+m_i);       // Recover betai
    nll += -0.5*m_i*log_lambda(i) + lambda(i)*GMRF(S_i,false)(beta_i);
    k++;
  }
  
  vector<Type> mu(Y.size());
  mu = beta0 + X*beta;
  for(int i=0; i<Y.size(); i++){
    nll -= dnorm(Y(i), mu(i), sigma, true);
  }
  //----------------------------------
  
  vector<Type> splineForReport = designMatrixForReport*beta;
  ADREPORT(splineForReport);
  ADREPORT(beta);

  return nll;
}
