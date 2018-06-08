#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density; //use GMRF
  using namespace R_inla; //use Q_spde
  using namespace Eigen; //Utilize sparse structures

  DATA_VECTOR(age); //The response
  DATA_SPARSE_MATRIX(X); //Design matrix for splines
  DATA_SPARSE_MATRIX(designMatrixForReport); //Design matrix for report of splines
  DATA_VECTOR(Sdims); //Dimensions of penalization matrices
  DATA_INTEGER(antObs); //Number of observation
  DATA_SPARSE_MATRIX(S); //Penalization matrix
  DATA_STRUCT(spde,spde_t); //Three matrices needed for representing the spatial field
  DATA_SPARSE_MATRIX(A);  //Matrix for interpolating points witin triangles
  DATA_INTEGER(flag); // if flag=0 the prior for x is calculated
  //-----------------------

  PARAMETER_VECTOR(beta0); //Intercepts
  PARAMETER_VECTOR(betaLength); //Regression parameters
  PARAMETER_VECTOR(log_lambda); //Penalization parameter
  PARAMETER_VECTOR(x1);//spatial field age 1
  PARAMETER_VECTOR(x2);//spatial field age 2
  PARAMETER_VECTOR(x3);//spatial field age 3
  PARAMETER_VECTOR(x4);//spatial field age 4
  PARAMETER_VECTOR(x5);//spatial field age 5
  PARAMETER(logKappa);//Spatial scale parameter in Matern covariance structures
  PARAMETER(logTau);//Precision parameter in Matern covariance structures
  //-----------------------

  Type kappa = exp(logKappa);
  Type tau = exp(logTau);
  vector<Type> lambda = exp(log_lambda);

  SparseMatrix<Type> Q1 = Q_spde(spde,kappa);
  SparseMatrix<Type> Q2 = Q_spde(spde,kappa);
  SparseMatrix<Type> Q3 = Q_spde(spde,kappa);
  SparseMatrix<Type> Q4 = Q_spde(spde,kappa);
  SparseMatrix<Type> Q5 = Q_spde(spde,kappa);

  Type nll = 0;

  nll += GMRF(Q1,false)(x1)+
  GMRF(Q2,false)(x2)+
  GMRF(Q3,false)(x3)+
  GMRF(Q4,false)(x4)+
  GMRF(Q5,false)(x5);

  //Construct the penelization precision matrix----------------
  matrix<Type> lambdaDiag(CppAD::Integer(sum(Sdims)),CppAD::Integer(sum(Sdims)));
  lambdaDiag.setZero();
  int counter = 0;
  for(int i =0; i<Sdims.size(); i++){
   for(int j=0; j<Sdims(i); j++){
     lambdaDiag(counter,counter) = lambda(i);
     counter++;
   }
  }
  SparseMatrix<Type> sparse = lambdaDiag.sparseView();
  SparseMatrix<Type> lambdaS = S*sparse;
  //---------------------------------------------------------

  //Add prior of spline regression prameters-----------------
  nll += GMRF(lambdaS,false)(betaLength);
  //---------------------------------------------------------

  //Return prior on request----------------------------------
  if (flag == 0) return nll;
  //---------------------------------------------------------

  //Calculates linear predictors-----------------------------
  vector<Type> ageLength = X* betaLength;

  vector<Type> delta1 = (A*x1)/tau;
  vector<Type> delta2 = (A*x2)/tau;
  vector<Type> delta3 = (A*x3)/tau;
  vector<Type> delta4 = (A*x4)/tau;
  vector<Type> delta5 = (A*x5)/tau;

  vector<Type> nu1 = exp(beta0(0) + ageLength.segment(0,antObs)+  delta1);
  vector<Type> nu2 = exp(beta0(1) + ageLength.segment(antObs,antObs)+  delta2);
  vector<Type> nu3 = exp(beta0(2) + ageLength.segment(2*antObs,antObs)+  delta3);
  vector<Type> nu4 = exp(beta0(3) + ageLength.segment(3*antObs,antObs)+  delta4);
  vector<Type> nu5 = exp(beta0(4) + ageLength.segment(4*antObs,antObs)+  delta5);
  //---------------------------------------------------------

  //Add the likelihood contribution-----------------------------
  vector<Type> prob(7);
  Type sum;
  for(int i=0;i<antObs; ++i){

    sum = nu1(i) + nu2(i) + nu3(i) + nu4(i)+ nu5(i);
    prob(1) = nu1(i)/(1+sum);
    prob(2) = nu2(i)/(1+sum);
    prob(3) = nu3(i)/(1+sum);
    prob(4) = nu4(i)/(1+sum);
    prob(5) = nu5(i)/(1+sum);
    prob(6) = 1/(1+sum);

    if(age(i)==1){
      nll += -log(prob(1));
    }
    if(age(i)==2){
      nll += -log(prob(2));
    }
    if(age(i)==3){
      nll += -log(prob(3));
    }
    if(age(i)==4){
      nll += -log(prob(4));
    }
    if(age(i)==5){
      nll += -log(prob(5));
    }
    if(age(i)>5){
      nll += -log(prob(6));
    }
  }
  //--------------------------------------------------------------------


  //Report the splines for length effect---------------
  vector<Type> repLength = designMatrixForReport* betaLength;
  ADREPORT(repLength);
  //--------------------------------------------------------------------
 return nll;
 }

