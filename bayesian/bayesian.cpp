#include <TMB.hpp>

// Hand-coded Cauchy distribution
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(Ndata);
  DATA_INTEGER(Nstage);
  DATA_INTEGER(Nyear);
  DATA_INTEGER(Nplant);
  DATA_IVECTOR(year);
  DATA_IVECTOR(plant);
  DATA_IVECTOR(stage);
  DATA_VECTOR(Pods);
  DATA_IVECTOR(toF);

  // fixed effects -- bounds added in R
  PARAMETER(yearInterceptSD);
  PARAMETER(plantInterceptSD);
  PARAMETER(plantSlopeSD);
  PARAMETER_VECTOR(intercept);
  PARAMETER(slope);

  // non-centered random effects
  PARAMETER_VECTOR(yearInterceptEffect_raw);
  PARAMETER_VECTOR(plantInterceptEffect_raw);
  PARAMETER_VECTOR(plantSlopeEffect_raw);

  Type nlp=0.0; // negative log prior
  Type nll=0.0; // negative log likelihood

  // Postive transformations, jacobians below
  Type yearInterceptSD2=exp(yearInterceptSD);
  Type plantInterceptSD2=exp(plantInterceptSD);
  Type plantSlopeSD2=exp(plantSlopeSD);

  // priors
  nlp-= dcauchy(yearInterceptSD2, Type(0), Type(5));
  nlp-= dcauchy(plantInterceptSD2, Type(0), Type(5));
  nlp-= dcauchy(plantSlopeSD2, Type(0), Type(5));
  nlp-= dnorm(slope, Type(0.0), Type(10.0));
  nlp-= dnorm(intercept, Type(0.0), Type(10.0)).sum();

  Type ypred;
  // model predictions
  for(int i=0; i<Ndata; i++){
    // prediction logit scale
    ypred= intercept(stage(i)-1) +
      yearInterceptEffect_raw(year(i)-1)*yearInterceptSD2 +
      plantInterceptEffect_raw(plant(i)-1)*plantInterceptSD2+
      Pods(i) * plantSlopeEffect_raw(plant(i)-1)*plantSlopeSD2+
      Pods(i) * slope;
    // likelihood contribution
    if(toF(i)==1){
      nll+= log(1+exp(-ypred));
    } else {
      nll+= ypred+log(1+exp(-ypred));
    }
  }

  // random effects; non-centered
  nll-=dnorm(yearInterceptEffect_raw, Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(plantInterceptEffect_raw,Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(plantSlopeEffect_raw, Type(0.0), Type(1.0), true).sum();

  // Jacobian adjustments
  nll-= yearInterceptSD + plantInterceptSD + plantSlopeSD;
  Type nld=nll+nlp; // negative log density
  return(nld);
}
