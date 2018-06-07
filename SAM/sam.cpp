//  --------------------------------------------------------------------------
  // Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>,
// Casper Berg <cbe@aqua.dtu.dk>, and Kasper Kristensen <kkr@aqua.dtu.dk>.
// All rights reserved.
//
  // Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
  //   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN, CASPER BERG OR KASPER
// KRISTENSEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
                                        // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
                                        // OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
                                                   // OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------

#include <TMB.hpp>
#include <iostream>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type>
Type square(Type x){return x*x;}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(fleetTypes);    //Indicator variables for data types
  DATA_VECTOR(sampleTimes);   //Time in year the time series is collected
  DATA_INTEGER(nobs);         //Number of observations
  DATA_ARRAY(obs);            //Observations from every year, fleet and age.
  DATA_ARRAY(propMat);        //Proportion mature
  DATA_ARRAY(stockMeanWeight);//Estimated mean weight in population
  DATA_ARRAY(catchMeanWeight);//Estimated mean weight in catch
  DATA_ARRAY(natMor);         //Natural mortality
  DATA_INTEGER(minAge);       //Minimum age
  DATA_INTEGER(maxAge);       //Maximum age
  DATA_INTEGER(maxAgePlusGroup); //Indicator variable stating if we have defined a pluss group
  DATA_IARRAY(keyLogFsta);    //Used for defining fishing mortality equal for several ages
  DATA_INTEGER(corFlag);      //Variable defining correlation strucutre for F
  DATA_IARRAY(keyLogFpar);    //Variable used for catchability
  DATA_IARRAY(keyVarF);       //Variable used for setting variances in F equal
  DATA_IARRAY(keyVarLogN);    //Variable used for setting variances in \eta equal
  DATA_IARRAY(keyVarObs);     //Variable used for setting variances in observations equal
  DATA_INTEGER(stockRecruitmentModelCode); //Variable used for choosing recrutiment model
  DATA_INTEGER(noScaledYears);    //Number of years we just know the scaled catch
  DATA_IVECTOR(keyScaledYears);   //Which year do we just know the scaled catch
  DATA_IMATRIX(keyParScaledYA);   //Variable used for estimating scaling
  DATA_IVECTOR(fbarRange);        //Age range of fbar

  PARAMETER_VECTOR(logFpar);      //Catchaility (Q_s^a)
  PARAMETER_VECTOR(logSdLogFsta); //Log standatd deviation of F
  PARAMETER_VECTOR(logSdLogN);    //Log standatd deviation of \eta
  PARAMETER_VECTOR(logSdLogObs);  //Log standatd deviation of observations
  PARAMETER_VECTOR(rec_loga); //Recruitment parameter
  PARAMETER_VECTOR(rec_logb); //Recruitment parameter
  PARAMETER(logit_rho);       //Used in correlation structure for F
  PARAMETER_VECTOR(logScale); //Scaling parameters in those year we just know the scaled catch
  PARAMETER_ARRAY(logF);      //Random effect for F
  PARAMETER_ARRAY(logN);      //Random effect for \eta

  //Define some variables needed--------------------
  int timeSteps=logF.dim[1];
  int stateDimF=logF.dim[0];
  int stateDimN=logN.dim[0];
  Type rho=f(logit_rho);
  vector<Type> sdLogFsta=exp(logSdLogFsta);
  vector<Type> varLogN=exp(logSdLogN*Type(2.0));
  vector<Type> varLogObs=exp(logSdLogObs*Type(2.0));
  vector<Type> ssb(timeSteps);
  vector<Type> logssb(timeSteps);
  vector<Type> fbar(timeSteps);
  vector<Type> logfbar(timeSteps);
  vector<Type> cat(catchMeanWeight.dim(0));
  vector<Type> logCatch(catchMeanWeight.dim(0));
  vector<Type> tsb(timeSteps);
  vector<Type> logtsb(timeSteps);
  vector<Type> logR(timeSteps);
  vector<Type> R(timeSteps);

  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);
  //--------------------------------------


  Type nll=0; //negative log-likelihood


  //Calculates likelihood contribution by F-process----------
  for(int i=0; i<stateDimF; ++i){
    fcor(i,i)=1.0;
  }

  if(corFlag==1){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=rho;
        fcor(j,i)=fcor(i,j);
      }
    }
  }
  if(corFlag==2){
    for(int i=0; i<stateDimF; ++i){
      for(int j=0; j<i; ++j){
        fcor(i,j)=pow(rho,abs(Type(i-j)));
        fcor(j,i)=fcor(i,j);
      }
    }
  }

  int i,j;
  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimN; ++j){
      if(keyLogFsta(0,j)==i)break;
    }
    fsd(i)=sdLogFsta(keyVarF(0,j));
  }

  for(i=0; i<stateDimF; ++i){
    for(j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }

  using namespace density;
  MVNORM_t<Type> neg_log_densityF(fvar);
  for(int i=1;i<timeSteps;i++){
    nll+=neg_log_densityF(logF.col(i)-logF.col(i-1)); // F-Process likelihood
  }
  //---------------------------------------------------------

  //Calculates the SSB---------------------------------------
  for(int i=0;i<timeSteps;i++){ // calc ssb
    ssb(i)=0.0;
    for(int j=0; j<stateDimN; ++j){
        ssb(i)+=exp(logN(j,i))*propMat(i,j)*stockMeanWeight(i,j);
    }
    logssb(i)=log(ssb(i));
  }
  //----------------------------------------------------------

  //Calculates likelihood contribution from \eta--------------
  matrix<Type> nvar(stateDimN,stateDimN);
  for(int i=0; i<stateDimN; ++i){
    for(int j=0; j<stateDimN; ++j){
      if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(keyVarLogN(0,i));}
    }
  }
  MVNORM_t<Type> neg_log_densityN(nvar);
  vector<Type> predN(stateDimN);
  for(int i=1;i<timeSteps;i++){
    if(stockRecruitmentModelCode==0){ // straight RW
      predN(0)=logN(0,i-1);
    }else{
      if(stockRecruitmentModelCode==1){//ricker
        predN(0)=rec_loga(0)+log(ssb(i-1))-exp(rec_logb(0))*ssb(i-1);
      }else{
        if(stockRecruitmentModelCode==2){//BH
          predN(0)=rec_loga(0)+log(ssb(i-1))-log(1.0+exp(rec_logb(0))*ssb(i-1));
        }else{
          error("SR model code not recognized");
        }
      }
    }

    for(int j=1; j<stateDimN; ++j){
      if(keyLogFsta(0,j-1)>(-1)){
        predN(j)=logN(j-1,i-1)-exp(logF(keyLogFsta(0,j-1),i-1))-natMor(i-1,j-1);
      }else{
        predN(j)=logN(j-1,i-1)-natMor(i-1,j-1);
      }
    }
    if(maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF(keyLogFsta(0,stateDimN-2),i-1))-natMor(i-1,stateDimN-2))+
                               exp(logN(stateDimN-1,i-1)-exp(logF(keyLogFsta(0,stateDimN-1),i-1))-natMor(i-1,stateDimN-1)));
    }
    nll+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood
  }
  //----------------------------------------------------------------

  //Calculates the likelihood contribution from the observations----
  int f, ft, a, y,yy, scaleIdx;
  int minYear=CppAD::Integer((obs(0,0)));
  Type zz;
  vector<Type> predObs(nobs);
  vector<Type> predSd(nobs);
  for(int i=0;i<nobs;i++){
    y=CppAD::Integer(obs(i,0))-minYear;
    f=CppAD::Integer(obs(i,1));
    ft=CppAD::Integer(fleetTypes(f-1));
    a=CppAD::Integer(obs(i,2))-minAge;
    zz=exp(logF(keyLogFsta(0,a),y))+natMor(y,a);

    switch(ft){
      case 0:
        predObs(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(keyLogFsta(f-1,a)>(-1)){
          predObs(i)+=logF(keyLogFsta(0,a),y);
        }
        scaleIdx=-1;
        yy=CppAD::Integer(obs(i,0));
        for(int j=0; j<noScaledYears; ++j){
          if(yy==keyScaledYears(j)){
            scaleIdx=keyParScaledYA(j,a);
            if(scaleIdx>=0){
              predObs(i)-=logScale(scaleIdx);
            }
            break;
          }
        }
        break;

        case 2:
          predObs(i)=logN(a,y)-zz*sampleTimes(f-1);
          if(keyLogFpar(f-1,a)>(-1)){
            predObs(i)+=logFpar(keyLogFpar(f-1,a));
          }

        break;

        default:
          std::cerr<<"Unknown fleet code: "<<ft<<std::endl;
        return 0 ;
        break;
    }
    predSd(i)=sqrt(varLogObs(keyVarObs(f-1,a)));
    nll+=-dnorm(log(obs(i,3)),predObs(i),predSd(i),true);
  }
  //-----------------------------------------------------------------

  //Calculates the average fishery mortality----
  for(int y=0;y<timeSteps;y++){
    fbar(y)=Type(0);
    for(int a=fbarRange(0);a<=fbarRange(1);a++){
      fbar(y)+=exp(logF(keyLogFsta(0,a-minAge),y));
    }
    fbar(y)/=Type(fbarRange(1)-fbarRange(0)+1);
    logfbar(y)=log(fbar(y));
  }
  //----------------------------------------------------------------

  //Calculates the estimated catch------------
  for(int y=0;y<catchMeanWeight.dim(0);y++){
    cat(y)=Type(0);
    for(int a=minAge;a<=maxAge;a++){
      Type z=exp(logF(keyLogFsta(0,a-minAge),y))+natMor(y,a-minAge);
      cat(y)+=exp(logF(keyLogFsta(0,a-minAge),y))/z*exp(logN(a-minAge,y))*(Type(1.0)-exp(-z))*catchMeanWeight(y,a-minAge);
    }
    logCatch(y)=log(cat(y));
  }
  //----------------------------------------------------------------

  //Calculates the recruitment--------------------------------------
  for(int y=0;y<timeSteps;y++){
    logR(y)=logN(0,y);
    R(y)=exp(logR(y));
  }
  //----------------------------------------------------------------

  ADREPORT(logssb);
  ADREPORT(logfbar);
  ADREPORT(logCatch);
  ADREPORT(logR);

  return nll;
}

