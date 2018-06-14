// Debugger (gdb) tutorial
#include <TMB.hpp>
#include "debug_print.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  PARAMETER(a);
  
  X(0,0);
  for(int i=1; i<=y.size(); i++){
    y(i);  // Out-of-range error when i == y.size()
  }
  
  Type nll = a;
  return nll;
}

