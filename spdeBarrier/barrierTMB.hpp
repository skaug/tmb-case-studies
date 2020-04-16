
template<class Type>
struct spde_barrier_t{
  vector<Type> C0;
  vector<Type> C1;
  SparseMatrix<Type> D0;
  SparseMatrix<Type> D1;
  SparseMatrix<Type> I;
  spde_barrier_t(SEXP x){  /* x = List passed from R */
    C0 = asVector<Type>(getListElement(x,"C0"));
    C1 = asVector<Type>(getListElement(x,"C1"));
    D0 = asSparseMatrix<Type>(getListElement(x,"D0"));
    D1 = asSparseMatrix<Type>(getListElement(x,"D1"));
    I = asSparseMatrix<Type>(getListElement(x,"I"));
  }
};

template<class Type>
SparseMatrix<Type> Q_spde(spde_barrier_t<Type> spde, Type kappa, vector<Type> c){
  vector <Type> range(2);
  range(0) = sqrt(8)/kappa*c(0);
  range(1) = range(0)*c(1);
  Type pi = 3.141592;

  int dimLatent = spde.D0.row(0).size();
  vector<Type> Cdiag(dimLatent);
  SparseMatrix<Type > Cinv(dimLatent,dimLatent);

  Cdiag = spde.C0*pow(range(0),2) + spde.C1*pow(range(1),2);
  for(int i =0; i<dimLatent; ++i){
    Cinv.coeffRef(i,i) = 1/Cdiag(i);
  }

  SparseMatrix<Type >A = spde.I;
  A = A + (pow(range(0),2)/8) * spde.D0 + (pow(range(1),2)/8) * spde.D1;

  Eigen::SparseMatrix<Type> Q = A.transpose() * Cinv * A/pi *2 * 3;

  return Q;
}

