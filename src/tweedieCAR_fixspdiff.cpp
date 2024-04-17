#include <TMB.hpp>

using namespace density;

template<class Type>
Type objective_function<Type>::operator() () {
  
  //Data block
  DATA_VECTOR(y);           //observations
  DATA_MATRIX(X);           //covariates
  DATA_FACTOR(group);      // Integer indicator vector for random effects
  DATA_MATRIX(B);          // neighbors matrix
  
  
  // Transformed data block --------------------------------------------------//
  int n = y.size();             // y is the vector of holding observations
  int S = B.cols();
  matrix<Type> D(S, S);
  for(int i = 0; i < S; i++){
    for(int j = 0; j < S; j++){
      if(i == j){
        D(i, j) = B.row(i).sum();
      } else {
        D(i, j) = Type(0.0);
      }
    }
  }
  
  
  // Parameters block -------------------------------------------------------//
  PARAMETER_VECTOR(s);         // spatial random effect
  PARAMETER(log_sigma_s);      // log sd of random effects
  PARAMETER(logit_rho);        // logit of spatial autocorrelation
  PARAMETER_VECTOR(beta);      // regression coefficients
  PARAMETER(log_phi);          // error deviance
  PARAMETER(logit_xi);         // power parameter
  
  
  // Transformed parameters block ----------------------------------------- //
  Type phi = exp(log_phi);
  Type xi = exp(logit_xi)/(Type(1) + exp(logit_xi)) + Type(1);
  Type sigma_s = exp(log_sigma_s);
  Type rho = exp(logit_rho) / (Type(1) + exp(logit_rho));
  vector<Type> eta = X * beta;
  
  // construct precision matrix based on neighbors
  matrix<Type> Q = (D - rho * B);
  Eigen::SparseMatrix<Type> Q_sparse = tmbutils::asSparseMatrix(Q);
  
  
  // Likelihood ---------------------------------------------------------- //
  Type nll = Type(0);
  
  // random effects contribution
  nll += SCALE(GMRF(Q_sparse), sigma_s)(s);
  
  // simulation of random effects
  SIMULATE{
    SCALE(GMRF(Q_sparse), sigma_s).simulate(s);
    REPORT(s);
  }
  
  for(int i=0; i < n; i++){
    // conditional effects
    eta[i] += s[group[i]];
    nll -= dtweedie(y[i], exp(eta[i]), phi, xi, true);
  }
  
  // simulate response for residuals checking
  SIMULATE{
    for(int i = 0; i < n; i++){
      y[i] = rtweedie(exp(eta[i]), phi, xi);
    }
    REPORT(y);
  }
  
  ADREPORT(phi);
  ADREPORT(xi);
  ADREPORT(rho);
  ADREPORT(sigma_s);
  return nll;
}
