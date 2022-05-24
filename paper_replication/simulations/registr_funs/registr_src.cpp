// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

//' Calculate expected score and score variance for the current subject.
//'
//' Calculations derived using maximum likelihood estimation.
//'
//' @param Y vector of observations for the current subject.
//' @param mu vector of spline coefficients for the population mean.
//' @param psi matrix of spline coefficients for the principal component basis functions.
//' @param theta spline basis functions for the current subject.
//' @param theta_quad quadratic form of theta for the current subject.
//' @return A list with expected score mean and variance for the current subject.
// [[Rcpp::export]]
List expectedScores(arma::vec Y, arma::vec mu, arma::mat psi, arma::mat theta, arma::mat theta_quad){
  int npc = psi.n_cols;
  arma::mat Ci_inner = arma::eye(npc, npc) - 2 * trans(psi) * theta_quad * psi;
  arma::mat Ci = inv(Ci_inner);
  
  arma::mat mi_inner = trans(Y - 0.5) * theta * psi + 2 * trans(mu) * theta_quad * psi;
  arma::mat mi = Ci * trans(mi_inner);
  
  List result;
  result["mi"] = mi;
  result["Ci"] = Ci;
  return result;
}


//' Estimate variational parameter for the current subject.
//'
//' Function calculates value of variational parameter using maximum likelihood.
//'
//' @param theta spline basis functions for the current subject.
//' @param mu vector of spline coefficients for the population mean.
//' @param mi vector of expected mean scores for the current subject.
//' @param psi matrix of spline coefficients for the principal component basis functions.
//' @param Ci expected covariance matrix of scores for the current subject.
//' @return A vector of variational parameters for the current subject.
// [[Rcpp::export]]
std::vector<double> expectedXi(arma::mat theta, arma::vec mu, arma::vec mi, arma::mat psi, arma::mat Ci){
  int Di = theta.n_rows;
  std::vector<double> xi;
  arma::mat theta_sq, theta_sq_psi;
  double xi_sq;
  
  for(int t = 0; t < Di; t++){
    theta_sq = trans(theta.row(t)) * theta.row(t);
    theta_sq_psi = trans(psi) * theta_sq * psi;
    xi_sq = as_scalar(trans(mu) * theta_sq * mu) + 
      2.0 * as_scalar(trans(mu) * theta_sq * psi * mi) +
      trace(theta_sq_psi * Ci) + as_scalar(trans(mi) * theta_sq_psi * mi );
    
    xi.push_back(sqrt(xi_sq));
  } // end for loop
  
  return(xi);
} 

//' Apply lambda transformation of variational parameter.
//'
//' Simple function for use within other C++ functions.
//'
//' @param x The value to which you apply the function
//' @return A numeric value that has been transformed.
// [[Rcpp::export]]
double lambdaF(double x){
  double y = (0.5- pow(1.0+exp(-x), -1))/2.0/x;
  return y;
}


//' Calculate quadratic form of spline basis functions for the current subject.
//'
//' Calculations quadratic form of theta with diagonalized variational parameter in the center.
//'
//' @param xi vector of variational parameters for the current subject.
//' @param theta spline basis functions for the current subject.
//' @return A matrix of the quadratic form of theta for the current subject.
// [[Rcpp::export]]
arma::mat squareTheta(arma::vec xi, arma::mat theta){
  xi.transform(lambdaF); 
  arma::mat temp_xi_theta = diagmat(xi) * theta;
  arma::mat theta_quad =  trans(theta) * temp_xi_theta;
  
  return(theta_quad);
}