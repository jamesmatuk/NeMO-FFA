#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;

/* *********** POST PROCESSING FUNCTIONS *********** 
* The following functions resolve label/sign switching ambiguity,
* and return the sign and column indices used to align samples.
* forked from the infinitefactor R package (https://github.com/poworoznek/infinitefactor)
*/

// [[Rcpp::export]]
Rcpp::NumericMatrix msf(arma::mat lambda, arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  arma::uvec ind = regspace<arma::uvec> (0, k-1);
  arma::uvec perm(k);
  arma::vec signs(k);
  arma::rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff, diffsq;
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    diffsq = square(diff);
    norms = sum(diffsq);
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
  }
  
  arma::mat permmat = zeros<mat>(k,k);
  for(int i=0; i<k; i++){
    permmat(perm(i), i) = signs(i);
  }
  
  lambda *= permmat;
  return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericVector msfOUT_functional(arma::mat lambda, arma::mat W,arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  arma::uvec ind = regspace<arma::uvec> (0, k-1);
  arma::uvec perm(k);
  arma::vec signs(k);
  arma::rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff;
  arma::mat temp_norm; 
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    norms = sqrt(diagvec(diff.t()*W*diff).t()); 
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
  }
  
  arma::vec out = (perm + ones<vec>(k)) % signs;
  
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
Rcpp::NumericVector msfOUT(arma::mat lambda, arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  arma::uvec ind = regspace<arma::uvec> (0, k-1);
  arma::uvec perm(k);
  arma::vec signs(k);
  arma::rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff, diffsq;
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    diffsq = square(diff);
    norms = sum(diffsq);
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
  }
  
  arma::vec out = (perm + ones<vec>(k)) % signs;
  
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix aplr(arma::mat matr, arma::vec perm){
  int k = matr.n_cols;
  arma::mat permmat = zeros<arma::mat>(k,k);
  arma::vec perms = abs(perm) - ones<arma::vec>(k);
  arma::vec signs = sign(perm);
  
  for(int i=0; i<k; i++){
    permmat(perms(i), i) = signs(i);
  }
  
  matr *= permmat;
  return(Rcpp::wrap(matr));
}