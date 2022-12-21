#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <truncnorm.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
arma::mat make_cov(const arma::vec grid,const double length_param,const double scale_param) {
  arma::vec one_vec(arma::size(grid),arma::fill::ones);
  return scale_param*exp(-pow(one_vec*grid.t() - grid*one_vec.t(),2)/(2*length_param)) ;
}

// [[Rcpp::export]]
arma::mat compute_y_st(const arma::mat y_common,const arma::mat f,const arma::mat obs_grid) {
  return y_common - f%obs_grid;
}

// [[Rcpp::export]]
arma::mat sample_mu_sparse(const arma::mat y_st, const double sigma_sq_current,
                           const arma::vec grid, const double l_param_cur, const double scale_param_cur,
                           const arma::mat obs_grid) {
  
  arma::vec mu_cond_mean;
  arma::mat covariance_matrix = make_cov(grid,l_param_cur,1);
  arma::mat mu_cond_var = covariance_matrix;
  arma::mat eye_m = arma::eye(arma::size(mu_cond_var));
  
  mu_cond_var = arma::inv(arma::inv(covariance_matrix + eye_m*1e-10)/scale_param_cur
                                  + arma::diagmat(arma::sum(obs_grid,1))/sigma_sq_current);
  mu_cond_mean = mu_cond_var*(arma::sum(y_st,1)/sigma_sq_current);
  
  return arma::mvnrnd(mu_cond_mean,mu_cond_var);
}

// [[Rcpp::export]]
double evaluate_mu_cov(arma::mat y_st, const arma::mat covariance_mat, const double sigma_sq_cur,const arma::mat obs_grid) {
  
  arma::vec ret;
  arma::uvec t_i;
  arma::mat marginal_cov;
  double m_i;
  arma::mat I_m_i;
  arma::vec y_st_i;
  arma::mat temp_mat_mult;
  
  for(int i = 0;i < int(y_st.n_cols);i++){
    
    t_i = arma::find(obs_grid.col(i) == 1);
    m_i = arma::sum(obs_grid.col(i));
    I_m_i = arma::eye(m_i,m_i);
    
    marginal_cov = covariance_mat.rows(t_i);
    marginal_cov = marginal_cov.cols(t_i);
    marginal_cov = marginal_cov + sigma_sq_cur*I_m_i;
    
    y_st_i = y_st.col(i);
    y_st_i = y_st_i.rows(t_i);

    ret = -.5*arma::log_det(marginal_cov).real() +
      -.5*y_st_i.t()*arma::inv(marginal_cov)*y_st_i;
      
  }
  
  return ret(0);
  
}

// [[Rcpp::export]]
double evaluate_mu_prior_cov(const arma::vec mu_cur, const arma::mat covariance_mat) {
  
  arma::vec ret;
  arma::mat eye_m = arma::eye(arma::size(covariance_mat));
  
  ret = -.5*arma::log_det(covariance_mat + eye_m*1e-10).real() +
    -.5*mu_cur.t()*arma::inv(covariance_mat + eye_m*1e-10)*mu_cur;
    
    return ret(0);
    
}


// [[Rcpp::export]]
arma::vec l_mu_MH(arma::mat y_st, const arma::vec grid, const arma::mat obs_grid, double l_param_cur, 
                  const double scale_param_cur, const arma::vec mu_cur, const double sigma_sq_cur,
                  double proposal_l,const bool adapt,const double mcmc_iter,const double prior_a,const double prior_b) {
  
  double l_param_can = std::abs(l_param_cur + arma::randn()*proposal_l);
  
  arma::mat cov_mat_cur = make_cov(grid,l_param_cur,scale_param_cur);
  arma::mat cov_mat_can = make_cov(grid,l_param_can,scale_param_cur);
  
  double l_cur = evaluate_mu_prior_cov(mu_cur, cov_mat_cur);
  double l_can = evaluate_mu_prior_cov(mu_cur, cov_mat_can);
  
  double pi_cur = -(prior_a + 1)*log(l_param_cur) - prior_b/l_param_cur;
  double pi_can = -(prior_a + 1)*log(l_param_can) - prior_b/l_param_can;
  
  double l_acceptance = (l_can + pi_can) - (l_cur + pi_cur);
  
  double accept = 0;
  if(log(arma::randu())<l_acceptance){
    l_param_cur = l_param_can;
    accept += 1;
  }
  
  if(adapt){
    proposal_l = proposal_l*exp(sqrt(1/mcmc_iter)*(std::min(1.0,exp(l_acceptance)) - .234));
  }
  
  arma::vec ret(3);
  ret(0) = l_param_cur;
  ret(1) = accept;
  ret(2) = proposal_l;
  
  return ret;
}


// [[Rcpp::export]]
arma::mat scale_mu_MH(const arma::mat y_st, const arma::vec grid, const arma::mat obs_grid, const double l_param_cur, double scale_param_cur,
                      const arma::vec mu_cur, const double sigma_sq_cur, double proposal_scale,const bool adapt,const double mcmc_iter,const double prior_var) {
  
  
  double scale_param_can = std::abs(scale_param_cur + arma::randn()*proposal_scale);
  
  arma::mat cov_mat_cur = make_cov(grid,l_param_cur,scale_param_cur);
  arma::mat cov_mat_can = make_cov(grid,l_param_cur,scale_param_can);
  
  double l_cur = evaluate_mu_prior_cov(mu_cur, cov_mat_cur);
  double l_can = evaluate_mu_prior_cov(mu_cur, cov_mat_can);
  
  double pi_cur = -pow(scale_param_cur,1)/(2*prior_var);
  double pi_can = -pow(scale_param_can,1)/(2*prior_var);
  
  double l_acceptance = (l_can + pi_can) - (l_cur + pi_cur);
  
  double accept = 0;
  if(log(arma::randu())<l_acceptance){
    scale_param_cur = scale_param_can;
    accept += 1;
  }
  
  if(adapt){
    proposal_scale = proposal_scale*exp(sqrt(1/mcmc_iter)*(std::min(1.0,exp(l_acceptance)) - .234));
  }
  
  arma::vec ret(3);
  ret(0) = scale_param_cur;
  ret(1) = accept;
  ret(2) = proposal_scale;
  
  return ret;
}

// [[Rcpp::export]]
arma::mat sample_lambda_sparse(const arma::mat y_data,const arma::mat W,const arma::vec grid,
                               const arma::vec cov_l_cur,const arma::vec cov_scale_cur,
                               arma::mat lambda_current,const arma::mat eta_current,
                               const double sig_sq_current, const double nu, arma::mat obs_grid) {
  
  arma::mat lambda_temp = lambda_current;
  arma::mat eta_temp = eta_current;
  arma::mat y_star = y_data;
  arma::mat covariance_mat_current;
  arma::mat covariance_mat_inv;
  arma::mat A_inv;
  
  arma::vec lambda_cond_mean;
  arma::mat lambda_cond_var;
  
  const arma::mat id_pen = arma::eye(lambda_current.n_cols - 1,lambda_current.n_cols - 1)*nu;
  const arma::mat I_m = arma::eye(lambda_current.n_rows,lambda_current.n_rows);
  const arma::vec one_col(y_data.n_rows,arma::fill::ones);
  
  for(int k=0;k<int(lambda_current.n_cols);k++){
    lambda_temp.shed_col(k);
    eta_temp.shed_row(k);
    covariance_mat_current = make_cov(grid,cov_l_cur[k],1);
    
    y_star = compute_y_st(y_data,lambda_temp*eta_temp,obs_grid);
    
    covariance_mat_inv = arma::inv(covariance_mat_current + I_m*1e-10)/cov_scale_cur[k];
    
    lambda_cond_var = arma::inv(covariance_mat_inv + 
      W*lambda_temp*lambda_temp.t()*W/nu+ 
      arma::diagmat(arma::sum(one_col*arma::pow(eta_current.row(k),2.0)%obs_grid,1))/sig_sq_current); 
    
    lambda_cond_mean = lambda_cond_var*(arma::sum(one_col*eta_current.row(k)%y_star,1)/sig_sq_current);
    
    lambda_current.col(k) = arma::mvnrnd(lambda_cond_mean,lambda_cond_var); 
    
    lambda_temp = lambda_current;
    eta_temp = eta_current;
  }
  
  return lambda_current;
}


// [[Rcpp::export]]
double evaluate_lambda_cov(const arma::mat covariance_matrix, const arma::mat y_data_minus_k,
                           const arma::mat lambda_minus_k, const arma::vec eta_k, 
                           const double sigma_sq_cur,const double nu_param, const arma::mat obs_grid) {
  
  arma::mat I_k = arma::eye(lambda_minus_k.n_cols,lambda_minus_k.n_cols);
  
  arma::mat prior_var = covariance_matrix - covariance_matrix*lambda_minus_k*
    arma::inv(nu_param*I_k + lambda_minus_k.t()*covariance_matrix*lambda_minus_k)*
    lambda_minus_k.t()*covariance_matrix;
  prior_var = .5*(prior_var + prior_var.t());
  
  arma::vec ret;
  arma::uvec t_i;
  arma::mat marginal_cov;
  double m_i;
  arma::mat I_m_i;
  arma::vec y_st_i;
  
  for(int i = 0;i<int(y_data_minus_k.n_cols);i++){
    
    t_i = arma::find(obs_grid.col(i) == 1);
    m_i = arma::sum(obs_grid.col(i));
    I_m_i = arma::eye(m_i,m_i);
    
    marginal_cov = prior_var.rows(t_i);
    marginal_cov = marginal_cov.cols(t_i);
    marginal_cov = pow(eta_k(i),2.0)*marginal_cov + sigma_sq_cur*I_m_i;
    
    y_st_i = y_data_minus_k.col(i);
    y_st_i = y_st_i.rows(t_i);
    
    ret = -.5*arma::log_det(marginal_cov).real() +
      -.5*y_st_i.t()*arma::inv(marginal_cov)*y_st_i;
      
  }
  
  return ret(0);
}

// [[Rcpp::export]]
double evaluate_lambda_prior_cov(const arma::vec lambda_k, const arma::mat lambda_minus_k, 
                                 const arma::mat covariance_mat,const double nu_param) {
  
  arma::vec ret;
  arma::mat eye_m = arma::eye(arma::size(covariance_mat));
  
  arma::mat I_k = arma::eye(lambda_minus_k.n_cols,lambda_minus_k.n_cols);
  
  arma::mat prior_var = covariance_mat - covariance_mat*lambda_minus_k*
    arma::inv(nu_param*I_k + lambda_minus_k.t()*covariance_mat*lambda_minus_k)*
    lambda_minus_k.t()*covariance_mat;
  prior_var = .5*(prior_var + prior_var.t());
  
  ret = -.5*arma::log_det(prior_var + eye_m*1e-10).real() +
    -.5*lambda_k.t()*arma::inv(prior_var + eye_m*1e-10)*lambda_k;
    
    return ret(0);
    
}

// [[Rcpp::export]]
arma::mat l_param_MH(arma::vec l_param_cur,arma::vec proposal_l, const arma::vec cov_scale_current,const arma::vec grid, const arma::mat obs_grid,
                     const arma::mat W,const arma::mat y_data,const arma::mat lambda_current, const arma::mat eta_current, 
                     const double sigma_sq_current,const double nu_param,
                     const bool adapt, const double mcmc_iter,const double prior_a,const double prior_b) {
  
  arma::mat y_minus_k;
  arma::mat lambda_minus_k;
  arma::mat eta_minus_k;
  arma::vec eta_k;
  
  double l_cur = 0;
  double l_can = 0;
  double pi_cur = 0;
  double pi_can = 0;
  double proposal_cur_given_can = 0;
  double proposal_can_given_cur = 0;
  
  double l_acceptance = 1;
  arma::vec accept(size(l_param_cur),arma::fill::zeros);
  
  arma::vec l_param_can = arma::abs(l_param_cur + arma::randn(arma::size(l_param_cur)) % arma::sqrt(proposal_l));
  
  for(int k = 0;k<int(l_param_cur.n_rows);k++){
    
    arma::mat cov_mat_cur = make_cov(grid,l_param_cur[k],cov_scale_current[k]);
    arma::mat cov_mat_can = make_cov(grid,l_param_can[k],cov_scale_current[k]);
    
    lambda_minus_k = lambda_current;
    lambda_minus_k.shed_col(k);
    eta_k = eta_current.row(k).t();
    eta_minus_k = eta_current;
    eta_minus_k.shed_row(k);
    y_minus_k = compute_y_st(y_data,lambda_minus_k*eta_minus_k,obs_grid);
    
    // l_cur = evaluate_lambda_cov(cov_mat_cur,y_minus_k,W*lambda_minus_k,eta_k,sigma_sq_current,nu_param,obs_grid);
    // l_can = evaluate_lambda_cov(cov_mat_can,y_minus_k,W*lambda_minus_k,eta_k,sigma_sq_current,nu_param,obs_grid);
    
    l_cur = evaluate_lambda_prior_cov(lambda_current.col(k), lambda_minus_k,cov_mat_cur,nu_param);
    l_can = evaluate_lambda_prior_cov(lambda_current.col(k), lambda_minus_k,cov_mat_can,nu_param);
    
    pi_cur = -(prior_a + 1)*log(l_param_cur[k]) - prior_b/l_param_cur[k];
    pi_can = -(prior_a + 1)*log(l_param_can[k]) - prior_b/l_param_can[k];
    
    l_acceptance = (l_can + pi_can + proposal_cur_given_can) - (l_cur + pi_cur + proposal_can_given_cur);
    
    if(log(arma::randu())<l_acceptance){
      l_param_cur[k] = l_param_can[k];
      accept[k] += 1;
    }
    
    if(adapt){
      proposal_l[k] = proposal_l[k]*exp(sqrt(1/mcmc_iter)*(std::min(1.0,exp(l_acceptance)) - .234));
    }
    
    
  }
  arma::mat ret(l_param_cur.n_rows,3);
  ret.col(0) = l_param_cur;
  ret.col(1) = accept;
  ret.col(2) = proposal_l;
  return ret;
}


// [[Rcpp::export]]
arma::mat scale_param_MH(arma::vec cov_scale_cur, arma::vec proposal_scale,const arma::vec l_param_cur,const arma::vec grid,
                         const arma::mat obs_grid, const arma::mat W,const arma::mat y_data,const arma::mat lambda_current, 
                         const arma::mat eta_current, const double sigma_sq_current,const double nu_param,
                         const bool adapt, const double mcmc_iter,const double prior_var) {
  
  arma::mat y_minus_k;
  arma::mat lambda_minus_k;
  arma::mat eta_minus_k;
  arma::vec eta_k;
  
  double l_cur = 0;
  double l_can = 0;
  double pi_cur = 0;
  double pi_can = 0;
  double proposal_cur_given_can = 0;
  double proposal_can_given_cur = 0;
  
  double l_acceptance = 1;
  arma::vec accept(size(cov_scale_cur),arma::fill::zeros);
  
  arma::vec cov_scale_can = arma::abs(cov_scale_cur + arma::randn(arma::size(cov_scale_cur)) % arma::sqrt(proposal_scale));
  
  for(int k = 0;k<int(l_param_cur.n_rows);k++){
    
    arma::mat cov_mat_cur = make_cov(grid,l_param_cur[k],cov_scale_cur[k]);
    arma::mat cov_mat_can = make_cov(grid,l_param_cur[k],cov_scale_can[k]);
    
    lambda_minus_k = lambda_current;
    lambda_minus_k.shed_col(k);
    eta_k = eta_current.row(k).t();
    eta_minus_k = eta_current;
    eta_minus_k.shed_row(k);
    y_minus_k = compute_y_st(y_data,lambda_minus_k*eta_minus_k,obs_grid);
    
    // l_cur = evaluate_lambda_cov(cov_mat_cur,y_minus_k,W*lambda_minus_k,eta_k,sigma_sq_current,nu_param,obs_grid);
    // l_can = evaluate_lambda_cov(cov_mat_can,y_minus_k,W*lambda_minus_k,eta_k,sigma_sq_current,nu_param,obs_grid);
    
    l_cur = evaluate_lambda_prior_cov(lambda_current.col(k), lambda_minus_k,cov_mat_cur,nu_param);
    l_can = evaluate_lambda_prior_cov(lambda_current.col(k), lambda_minus_k,cov_mat_can,nu_param);
    
    pi_cur = -pow(cov_scale_cur[k],1)/(2*prior_var);
    pi_can = -pow(cov_scale_can[k],1)/(2*prior_var);
    
    l_acceptance = (l_can + pi_can + proposal_cur_given_can) - (l_cur + pi_cur + proposal_can_given_cur);
    
    if(log(arma::randu())<l_acceptance){
      cov_scale_cur[k] = cov_scale_can[k];
      accept[k] += 1;
    }
    
    if(adapt){
      proposal_scale[k] = proposal_scale[k]*exp(sqrt(1/mcmc_iter)*(std::min(1.0,exp(l_acceptance)) - .234));
    }
    
    
  }
  arma::mat ret(cov_scale_cur.n_rows,3);
  ret.col(0) = cov_scale_cur;
  ret.col(1) = accept;
  ret.col(2) = proposal_scale;
  return ret;
}




// [[Rcpp::export]]
arma::mat sample_eta_sparse(const arma::mat y_data, arma::mat eta_current, arma::vec psi_current, 
                            const arma::mat lambda_current, const double sig_sq_current, 
                            const double nu_param,const arma::mat obs_grid) {
  double n = y_data.n_cols;
  arma::mat y_st = y_data;
  arma::mat lambda_minus_k = lambda_current;
  arma::mat eta_minus_k = eta_current;
  arma::vec one_vec_n(n,arma::fill::ones);
  arma::mat I_n = arma::eye(n,n);
  arma::vec D = one_vec_n;
  arma::vec D_inv = D;
  arma::mat J_n = one_vec_n*one_vec_n.t();
  arma::vec lambda_k;
  
  arma::vec temp;
  arma::mat eta_cond_var; 
  arma::vec eta_cond_mean;
  
  for(int k = 0;k<int(eta_current.n_rows);k++){
    lambda_minus_k.shed_col(k);
    eta_minus_k.shed_row(k);
    lambda_k = lambda_current.col(k);
    
    y_st = compute_y_st(y_data,lambda_minus_k*eta_minus_k,obs_grid);
    
    D = arma::diagvec((arma::sqrt(obs_grid)%arma::repmat(lambda_k,1,n)).t()*
      (arma::sqrt(obs_grid)%arma::repmat(lambda_k,1,n)))/sig_sq_current + one_vec_n/(one_vec_n*psi_current.row(k));
    D_inv  = 1/D;
    
    temp = psi_current.row(k)*nu_param + arma::sum(D_inv);
    eta_cond_var = arma::diagmat(D_inv) - D_inv*D_inv.t()/temp(0);
    
    eta_cond_mean = eta_cond_var*y_st.t()*lambda_k/sig_sq_current;
    
    temp = arma::mvnrnd(eta_cond_mean,eta_cond_var);
    eta_current.row(k) = temp.t();
    
    lambda_minus_k = lambda_current;
    eta_minus_k = eta_current;
  }
  return eta_current;
}

// [[Rcpp::export]]
arma::vec sample_psi_sparse(arma::mat eta_current, const double nu_param,const double px_a, const double px_b) {
  double n = eta_current.n_cols;
  int K = eta_current.n_rows;
  
  arma::vec one_vec_n(n,arma::fill::ones);
  arma::mat I_n = arma::eye(n,n);
  arma::mat J_n = one_vec_n*one_vec_n.t();
  
  arma::vec temp; 
  arma::vec psi_ret = arma::ones(K);
  
  double alpha_cond = px_a + n/2;
  double beta_cond = 1;
  
  for(int k = 0;k<K;k++){
    
    temp = arma::sum(arma::pow(eta_current.row(k),2.0),1) +  pow(arma::sum(eta_current.row(k)),2.0)/nu_param;
    beta_cond = px_b + temp(0)/2;
    
    psi_ret.row(k) = 1/arma::randg(arma::distr_param(alpha_cond,1/beta_cond));
  }
  return psi_ret;
}

// [[Rcpp::export]]
double sample_sig_sq_sparse(const arma::mat y_data,const arma::mat lambda_current,
                            const arma::mat eta_current,const arma::mat obs_grid) {
  double alpha_cond = .000001 + .5*arma::accu(obs_grid);
  arma::mat y_st = compute_y_st(y_data,lambda_current*eta_current,obs_grid);
  double beta_cond = .000001 + .5*arma::trace(arma::trans(y_st)*(y_st));
  return 1/arma::randg(arma::distr_param(alpha_cond,1/beta_cond));
}


// [[Rcpp::export]]
arma::mat cholupdate(arma::mat RMat, arma::vec xVec, double nu){
  // Adapted from code written by Piotr Orlowski
  // See https://github.com/piotrek-orlowski/ukfRcpp/blob/master/src/cholupdate.cpp) for documentation 
  
  // Take square root of nu
  nu = pow(nu,0.5);
  
  // Take size of problem
  int Nx = xVec.n_elem;
  
  // scale by nu
  xVec *= nu;
  
  double r, c, s;
  r=0;
  
  // algo
  for(int kk = 0; kk < Nx; kk++){
    
    r = sqrt(pow(RMat(kk,kk),2.0) - pow(xVec(kk),2.0));
    
    c = r / RMat(kk,kk);
    
    s = xVec(kk) / RMat(kk,kk);
    
    RMat(kk,kk) = r;
    if(kk < Nx-1){
      RMat.submat(kk+1,kk,Nx-1,kk) = (RMat.submat(kk+1,kk,Nx-1,kk) - s * xVec.subvec(kk+1,Nx-1)) / c;
      xVec.subvec(kk+1,Nx-1) = c * xVec.subvec(kk+1,Nx-1) - s * RMat.submat(kk+1,kk,Nx-1,kk);
    }
  }
  
  return RMat;
}


// [[Rcpp::export]]
arma::mat sample_eta_cholupdate(const arma::mat y_data, arma::mat eta_current, arma::vec psi_current, 
                                const arma::mat lambda_current, const double sig_sq_current, 
                                const double nu_param,const arma::mat obs_grid) {
  double n = y_data.n_cols;
  arma::mat y_st = y_data;
  arma::mat lambda_minus_k = lambda_current;
  arma::mat eta_minus_k = eta_current;
  arma::vec one_vec_n(n,arma::fill::ones);
  arma::mat I_n = arma::eye(n,n);
  arma::vec D = one_vec_n;
  arma::vec D_inv = D;
  arma::mat chol_decomp;
  arma::mat J_n = one_vec_n*one_vec_n.t();
  arma::vec lambda_k;
  
  arma::vec temp;
  arma::mat eta_cond_var; 
  arma::vec eta_cond_mean;
  
  for(int k = 0;k<int(eta_current.n_rows);k++){
    lambda_minus_k.shed_col(k);
    eta_minus_k.shed_row(k);
    lambda_k = lambda_current.col(k);
    
    y_st = compute_y_st(y_data,lambda_minus_k*eta_minus_k,obs_grid);
    
    D = arma::diagvec((arma::sqrt(obs_grid)%arma::repmat(lambda_k,1,n)).t()*
      (arma::sqrt(obs_grid)%arma::repmat(lambda_k,1,n)))/sig_sq_current + one_vec_n/(one_vec_n*psi_current.row(k));
    D_inv  = 1/D;
    
    temp = psi_current.row(k)*nu_param + arma::sum(D_inv);
    eta_cond_var = arma::diagmat(D_inv) - D_inv*D_inv.t()/temp(0);
    
    eta_cond_mean = eta_cond_var*y_st.t()*lambda_k/sig_sq_current;
    
    chol_decomp = cholupdate(arma::diagmat(arma::sqrt(D_inv)),D_inv,1/temp(0));
    
    temp = eta_cond_mean + chol_decomp*arma::randn(n);
    eta_current.row(k) = temp.t();
    
    lambda_minus_k = lambda_current;
    eta_minus_k = eta_current;
  }
  return eta_current;
}


// [[Rcpp::export]]
arma::mat sample_b(const arma::mat y_data,const arma::mat x_data,
                   arma::mat b_current,const arma::mat lambda_current, const arma::vec psi_current,
                   const double sig_sq_current,const arma::mat obs_grid) {
  arma::mat temp;
  arma::mat x_temp = x_data;
  arma::mat b_temp = b_current;
  arma::mat y_star = y_data;
  
  arma::vec b_cond_mean;
  arma::mat b_cond_var;
  
  const arma::mat prior_var_inv = arma::diagmat(1/psi_current);
  const arma::vec one_col(y_data.n_rows,arma::fill::ones);
  
  
  for(int q=0;q<int(b_current.n_cols);q++){
    x_temp.shed_row(q);
    b_temp.shed_col(q);
    
    y_star = compute_y_st(y_data,lambda_current*b_temp*x_temp,obs_grid);
    
    temp = arma::diagmat(arma::sum(one_col*arma::pow(x_data.row(q),2.0)%obs_grid,1))/sig_sq_current;
    b_cond_var = arma::inv(prior_var_inv + lambda_current.t()*temp*lambda_current);
    
    b_cond_mean = b_cond_var*lambda_current.t()*(arma::sum(one_col*x_data.row(q)%y_star,1)/sig_sq_current);
    
    b_current.col(q) = arma::mvnrnd(b_cond_mean,b_cond_var); 
    
    b_temp = b_current;
    x_temp = x_data;
  }
  
  return b_current;
}

// [[Rcpp::export]]
arma::vec sample_psi_regression(const arma::mat eta_current, const arma::mat b_current,const double nu_param,const double px_a, const double px_b) {
  double n = eta_current.n_cols;
  int K = eta_current.n_rows;
  int Q = b_current.n_cols;
  
  arma::vec one_vec_n(n,arma::fill::ones);
  arma::mat I_n = arma::eye(n,n);
  arma::mat J_n = one_vec_n*one_vec_n.t();
  
  arma::vec temp; 
  arma::vec psi_ret = arma::ones(K);
  
  double alpha_cond = px_a + (n + Q)/2;
  double beta_cond = 1;
  
  for(int k = 0;k<K;k++){
    
    temp = arma::sum(arma::pow(eta_current.row(k),2.0),1) +  
      pow(arma::sum(eta_current.row(k)),2.0)/nu_param + 
      arma::sum(arma::pow(b_current.row(k),2.0),1);
    beta_cond = px_b + temp(0)/2;
    
    psi_ret.row(k) = 1/arma::randg(arma::distr_param(alpha_cond,1/beta_cond));
  }
  return psi_ret;
}

// [[Rcpp::export]]
arma::mat sample_z_latent(const arma::mat z_data, const arma::mat f,const arma::mat obs_grid) {
  
  int n = z_data.n_cols;
  int M = z_data.n_rows;
  double trunc_bound = 1e5;
  
  arma::mat z_ret(M,n,arma::fill::zeros);
  
  for(int i = 0; i<n;i++){
    for(int m = 0; m<M; m++){
      if(obs_grid(m,i) == 1){ // sample latent z for observed values
        if(z_data(m,i) == 1){
          z_ret(m,i) = r_truncnorm(f(m,i),1,0,trunc_bound);
        }else{
          z_ret(m,i) = r_truncnorm(f(m,i),1,-trunc_bound,0);
        }
      }
    }
  }
  return z_ret;
}

// [[Rcpp::export]]
arma::mat sample_z_complete(const arma::mat z_data, const arma::mat f,const arma::mat obs_grid) {
  
  int n = z_data.n_cols;
  int M = z_data.n_rows;
  
  arma::mat z_ret(M,n,arma::fill::zeros);
  
  for(int i = 0; i<n;i++){
    for(int m = 0; m<M; m++){
      if(obs_grid(m,i) == 1){ // sample only for unobserved values of z
        z_ret(m,i) = z_data(m,i);
      }else{
        z_ret(m,i) = arma::randu() < arma::normcdf(f(m,i)); 
      }
    }
  }
  return z_ret;
}



// [[Rcpp::export]]
arma::mat sample_b_marginal(const arma::mat y_data,const arma::mat x_data,
                            arma::mat b_current,const arma::mat lambda_current, const arma::vec psi_current,
                            const double sig_sq_current,const arma::mat obs_grid, const double nu_param) {
  int n = y_data.n_cols;
  int K = lambda_current.n_cols;
  arma::mat temp;
  arma::mat x_temp = x_data;
  arma::mat b_temp = b_current;
  arma::mat y_star = y_data;
  arma::mat OO_i;
  arma::mat lambda_prod;
  arma::rowvec y_lambda_prod(K);
  arma::rowvec term21(K);
  arma::rowvec term22(K);
  arma::rowvec term2(K);  
  arma::rowvec y_prod_inv(K);
  arma::mat lambda_prod_inv(K,K,arma::fill::zeros);
  arma::mat term11(K,K,arma::fill::zeros);
  arma::mat term12(K,K,arma::fill::zeros);
  arma::mat term13(K,K,arma::fill::zeros);
  arma::mat term1;
  arma::mat zero_K(K,K,arma::fill::zeros);
  arma::mat eye_K = arma::eye(K,K);
  arma::rowvec zero_row_K(K,arma::fill::zeros);
  arma::mat temp_inv;
  
  arma::vec b_cond_mean;
  arma::mat b_cond_var;
  arma::rowvec x_q;
  double x_q_i;
  
  const arma::mat prior_var_inv = arma::diagmat(1/psi_current);
  const arma::vec one_col(y_data.n_rows,arma::fill::ones);
  
  for(int q=0;q<int(b_current.n_cols);q++){
    x_temp.shed_row(q);
    b_temp.shed_col(q);
    x_q = x_data.row(q);
    y_star = compute_y_st(y_data,lambda_current*b_temp*x_temp,obs_grid);
    
    term11 = zero_K;
    term12 = zero_K;
    term13 = zero_K;
    term21 = zero_row_K;
    term22 = zero_row_K;
    for(int i = 0; i<n;i++){
      x_q_i = x_q(i);
      OO_i = arma::diagmat(obs_grid.col(i));
      lambda_prod = lambda_current.t()*OO_i*lambda_current;
      y_lambda_prod = y_star.col(i).t()*lambda_current;
      temp_inv = eye_K - arma::inv(sig_sq_current*eye_K + lambda_prod)*lambda_prod;
      lambda_prod_inv = lambda_prod/sig_sq_current*temp_inv;
      y_prod_inv = y_lambda_prod/sig_sq_current*temp_inv;
      term13 += lambda_prod_inv;
      term12 += x_q_i*lambda_prod_inv;
      term11 += pow(x_q_i,2.0)*lambda_prod_inv;
      term22 += y_prod_inv;
      term21 += x_q_i*y_prod_inv;  
    }
    
    temp_inv = arma::inv(-(nu_param + n)*eye_K + term13)*term12;
    term1 = term11 - term12*temp_inv;
    term2 = term21 - term22*temp_inv;
    
    b_cond_var = arma::inv(term1 + prior_var_inv);
    b_cond_mean = b_cond_var*term2.t();
    
    b_current.col(q) = arma::mvnrnd(b_cond_mean,b_cond_var); 
    
    b_temp = b_current;
    x_temp = x_data;
  }
  
  return b_current;
}

// [[Rcpp::export]]
arma::cube compute_fun_prod(const arma::mat z_data,const arma::mat basis_funs) {
  int n = z_data.n_cols;
  int P = basis_funs.n_cols;
  int M = z_data.n_rows;
  
  arma::cube ret_prod(M,P,n);
  
  for(int i = 0; i < n; i++){
    for(int p = 0;p<P; p++){
      ret_prod.slice(i).col(p) = z_data.col(i)%basis_funs.col(p);
    }
  }
  
  return ret_prod;
}



// [[Rcpp::export]]
arma::cube compute_historical_integrals(const arma::mat z_data,const arma::cube basis_funs,const arma::mat W) {
  int n = z_data.n_cols;
  int P = basis_funs.n_slices;
  int M = basis_funs.n_rows;
  
  arma::cube ret_integrals(M,P,n);
  
  for(int i = 0; i < n; i++){
    for(int p = 0; p < P; p++){
      ret_integrals.slice(i).col(p) = basis_funs.slice(p)*W*z_data.col(i);
    }
  }
  
  return ret_integrals;
}

// [[Rcpp::export]]
arma::mat compute_fun_reg_fit(const arma::cube historical_integrals,const arma::vec beta_current) {
  int n = historical_integrals.n_slices;
  int M = historical_integrals.n_rows;
  
  arma::mat ret_fits(M,n);
  
  for(int i = 0; i < n; i++){
    ret_fits.col(i) = historical_integrals.slice(i)*beta_current;
  }
  
  return ret_fits;
}

// [[Rcpp::export]]
arma::vec sample_beta_historical(const arma::mat y_data, const arma::cube historical_integrals,
                      const double sig_sq_current,const arma::mat obs_grid,const double prior_var) {
  int n = historical_integrals.n_slices;
  int P = historical_integrals.n_cols;
  
  arma::mat sum_cov(P,P,arma::fill::zeros);
  arma::vec sum_mean(P,arma::fill::zeros);

  for(int i = 0; i < n; i++){
    sum_cov += historical_integrals.slice(i).t()*diagmat(obs_grid.col(i))*historical_integrals.slice(i);
    sum_mean += historical_integrals.slice(i).t()*y_data.col(i);
  }
  
  arma::mat beta_cond_var = arma::inv(sum_cov/sig_sq_current + arma::eye(P,P)/prior_var);
  arma::vec beta_cond_mean = beta_cond_var*sum_mean/sig_sq_current;
  
  return arma::mvnrnd(beta_cond_mean,beta_cond_var);
}

// [[Rcpp::export]]
arma::mat compute_concurrent_fit(const arma::cube z_data, const arma::mat basis_funs,const arma::mat beta_current) {
  
  int Q = z_data.n_slices;
  int M = z_data.n_rows;
  int n = z_data.n_cols;
  
  arma::cube temp_fun_product;
  arma::mat ret_fit(M,n,arma::fill::zeros);
  
  for(int q = 0; q  < Q ; q++){
    temp_fun_product = compute_fun_prod(z_data.slice(q),basis_funs);
    ret_fit += compute_fun_reg_fit(temp_fun_product,beta_current.col(q));
  }
  
  return ret_fit;
  
}

// [[Rcpp::export]]
arma::mat sample_beta_concurrent(const arma::mat y_data, const arma::cube z_data, const arma::mat basis_funs,
                                 arma::mat beta_current, const double sig_sq_current, const arma::mat obs_grid,const double prior_var) {
  
  int Q = z_data.n_slices;
  arma::field<arma::cube> temp_fun_product(Q);
  
  for(int q = 0; q  < Q ; q++){
    temp_fun_product(q) = compute_fun_prod(z_data.slice(q),basis_funs);
  }
  
  int M = y_data.n_rows;
  int n = y_data.n_cols;
  int P = basis_funs.n_cols;
  
  arma::mat sum_cov(P,P);
  arma::vec sum_mean(P);
  arma::mat beta_cond_var(P,P);
  arma::vec beta_cond_mean(P);
  
  arma::mat temp_fit(M,n);
  arma::mat y_star = y_data;
  
  for(int q = 0; q < Q; q++){
    temp_fit.zeros();
    for(int q_star = 0; q_star < Q; q_star++){
     if(q_star != q){
       temp_fit += compute_fun_reg_fit(temp_fun_product(q_star),beta_current.col(q_star));
     }
    }
    y_star = compute_y_st(y_data,temp_fit,obs_grid);
  
    sum_cov.zeros();
    sum_mean.zeros();
    
    for(int i = 0; i < n; i++){
      sum_cov += temp_fun_product(q).slice(i).t()*diagmat(obs_grid.col(i))*temp_fun_product(q).slice(i);
      sum_mean += temp_fun_product(q).slice(i).t()*y_star.col(i);
    }
    
    beta_cond_var = arma::inv(sum_cov/sig_sq_current + arma::eye(P,P)/prior_var);
    beta_cond_mean = beta_cond_var*sum_mean/sig_sq_current;
    
    beta_current.col(q) = arma::mvnrnd(beta_cond_mean,beta_cond_var);
    
    y_star = y_data; 
  }
  
  return beta_current;
}


// [[Rcpp::export]]
arma::vec regress_basis(const arma::cube basis_ints,const arma::mat y_data, const arma::mat obs_grid) {
  int n = y_data.n_cols;
  int P = basis_ints.n_cols;
  
  arma::mat temp_mat(P,P,arma::fill::zeros);
  arma::vec temp_vec(P,arma::fill::zeros);

  for(int i = 0; i < n; i++){
    temp_mat += basis_ints.slice(i).t()*arma::diagmat(obs_grid.col(i))*basis_ints.slice(i);
    temp_vec += basis_ints.slice(i).t()*y_data.col(i);
  }
  
  return arma::inv(temp_mat)*temp_vec;
}
