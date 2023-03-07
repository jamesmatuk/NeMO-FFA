library(tidyverse)
library(lubridate)
library(reshape2)
library(haven)
library(Rcpp)
library(ggplot2)
library(splines)
library(rootSolve)

##### Import data #####

set.seed(1234)

sourceCpp("./src/nemo_ffa.cpp")
sourceCpp("./src/msf.cpp")
source("./src/nemo_ffa.R")

cebu_long <- read_dta("data/mlong.dta")
cebu_mort <- read_dta("data/mmort.dta")
cebu_birth <- read_dta("data/mbirth2.dta")

num_subjects <-dim(unique( cebu_long[,c("basebrgy","basewman")]))[1]

##### import relevant data #####
# response - birthweight and weight every two months until 2 yrs old 

visit_num <- NULL
for(i in 1:12){
  visit_num <- c(visit_num,rep(i,num_subjects))
}
cebu_long$visit_num = visit_num

cebu_weight <- pivot_wider(cebu_long,id_cols = c(basebrgy,basewman),names_from = visit_num,values_from = c(weight))
cebu_birthweight <- cebu_birth$WEIGHT1

# scalar covariates

cebu_scalar_covariates <-  cebu_birth %>% 
  dplyr::select(basebrgy,basewman,height,sexchild,stratum,mnthbrth)

# logitudinal covariates

# breastfeading status - at birth and every two months afterwards

cebu_bresys <- pivot_wider(cebu_long,id_cols = c(basebrgy,basewman),names_from = visit_num,values_from = c(bresys))
cebu_brest7 <- pivot_wider(cebu_long,id_cols = c(basebrgy,basewman),names_from = visit_num,values_from = c(BREST7))

cebu_bf <- array(NA,dim = c(num_subjects,13))
cebu_bf[,1] <- pmax(cebu_birth$BRESMLK1,cebu_birth$BRETMKL2,na.rm = T) # 0 = no, 1 = yes, 2 = yes from wet nurse
cebu_bf[cebu_bf[,1] >= 1,1] <- 1


for(vis_num in 1:12){
  cebu_bf[,vis_num+1] <-  pmax(as.matrix(cebu_bresys[,vis_num+2]),as.matrix(cebu_brest7[,vis_num+2]),na.rm = T) 
}

# illnessnes recorded 

cebu_diar <- pivot_wider(cebu_long,id_cols = c(basebrgy,basewman),names_from = visit_num,values_from = c(DIA7DY))
cebu_fv <- pivot_wider(cebu_long,id_cols = c(basebrgy,basewman),names_from = visit_num,values_from = c(FEVER7))
cebu_cg <- pivot_wider(cebu_long,id_cols = c(basebrgy,basewman),names_from = visit_num,values_from = c(COUGH7))

# subjects to keep for analysis

#   conditions: 1) child alive the entire time during study
#               2) child born alive
#               3) was there at least 1 recording of their weight
#               4) mother's height is recorded

non_missing_ind <- apply(is.na(cebu_weight[,3:14]),1,sum) + is.na(cebu_birthweight) < 13

height_recorded <- !is.na(cebu_scalar_covariates$height)

alive_ind <- is.na(cebu_mort$agedied) & (cebu_birth$statusbb == 1)

keep_ind <- which(non_missing_ind & height_recorded & alive_ind)

cebu_weight <- cebu_weight[keep_ind,]
cebu_birthweight <- cebu_birthweight[keep_ind]

cebu_scalar_covariates <- cebu_scalar_covariates[keep_ind,]

cebu_bf <- cebu_bf[keep_ind,]

cebu_diar <- cebu_diar[keep_ind,]
cebu_fv <- cebu_fv[keep_ind,]
cebu_cg <- cebu_cg[keep_ind,]


##### format data for analysis #####

# response 
y_common <- rbind(cebu_birthweight,t(as.matrix(cebu_weight[,3:14])))
y_common[is.na(y_common)] = 0
obs_grid_y <- apply(y_common>0,2,as.numeric)

y_scale <- sqrt(norm(y_common))
y_common <- y_common/y_scale

# scalar covariates

x <- as.matrix(dplyr::select(cebu_scalar_covariates,height,sexchild,stratum,mnthbrth))
# transform covariates 
ht_mn <- mean(x[,1])
ht_sd <- sd(x[,1])
x[,1] <- (x[,1] - mean(x[,1]))/sd(x[,1]) # standardize heights of mothers
x[,2] <- (x[,2] - 1) # sex: 0 = male, 1 = female
x[,3] <- (x[,3] - 1) # strata: 0 = urban, 1 = rural
x[,4] <- as.numeric(x[,4] <= 11 & x[,4] >= 6)  # season of birth: 0 = rainy, 1 = dry

x <- t(x)

# breastmilk covariates

z_bf <- t(cebu_bf)
obs_grid_bf <- as.matrix(!is.na(z_bf))

# illness covariates

num_covs <- 3
z_ill <- array(0,dim = c(dim(y_common) - c(1,0),num_covs))

z_ill[,,1] <- t(as.matrix(cebu_diar[,3:14]))
z_ill[,,2] <- t(as.matrix(cebu_fv[,3:14]))
z_ill[,,3] <- t(as.matrix(cebu_cg[,3:14]))

obs_grid_ill <- array(0,dim = c(dim(y_common) - c(1,0),num_covs))
for(cov_ind in 1:num_covs){
  obs_grid_ill[,,cov_ind] <- as.matrix(!is.na(z_ill[,,cov_ind]))
}

# time scale 
age_grid <- 0:12/6
M <- length(age_grid)
w <- diff(c(age_grid[1],(age_grid[2:M]+age_grid[1:(M-1)])/2,age_grid[M]))
W <- diag(w)

##### Initialization for MCMC #####

N <- dim(y_common)[2]
K <- 10
K_bf <- 10
Q <- dim(x)[1]
one_vec_N <- rep(1,N)
nu_eta <- 1e-4
nu_param <- 1e-4

PX <- TRUE
px_b <- 10
px_a <- px_b + 1

inv_g_hyper <- length_scale_hyper(age_grid)
prior_a <- inv_g_hyper[1]
prior_b <- inv_g_hyper[2]

prior_var <- 1

# initilize ortho gp 
lambda_y_cur <- matrix(0,nrow = M,ncol = K)
xi_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
b_cur <- matrix(0,nrow = K,ncol = Q)
eta_y_cur <- b_cur%*%x + xi_cur
psi_y_cur <- rep(1,K)
sig_sq_cur <- var((y_common-lambda_y_cur%*%eta_y_cur)[obs_grid_y>0])
l_mu_y_current <- 1
scale_mu_y_current <- 1
l_y_param_cur <- rep(1,K)
scale_y_param_cur <- rep(1,K)

# y-specify proposal kernel parameters
proposal_l_mu_y <- .0001
proposal_scale_mu_y <- .0001
proposal_l_y <- rep(.001,K)
proposal_scale_y <- rep(.001,K)
accept_l_mu_y_count <- 0
accept_scale_mu_y_count <- 0
accept_l_y_count <- rep(0,K)
accept_scale_y_count <- rep(0,K)

# initialize historical regression
num_centers<- 4
source("create_tent_basis.R")

hist_basis <- basis_funs
P <- dim(hist_basis)[3]
beta_bf_cur <- rnorm(P,0,1)

mu_bf_cur <- rep(0,M)
lambda_bf_cur <- matrix(0,nrow = M,ncol = K_bf)
eta_bf_cur <- matrix(rnorm(K*N,0,1),nrow = K_bf,ncol = N)
psi_bf_cur <- rep(1,K_bf)
l_mu_bf_current <- 1
scale_mu_bf_current <- 1
l_bf_param_cur <- rep(1,K_bf)
scale_bf_param_cur <- rep(1,K_bf)

z_bf_latent_cur <- array(0,dim = dim(z_bf))
z_bf_complete <- array(0,dim = dim(z_bf))

historical_integrals <- compute_historical_integrals(z_bf_complete,hist_basis,W)
historical_fit <- rbind(0*one_vec_N,compute_fun_reg_fit(historical_integrals,beta_bf_cur))

# z_bf proposal kernel parameters
proposal_l_mu_bf <- .0001
proposal_scale_mu_bf <- .0001
proposal_l_bf <- rep(.001,K_bf)
proposal_scale_bf <- rep(.001,K_bf)
accept_l_mu_bf_count <- 0
accept_scale_mu_bf_count <- 0
accept_l_bf_count <- rep(0,K_bf)
accept_scale_bf_count <- rep(0,K_bf)

# initialize concurrent regression

mu_ill_cur <- array(0,dim = c(M-1,num_covs))
l_mu_ill_current <- rep(1,num_covs)
scale_mu_ill_current <- rep(1,num_covs)
z_ill_latent_cur <- array(0,dim = dim(z_ill))
z_ill_complete <- array(0,dim = dim(z_ill))

proposal_l_mu_ill <- rep(.0001,num_covs)
proposal_scale_mu_ill <- rep(.0001,num_covs)
accept_l_mu_ill_count <- rep(0,num_covs)
accept_scale_mu_ill_count <- rep(0,num_covs)

beta_ill_basis_size <- 4
beta_ill_basis <- bs(age_grid[2:M],df = beta_ill_basis_size,intercept = T)
beta_ill_basis <- matrix(beta_ill_basis,nrow = M-1)
beta_ill_coefs_cur <- array(0,dim = c(beta_ill_basis_size,num_covs))

concurrent_fit <- rbind(0*one_vec_N,compute_concurrent_fit(z_ill_complete,beta_ill_basis,beta_ill_coefs_cur))

# pre-allocate memory for MCMC samples
n_iter <- 1000
lag <- 1
n_save <- n_iter/lag
mu_y_save <- matrix(0,nrow = M,ncol = n_save)
scale_mu_y_save <- rep(0,n_save)
l_mu_y_save <- rep(0,n_save)
lambda_y_save <- array(0,dim = c(M,K,n_save))
xi_save <- array(0,dim = c(K,N,n_save))
b_save <- array(0,dim = c(K,Q,n_save))
eta_y_save <- array(0,dim = c(K,N,n_save))
psi_y_save <- array(0,dim = c(K,n_save))
l_y_param_save <- array(0,dim = c(K,n_save))
scale_y_param_save <- array(0,dim = c(K,n_save))

mu_bf_save <- matrix(0,nrow = M,ncol = n_save)
scale_mu_bf_save <- rep(0,n_save)
l_mu_bf_save <- rep(0,n_save)
lambda_bf_save <- array(0,dim = c(M,K_bf,n_save))
eta_bf_save <- array(0,dim = c(K_bf,N,n_save))
psi_bf_save <- array(0,dim = c(K_bf,n_save))
l_bf_param_save <- array(0,dim = c(K_bf,n_save))
scale_bf_param_save <- array(0,dim = c(K_bf,n_save))
beta_bf_save <- array(0,dim = c(P,n_save))

mu_ill_save <- array(0,dim = c(M-1,num_covs,n_save))
l_mu_ill_save <- array(0,dim = c(num_covs,n_save))
scale_mu_ill_save <- array(0,dim = c(num_covs,n_save))
beta_ill_coefs_save <- array(0,dim = c(beta_ill_basis_size,num_covs,n_save))

sig_sq_save <- rep(0,n_save)

##### run MCMC #####
ptm <- proc.time()
for(iter in 1:n_iter){
  if(iter %% (n_iter/10) == 0){
    print(iter)
    print(proc.time() - ptm)
  }
  
  ##### historical regression #####
  
  # first sample latent continuous variables that underly binary response
  z_bf_latent_cur <- sample_z_latent(z_bf,mu_bf_cur%*%t(one_vec_N) +lambda_bf_cur%*%eta_bf_cur,obs_grid_bf)
  
  # sampler mean
  z_star_cur <- compute_y_st(z_bf_latent_cur,lambda_bf_cur%*%eta_bf_cur,obs_grid_bf)
  mu_bf_cur <- sample_mu_sparse(z_star_cur, 1, age_grid, l_mu_bf_current,scale_mu_bf_current,obs_grid_bf)
  
  out <- l_mu_MH(z_star_cur,age_grid,obs_grid_bf,l_mu_bf_current,scale_mu_bf_current,mu_bf_cur,1,proposal_l_mu_bf,1,iter,prior_a,prior_b) 
  l_mu_bf_current <- out[1]
  accept_l_mu_bf_count <- accept_l_mu_bf_count + out[2]
  proposal_l_mu_bf <- out[3]
  
  out <- scale_mu_MH(z_star_cur, age_grid,obs_grid_bf,l_mu_bf_current, scale_mu_bf_current,mu_bf_cur,1,proposal_scale_mu_bf,1,iter,prior_var)
  scale_mu_bf_current<- out[1]
  accept_scale_mu_bf_count <- accept_scale_mu_bf_count + out[2]
  proposal_scale_mu_bf <- out[3]
  
  z_star_cur <- compute_y_st(z_bf_latent_cur,mu_bf_cur%*%t(one_vec_N),obs_grid_bf)
  lambda_bf_cur <- sample_lambda_sparse(z_star_cur,W,age_grid,l_bf_param_cur,scale_bf_param_cur,lambda_bf_cur,eta_bf_cur,1,nu_param,obs_grid_bf)
  
  out <- l_param_MH(l_bf_param_cur,proposal_l_bf,scale_bf_param_cur,age_grid,obs_grid_bf,W,z_star_cur,lambda_bf_cur,eta_bf_cur,1,nu_param,1,iter,prior_a,prior_b)
  l_bf_param_cur <- out[,1]
  accept_l_bf_count <- accept_l_bf_count + out[,2]
  proposal_l_bf <- out[,3]
  
  out <- scale_param_MH(scale_bf_param_cur,proposal_scale_bf,l_bf_param_cur,age_grid,obs_grid_bf,W,z_star_cur,lambda_bf_cur,eta_bf_cur,1,nu_param,1,iter,prior_var)
  scale_bf_param_cur <- out[,1]
  accept_scale_bf_count <- accept_scale_bf_count + out[,2]
  proposal_scale_bf <- out[,3]
  
  eta_bf_cur <- sample_eta_cholupdate(z_star_cur,eta_bf_cur,psi_bf_cur,lambda_bf_cur,1,nu_eta,obs_grid_bf)
  if(PX){
    psi_bf_cur <- sample_psi_sparse(eta_bf_cur,nu_eta,px_a,px_b)
  }
  
  # sample completed version bf covariate 
  z_bf_complete <- sample_z_complete(z_bf,mu_bf_cur%*%t(one_vec_N) + lambda_bf_cur%*%eta_bf_cur,obs_grid_bf)
  
  
  ##### save current state of Markov chain #####
  if(iter%%lag == 0){
    
    iter_save <- iter/lag
    #mu_y_save[,iter_save] <- mu_y_cur
    scale_mu_y_save[iter_save] <- scale_mu_y_current
    l_mu_y_save[iter_save] <- l_mu_y_current
    lambda_y_save[,,iter_save] <- lambda_y_cur
    xi_save[,,iter_save] <- xi_cur
    eta_y_save[,,iter_save] <- eta_y_cur
    psi_y_save[,iter_save] <- psi_y_cur
    l_y_param_save[,iter_save] <- l_y_param_cur
    scale_y_param_save[,iter_save] <- scale_y_param_cur
    
    b_save[,,iter_save] <- b_cur
    
    mu_bf_save[,iter_save] <- mu_bf_cur
    scale_mu_bf_save[iter_save] <- scale_mu_bf_current
    l_mu_bf_save[iter_save] <- l_mu_bf_current
    lambda_bf_save[,,iter_save] <- lambda_bf_cur
    eta_bf_save[,,iter_save] <- eta_bf_cur
    psi_bf_save[,iter_save] <- psi_bf_cur
    l_bf_param_save[,iter_save] <- l_bf_param_cur
    scale_bf_param_save[,iter_save] <- scale_bf_param_cur
    beta_bf_save[,iter_save] <- beta_bf_cur
    
    mu_ill_save[,,iter_save] <- mu_ill_cur
    l_mu_ill_save[,iter_save] <- l_mu_ill_current
    scale_mu_ill_save[,iter_save] <- scale_mu_ill_current
    beta_ill_coefs_save[,,iter_save] <- beta_ill_coefs_cur
    
    sig_sq_save[iter_save] <- sig_sq_cur  
  }
  
}
print(proc.time() - ptm)


# resolve label switching an sign ambiguity
one_M <- rep(1,M)
one_N <- rep(1,N)
lambda_pivot <- lambda_bf_save[,,n_save]*(one_M%*%t(sqrt(psi_bf_save[,n_save])))
eta_pivot <- eta_bf_save[,,n_save]/(sqrt(psi_bf_save[,n_save])%*%t(one_N))
eta_pivot_sd <- sqrt(apply(eta_pivot,1,var))

lambda_pivot <- lambda_pivot*(one_M%*%t(eta_pivot_sd))
eta_pivot <- eta_pivot/(eta_pivot_sd%*%t(one_N))
# after selecting pivot, align mcmc samples with pivot
post_burn_iter <- 1:n_save
lambda_save_processed <- array(0,dim = c(M,K,length(post_burn_iter)))
for(iter in 1:length(post_burn_iter)){
  lambda_post <- lambda_bf_save[,,iter]*(one_M%*%t(sqrt(psi_bf_save[,iter])))
  eta_post <- eta_bf_save[,,iter]/(sqrt(psi_bf_save[,iter])%*%t(one_N))
  eta_sd <- sqrt(apply(eta_post,1,var))
  
  lambda_post <- lambda_post*(one_M%*%t(eta_sd))
  eta_post <- eta_post/(eta_sd%*%t(one_N))
  
  perm_out <- msfOUT_functional(lambda_post,W,lambda_pivot)
  lambda_save_processed[,,iter] <- aplr(lambda_post,perm_out)
}

# Compute family wise 95% credible intervales and see if 0 function is contained
level <- .05/K
lambda_cred_band <- array(0,dim = c(M,2,K)) 
post_burn_iter <- round(seq(n_save/2+1,n_save,1))
factor_active <- rep(0,K)
for(k in 1:K){
  lambda_cred_band[,,k] <- credBands(t(lambda_save_processed[,k,post_burn_iter]),alpha = level)
  factor_active[k] <- !(all(lambda_cred_band[,1,k] <0) & all(lambda_cred_band[,2,k] >0))
  
}
factor_active <- as.logical(factor_active)

for(k in 1:K){
  lambda_cred_band[,,k] <- credBands(t(lambda_save_processed[,k,post_burn_iter]),alpha = level)
  factor_active[k] <- !(all(lambda_cred_band[,1,k] <0) & all(lambda_cred_band[,2,k] >0))
}

print(sum(factor_active))

# generate scree plots
load(paste0("cebu_bffa_select_K.Rdata"))
bf_scale_norm <- array(0,dim = c(K,n_save))
for(iter in 1:n_save){
  lambda_temp <- lambda_save_processed[,,iter]
  bf_scale_norm[,iter] <- sort(sqrt(diag(t(lambda_temp)%*%W%*%lambda_temp)),decreasing = T)
}
norm_df <- melt(bf_scale_norm[,(n_save/2):n_save])
#scree plot
ggplot(norm_df) + geom_boxplot(aes(x = as.factor(Var1),y = value))+ 
  ggtitle("") + xlab("k") + ylab("norm") + theme(text = element_text(size = 24),legend.position = "none") 

