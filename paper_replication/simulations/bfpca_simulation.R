library(MASS)
library(Rcpp)
library(pracma)
library(boot)
library(SLFPCA)
library(fda)


sourceCpp("../src/remo_fpca.cpp")
sourceCpp("./registr_funs/registr_src.cpp")
source("./registr_funs/registr_funs.R")

##### Parameters that can be varied for simulation experiments #####

set.seed(1234)

n_rep <- 100
MISE_mat <- array(0,dim = c(3,3,n_rep))
sparsity_levels <- c(.25,.5,.75)
for(replicate in 1:n_rep){

N <- 100
K_true <- 2
lambda_l <- .4*rep(1,K_true)
lambda_scale <- 1*rep(1,K_true)
mu_l <- .4
mu_scale <- 1 # low variability
# mu_scale <- 1 # high variability

sig_true <- 1

##### generate underlying observations #####

M <- 30 
time <- seq(0,2,length.out = M)
w <- diff(c(time[1],(time[2:M]+time[1:(M-1)])/2,time[M]))
W <- diag(w)

mu_cov <- make_cov(time,mu_l,mu_scale) 
mu_true <-  mvrnorm(mu = rep(0,length(time)),Sigma = mu_cov)

prior_cov <- make_cov(time,mu_l,mu_scale) 
lambda_true <- t(mvrnorm(K_true,mu = rep(0,length(time)),Sigma = mu_cov))

nu_param <- 1e-4
for(iter in 1:1000){
  for(k in 1:K_true){
    prior_cov<- make_cov(time,lambda_l[k],lambda_scale[k]) 
    lambda_minus_k <- lambda_true[,-k,drop = F]
    ortho_cov <- prior_cov - 
      prior_cov%*%W%*%lambda_minus_k%*%
      inv(nu_param + t(lambda_minus_k)%*%W%*%prior_cov%*%W%*%lambda_minus_k)%*%
      t(lambda_minus_k)%*%W%*%prior_cov
    lambda_true[,k] <- mvrnorm(mu = rep(0,length(time)),Sigma = ortho_cov)
  }
}

eta_true <- matrix(rnorm(K_true*N),nrow = K_true,ncol = N)

f <- mu_true%*%t(rep(1,N)) + lambda_true%*%eta_true
f_prob <- pnorm(f)
y <- f + matrix(rnorm(N*M,0,sig_true),nrow = M,ncol = N)
z <- apply(y>0,2,as.numeric)


for(sparsity_ind in 1:3){
sparsity_level <-sparsity_levels[sparsity_ind]

obs_grid <- matrix(sample(c(0,1),M*N,replace = T,prob = c(sparsity_level,1-sparsity_level)),nrow = M,ncol = N)
while(prod(!apply(obs_grid,2,sum) <= 1) == 0){
  obs_grid <- matrix(sample(c(0,1),M*N,replace = T,prob = c(sparsity_level,1-sparsity_level)),nrow = M,ncol = N)
}

z_common <- z*obs_grid

z_long <- NULL
for(i in 1:N){
  id_temp <- rep(i,sum(obs_grid[,i]))
  t_temp <- time[which(obs_grid[,i]>0)]
  z_temp <- z_common[which(obs_grid[,i]>0),i]
  z_long <- rbind(z_long,cbind(id_temp,t_temp,z_temp))
}

##### Run SLFPCA#####
Ly <- list()
Lt <- list()
for(i in 1:N){
  inds_obs <- which(obs_grid[,i] > 0)
  Ly[[i]] <-z_common[inds_obs,i]
  Lt[[i]] <- time[inds_obs]
}

slfpca_results <- SLFPCA(Ly,Lt,interval = c(min(time),max(time)),npc = 5, L_list = 13,
                         norder = 4, kappa_theta = 0.2, sparse_pen = 0,
                         nRegGrid = length(time), stepmu = 0.005)

slfpca_mu <- eval.fd(time,slfpca_results$mufd)
K <- 5
slfpca_eig_vec <- array(0,dim = c(M,K)) 
for(k in 1:K){
  slfpca_eig_vec[,k] <- eval.fd(time,slfpca_results$eigfd_list[[k]])
}
slfpca_score <- slfpca_results$score

fit_slfpca <- inv.logit(slfpca_mu%*%t(rep(1,N)) + slfpca_eig_vec%*%t(slfpca_score))


##### run bfpca from the registr package #####
z_bfpca <- data.frame(id = z_long[,1], value = z_long[,3],index = z_long[,2],4)
bfpca_results <- bfpca2(z_bfpca,npc = 5)

fit_bfpca <- inv.logit(bfpca_results$mu%*%t(rep(1,N)) + bfpca_results$efunctions%*%t(bfpca_results$scores))

##### remo bfpca #####

K <- 5
N <- dim(z_common)[2]
one_vec_N <- rep(1,N)
nu_eta <- 1e-4
nu_param <- 1e-4

age_grid <- time
M <- length(age_grid)
w <- diff(c(age_grid[1],(age_grid[2:M]+age_grid[1:(M-1)])/2,age_grid[M]))
W <- diag(w)

prior_a <- 12
prior_b <- 4
prior_var <- 1

mu_cur <- rep(0,M)
lambda_cur <- matrix(0,nrow = M,ncol = K)
eta_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
psi_cur <- rep(1,K)
sig_sq_cur <- 1
l_mu_current <- 1
scale_mu_current <- 1
l_param_cur <- rep(1,K)
scale_param_cur <- rep(1,K)

# pre-specify proposal kernel parameters
proposal_l_mu <- .0001
proposal_scale_mu <- .0001
proposal_l <- rep(.001,K)
proposal_scale <- rep(.001,K)
accept_l_mu_count <- 0
accept_scale_mu_count <- 0
accept_l_count <- rep(0,K)
accept_scale_count <- rep(0,K)

# pre-allocate memory for MCMC samples
n_iter <- 1000
lag <- 1
n_save <- n_iter/lag
mu_save <- matrix(0,nrow = M,ncol = n_save)
scale_mu_save <- rep(0,n_save)
l_mu_save <- rep(0,n_save)
lambda_save <- array(0,dim = c(M,K,n_save))
eta_save <- array(0,dim = c(K,N,n_save))
psi_save <- array(0,dim = c(K,n_save))
l_param_save <- array(0,dim = c(K,n_save))
scale_param_save <- array(0,dim = c(K,n_save))
sig_sq_save <- rep(0,n_save)


for(iter in 1:n_iter){
  
  # first sample latent continuous variables that underly binary response
  z_latent_cur <- sample_z_latent(z_common,mu_cur%*%t(one_vec_N) +lambda_cur%*%eta_cur,obs_grid)
  
  # sampler mean
  z_star_cur <- compute_y_st(z_latent_cur,lambda_cur%*%eta_cur,obs_grid)
  mu_cur <- sample_mu_sparse(z_star_cur, sig_sq_cur, age_grid, l_mu_current,scale_mu_current,obs_grid)
  
  out <- l_mu_MH(z_star_cur,age_grid,obs_grid,l_mu_current,scale_mu_current,mu_cur,sig_sq_cur,proposal_l_mu,1,iter,prior_a,prior_b) 
  l_mu_current <- out[1]
  accept_l_mu_count <- accept_l_mu_count + out[2]
  proposal_l_mu <- out[3]
  
  out <- scale_mu_MH(z_star_cur, age_grid,obs_grid,l_mu_current, scale_mu_current,mu_cur,sig_sq_cur,proposal_scale_mu,1,iter,prior_var)
  scale_mu_current<- out[1]
  accept_scale_mu_count <- accept_scale_mu_count + out[2]
  proposal_scale_mu<- out[3]
  
  z_star_cur <- compute_y_st(z_latent_cur,mu_cur%*%t(one_vec_N),obs_grid)
  lambda_cur <- sample_lambda_sparse(z_star_cur,W,age_grid,l_param_cur,scale_param_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,obs_grid)
  
  out <- l_param_MH(l_param_cur,proposal_l,scale_param_cur,age_grid,obs_grid,W,z_star_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,1,iter,prior_a,prior_b)
  l_param_cur <- out[,1]
  accept_l_count <- accept_l_count + out[,2]
  proposal_l <- out[,3]
  
  out <- scale_param_MH(scale_param_cur,proposal_scale,l_param_cur,age_grid,obs_grid,W,z_star_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,1,iter,prior_var)
  scale_param_cur <- out[,1]
  accept_scale_count <- accept_scale_count + out[,2]
  proposal_scale <- out[,3]
  
  eta_cur <- sample_eta_cholupdate(z_star_cur,eta_cur,psi_cur,lambda_cur,sig_sq_cur,nu_eta,obs_grid)
  psi_cur <- sample_psi_sparse(eta_cur,nu_eta,11,10)
  
  if(iter%%lag == 0){
    iter_save <- iter/lag
    mu_save[,iter_save] <- mu_cur
    scale_mu_save[iter_save] <- scale_mu_current
    l_mu_save[iter_save] <- l_mu_current
    lambda_save[,,iter_save] <- lambda_cur
    eta_save[,,iter_save] <- eta_cur
    psi_save[,iter_save] <- psi_cur
    l_param_save[,iter_save] <- l_param_cur
    scale_param_save[,iter_save] <- scale_param_cur
    sig_sq_save[iter_save] <- sig_sq_cur  
  }
  
}

fit_save <- array(0,dim = c(M,N,n_iter/2))
for(iter in (n_iter/2 + 1):n_iter){
  fit_save[,,iter - n_iter/2] <- pnorm(mu_save[,iter]%*%t(one_vec_N) + lambda_save[,,iter]%*%eta_save[,,iter])
}
fit_mean <- apply(fit_save,c(1,2),mean)

##### compare estimated functions #####

l2_loss_slfpca <- diag(t(f_prob - fit_slfpca)%*%W%*%(f_prob - fit_slfpca))
l2_loss_bfpca <- diag(t(f_prob - fit_bfpca)%*%W%*%(f_prob - fit_bfpca))
l2_loss <- diag(t(f_prob - fit_mean)%*%W%*%(f_prob - fit_mean))

MISE_mat[sparsity_ind,1,replicate] <- mean(l2_loss)
MISE_mat[sparsity_ind,2,replicate] <- mean(l2_loss_bfpca)
MISE_mat[sparsity_ind,3,replicate] <- mean(l2_loss_slfpca)

}

print(replicate)
print(MISE_mat[,,replicate])

}

##### Visualize results of simulation #####

MISE_mean <- apply(MISE_mat,c(1,2),mean)
MISE_sd <- apply(MISE_mat,c(1,2),sd)
MISE_lower<- MISE_mean - 1.96*MISE_sd/sqrt(dim(MISE_mat)[3])
MISE_upper<- MISE_mean + 1.96*MISE_sd/sqrt(dim(MISE_mat)[3])

method_fact <- factor(c(rep("ReMO",3),rep("registr",3),rep("SLFPCA",3)),levels = c("ReMO","registr","SLFPCA"))
mean_df = data.frame(method = method_fact,value = c(MISE_mean))
lower_df = data.frame(method = method_fact,value = c(MISE_lower))
upper_df = data.frame(method = method_fact,value = c(MISE_upper))

sparsity_levels <- c(.25,.5,.75)
yl <- c(min(MISE_lower),max(MISE_upper))
ggplot(mean_df) + geom_line(aes(x = rep(sparsity_levels,3),y = value,color = method),size = 2) + 
  geom_ribbon(mapping = aes(x = rep(sparsity_levels,3),ymin = lower_df$value,ymax = upper_df$value ,fill = method),alpha = .2) + 
  ylab("MISE") + xlab("sparsity level") + ylim(yl) + theme(text = element_text(size = 24))+
  scale_x_continuous(breaks = sparsity_levels)

