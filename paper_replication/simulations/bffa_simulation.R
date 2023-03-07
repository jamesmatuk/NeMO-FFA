library(MASS)
library(Rcpp)
library(pracma)
library(boot)
library(SLFPCA)
library(fda)

sourceCpp("../src/nemo_ffa.cpp")
source("../src/nemo_ffa.R")
sourceCpp("../src/msf.cpp")
sourceCpp("./supp_funs/registr_src.cpp")
source("./supp_funs/registr_funs.R")

# scale setting - either "low" or "high"
scale_setting <- "high"

##### Parameters that can be varied for simulation experiments #####

set.seed(1234)

n_rep <- 100
MISE_mat <- array(0,dim = c(3,3,n_rep))
K_nemo <- array(0,dim = c(3,n_rep))
sparsity_levels <- c(.25,.5,.75)
for(replicate in 1:n_rep){

N <- 100
K_true <- 2
lambda_l <- .4*rep(1,K_true)
lambda_scale <- 1*rep(1,K_true)
mu_l <- .4
mu_scale <- 1 
sig_true <- 1

##### generate underlying observations #####

M <- 30 
one_M <- rep(1,M)
one_N <- rep(1,N)
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

nu_eta <- 1e-4
eta_cov <- diag(one_N) - one_N%*%t(one_N)/(nu_eta + N)
eta_true <- mvrnorm(K_true,mu = rep(0,N),Sigma = eta_cov)
eta_sd <- sqrt(apply(eta_true,1,var))

lambda_true <- lambda_true*(one_M%*%t(eta_sd))
eta_true <- eta_true/(eta_sd%*%t(one_N))

lambda_norm <- sqrt(diag(t(lambda_true)%*%W%*%lambda_true))
if(scale_setting == "low"){lambda_scale <- c(2,1)}
if(scale_setting == "high"){lambda_scale <- c(3,1.5)}

lambda_true <- lambda_true*(one_M%*%t(lambda_scale/lambda_norm))

lambda_norm_order <- order(lambda_norm,decreasing = T)
lambda_true <- lambda_true[,lambda_norm_order]
eta_true <- eta_true[lambda_norm_order,]

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

slfpca_results <- SLFPCA(Ly,Lt,interval = c(min(time),max(time)),npc = 2, L_list = 13,
                         norder = 4, kappa_theta = 0.2, sparse_pen = 0,
                         nRegGrid = length(time), stepmu = 0.005)

slfpca_mu <- eval.fd(time,slfpca_results$mufd)
K <- K_true
slfpca_eig_vec <- array(0,dim = c(M,K)) 
for(k in 1:K){
  slfpca_eig_vec[,k] <- eval.fd(time,slfpca_results$eigfd_list[[k]])
}
slfpca_score <- slfpca_results$score

fit_slfpca <- inv.logit(slfpca_mu%*%t(rep(1,N)) + slfpca_eig_vec%*%t(slfpca_score))


##### run bfpca from the registr package #####
z_bfpca <- data.frame(id = z_long[,1], value = z_long[,3],index = z_long[,2],4)
bfpca_results <- bfpca2(z_bfpca,npc = 2)

fit_bfpca <- inv.logit(bfpca_results$mu%*%t(rep(1,N)) + bfpca_results$efunctions%*%t(bfpca_results$scores))
##### nemo bffa #####
K_max <- 5

nu_eta <- 1e-4
nu_lambda <- 1e-4
inv_g_hyper <- length_scale_hyper(time)
prior_a <- inv_g_hyper[1]
prior_b <- inv_g_hyper[2]
prior_var <- 1
n_iter <- 1000

init_mcmc_params <- init_bffa_mcmc(z_common,obs_grid,time,
                                   K_max,nu_eta,nu_lambda,
                                   prior_a,prior_b,prior_var,
                                   n_iter)

if(dim(init_mcmc_params$lambda_cur)[2] == 0){
  
  init_mcmc_params$lambda_cur <- matrix(0,nrow = M,ncol = 1)
  init_mcmc_params$eta_cur <- matrix(rnorm(K*N,0,1),nrow = 1,ncol = N)
  init_mcmc_params$psi_cur <- rep(1,1)
  init_mcmc_params$l_param_cur <- rep(.1,1)
  init_mcmc_params$scale_param_cur <- rep(1,1)
  init_mcmc_params$proposal_l <- rep(.001,1)
  init_mcmc_params$proposal_scale <- rep(.001,1)
}

lag <- 1
bffa_mcmc <- run_bffa_mcmc(z_common,obs_grid,time,nu_eta,nu_lambda,
                           init_mcmc_params,prior_a,prior_b,prior_var,n_iter,lag)




K_nemo[sparsity_ind,replicate] <- dim(init_mcmc_params$lambda_cur)[2]

fit_save <- array(0,dim = c(M,N,n_iter/2))
for(iter in (n_iter/2 + 1):n_iter){
  if(K_nemo[sparsity_ind,replicate] == 1){
    fit_save[,,iter - n_iter/2] <- pnorm(bffa_mcmc$mu_save[,iter]%*%t(rep(1,N)) + bffa_mcmc$lambda_save[,,iter]%*%t(bffa_mcmc$eta_save[,,iter]))
  }else{
    fit_save[,,iter - n_iter/2] <- pnorm(bffa_mcmc$mu_save[,iter]%*%t(rep(1,N)) + bffa_mcmc$lambda_save[,,iter]%*%bffa_mcmc$eta_save[,,iter])  
  }
  
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
print(K_nemo[,replicate])
print(MISE_mat[,,replicate])

}

##### Visualize results of simulation #####
cube_diff <- function(cube_array){
  diff_array <- array(0,dim= c(3,2,100))
  for(method in 1:2){
    diff_array[,method,] <- cube_array[,1+method,] - cube_array[,1,]
  }
  return(diff_array)
}

reorg_results <- function(cube_array){
  s_levels <- c(.25,.5,.75)
  m_levels <- c("WZSG2019","ZLLZ2021")
  long_rep <- NULL
  long_s <- NULL
  long_m <- NULL
  long_value <- NULL
  for(s_ind in 1:3){
    for(m_ind in 1:2){
      long_rep <- c(long_rep,1:100)
      long_s <- c(long_s,rep(s_levels[s_ind],100))
      long_m <- c(long_m,rep(m_levels[m_ind],100))
      long_value <- c(long_value,cube_array[s_ind,m_ind,])
    }
  }
  long_df <- data.frame(replicate = long_rep,sparsity = long_s,method = long_m,value = long_value)
  return(long_df)
}

temp_df <- reorg_results(cube_diff(MISE_mat))
scale_breaks <- round(seq(0,quantile(temp_df$value, c(0.95)),length.out = 4),2)
scale_labels <- as.character(scale_breaks)  
scale_labels[1] <- "NeMO"

png(file=paste0('bffa_sim_',scale_setting,'.png'), width=400, height=400)
print(
ggplot(temp_df) + geom_boxplot(aes(x = as.factor(sparsity),y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(0,quantile(temp_df$value, c(0.95))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("MISE") + xlab("sparsity level") + theme(text = element_text(size = 24)) + ggtitle(expression(h(f)))
)
dev.off()