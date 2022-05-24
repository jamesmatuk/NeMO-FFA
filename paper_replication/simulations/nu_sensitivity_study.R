libery(reshape2)
library(Rcpp)
library(coda)

# source c++ files
sourceCpp("../src/remo_fpca.cpp")
 
# nu_lambda values used in sensitivity analysis
nu_vals <- c(100,1,1e-2,1e-4,1e-6,1e-8)
ess_nu <- array(0,dim = c(2,length(nu_vals)))
nu_title <- c(expression(nu[lambda] == 1e2),expression(nu[lambda] == 1e0),expression(nu[lambda] == 1e-2),
              expression(nu[lambda] == 1e-4),expression(nu[lambda] == 1e-6),expression(nu[lambda] == 1e-8))
for(nu_ind in 1:length(nu_vals)){

set.seed(1234)
  
# generate some data where the lambda are orthogonal

K_true <- 2
N <- 100

M <- 30
time_grid <- seq(-pi,pi,length.out = M)
time <- time_grid

M <- length(time_grid)
w <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
W <- diag(w)

lambda_true <- cbind(sin(time_grid),.5*sin(2*time_grid))

eta_true <- matrix(rnorm(K_true*N),nrow = K_true,ncol = N)

f <- lambda_true%*%eta_true

sig_true <- .5
noise <- matrix(rnorm(N*M,0,sig_true),nrow = M,ncol = N)
y <- f + noise

sparsity_level <- .5
obs_grid <- matrix(sample(c(0,1),M*N,replace = T,prob = c(sparsity_level,1-sparsity_level)),nrow = M,ncol = N)
while(prod(!apply(obs_grid,2,sum) <= 1) == 0){
  obs_grid <- matrix(sample(c(0,1),M*N,replace = T,prob = c(sparsity_level,1-sparsity_level)),nrow = M,ncol = N)
}
y_common <- y*obs_grid


##### run MCMC #####
one_vec_N <- rep(1,N)
K <- 5
nu_eta <- 1e-4
nu_param <- nu_vals[nu_ind]


lambda_cur <- matrix(rnorm(M*K),nrow = M,ncol = K)
eta_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
psi_cur <- rep(1,K)
sig_sq_cur <- var((y-lambda_cur%*%eta_cur)[obs_grid>0])
l_mu_current <- .1
scale_mu_current <- 1
l_param_cur <- rep(.1,K)
scale_param_cur <- rep(1,K)

prior_a <- 12
prior_b <- 12
prior_var <- 1

px_b <- 10
px_a <- px_b + 1

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
n_iter <- 10000
mu_save <- matrix(0,nrow = M,ncol = n_iter)
scale_mu_save <- rep(0,n_iter)
l_mu_save <- rep(0,n_iter)
lambda_save <- array(0,dim = c(M,K,n_iter))
eta_save <- array(0,dim = c(K,N,n_iter))
psi_save <- array(0,dim = c(K,n_iter))
l_param_save <- array(0,dim = c(K,n_iter))
scale_param_save <- array(0,dim = c(K,n_iter))
sig_sq_save <- rep(0,n_iter)

for(iter in 1:n_iter){
  
  y_st_cur <- compute_y_st(y_common,lambda_cur%*%eta_cur,obs_grid)
  mu_cur <- sample_mu_sparse(y_st_cur, sig_sq_cur, time, l_mu_current,scale_mu_current,obs_grid)
  
  out <- l_mu_MH(y_st_cur,time,obs_grid,l_mu_current,scale_mu_current,mu_cur,sig_sq_cur,proposal_l_mu,1,iter,prior_a,prior_b) 
  l_mu_current <- out[1]
  accept_l_mu_count <- accept_l_mu_count + out[2]
  proposal_l_mu <- out[3]
  
  out <- scale_mu_MH(y_st_cur, time,obs_grid,l_mu_current, scale_mu_current,mu_cur,sig_sq_cur,proposal_scale_mu,1,iter,prior_var)
  scale_mu_current<- out[1]
  accept_scale_mu_count <- accept_scale_mu_count + out[2]
  proposal_scale_mu<- out[3]
  
  y_st_cur <- compute_y_st(y_common,mu_cur%*%t(one_vec_N),obs_grid)
  lambda_cur <- sample_lambda_sparse(y_st_cur,W,time,l_param_cur,scale_param_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,obs_grid)
  
  out <- l_param_MH(l_param_cur,proposal_l,scale_param_cur,time,obs_grid,W,y_st_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,1,iter,prior_a,prior_b)
  l_param_cur <- out[,1]
  accept_l_count <- accept_l_count + out[,2]
  proposal_l <- out[,3]
  
  out <- scale_param_MH(scale_param_cur,proposal_scale,l_param_cur,time,obs_grid,W,y_st_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,1,iter,prior_var)
  scale_param_cur <- out[,1]
  accept_scale_count <- accept_scale_count + out[,2]
  proposal_scale <- out[,3]
  
  eta_cur <- sample_eta_sparse(y_st_cur,eta_cur,psi_cur,lambda_cur,sig_sq_cur,nu_eta,obs_grid)
  psi_cur <- sample_psi_sparse(eta_cur,nu_eta,px_a,px_b)
  
  sig_sq_cur <- sample_sig_sq_sparse(y_st_cur,lambda_cur,eta_cur,obs_grid)
  
  mu_save[,iter] <- mu_cur
  scale_mu_save[iter] <- scale_mu_current
  l_mu_save[iter] <- l_mu_current
  lambda_save[,,iter] <- lambda_cur
  eta_save[,,iter] <- eta_cur
  psi_save[,iter] <- psi_cur
  l_param_save[,iter] <- l_param_cur
  scale_param_save[,iter] <- scale_param_cur
  sig_sq_save[iter] <- sig_sq_cur
  
}

##### Figures #####

# plot data 
y_long <- NULL
for(i in 1:N){
  id_temp <- rep(i,sum(obs_grid[,i]))
  t_temp <- time[which(obs_grid[,i]>0)]
  y_temp <- y_common[which(obs_grid[,i]>0),i]
  y_long <- rbind(y_long,cbind(id_temp,t_temp,y_temp))
}


lambda_true_df <- as.data.frame(cbind(time,lambda_true*rep(1,M)%*%t(apply(eta_true,1,sd))))
colnames(lambda_true_df)[2:3] <- 1:2
lambda_true_df <- melt(lambda_true_df,id = "time")
colnames(lambda_true_df)[2] <- 'k'
lambda_true_df[,2] <- as.factor(lambda_true_df[,2])

# post - process samples
lambda_post <- lambda_save
eta_post <- eta_save
M <- dim(lambda_save[,,1])[1]
K <- dim(lambda_save[,,1])[2]
N <- dim(eta_save[,,1])[2]
one_M <- rep(1,M)
one_N <- rep(1,N)
for(iter in 1:n_iter){
  lambda_post[,,iter] <- lambda_save[,,iter]*(one_M%*%t(sqrt(psi_save[,iter])))
  eta_post[,,iter] <- eta_save[,,iter]/(sqrt(psi_save[,iter])%*%t(one_N))
}

# resolve label switching an sign ambiguity
sourceCpp("../src/msf.cpp")
post_burn_iter <- 1:n_iter
lambda_save_processed <- array(0,dim = c(M,K,length(post_burn_iter)))
eta_save_processed <- array(0,dim = c(K,N,length(post_burn_iter)))
l_param_save_processed <- array(0,dim = c(K,length(post_burn_iter)))
scale_param_save_processed <- array(0,dim = c(K,length(post_burn_iter)))
for(iter in 1:length(post_burn_iter)){
  perm_out <- msfOUT_functional(lambda_post[,,post_burn_iter[iter]],W,lambda_post[,,round(n_iter/(1.1))])
  #perm_out <- 1:5
  lambda_save_processed[,,iter] <- aplr(lambda_post[,,post_burn_iter[iter]],perm_out)
  eta_save_processed[,,iter] <- t(aplr(t(eta_post[,,post_burn_iter[iter]]),perm_out))
  l_param_save_processed[,iter] <- l_param_save[abs(perm_out),post_burn_iter[iter]]
  scale_param_save_processed[,iter] <- scale_param_save[abs(perm_out),post_burn_iter[iter]]
}

# arrange lambda related parameters by magnitude of lambda
lambda_post_mean <- apply(lambda_save_processed,c(1,2),mean)
lambda_mean_norm <- sqrt(diag(t(lambda_post_mean)%*%W%*%lambda_post_mean))
lambda_mean_order <- order(lambda_mean_norm,decreasing = T)

lambda_save_processed <- lambda_save_processed[,lambda_mean_order,]
eta_save_processed <- eta_save_processed[lambda_mean_order,,]
l_param_save_processed <- l_param_save_processed[lambda_mean_order,]
scale_param_save_processed <- scale_param_save_processed[lambda_mean_order,]


# plot norms 
lambda_save_norm <- array(0,dim = c(K,n_iter))
for(iter in 1:n_iter){
  temp_eta_sd <- apply(eta_save_processed[,,iter],1,sd)
  lambda_temp <- lambda_save_processed[,,iter]*rep(1,M)%*%t(temp_eta_sd)
  lambda_save_norm[,iter] <- sort(sqrt(diag(t(lambda_temp)%*%W%*%lambda_temp)),decreasing = T)
}

norm_df <- t(lambda_save_norm)
colnames(norm_df) <- 1:K
norm_df <- melt(norm_df)
colnames(norm_df)[2] <- "k"
norm_df[,2] <- as.factor(norm_df[,2] )

yl <- c(0,2)
print(
ggplot(norm_df) + geom_line(aes(x = Var1,y = value,color = k)) + 
  ylab("norm") + xlab("iteration") + theme(text = element_text(size = 24)) + 
  ggtitle(nu_title[nu_ind])+ylim(yl) + 
  geom_hline(yintercept = sqrt(diag(t(lambda_true)%*%W%*%lambda_true))*apply(eta_true,1,sd))
)

ess_nu[1,nu_ind] <- effectiveSize(lambda_save_norm[1,(n_iter/2):n_iter])
ess_nu[2,nu_ind] <- effectiveSize(lambda_save_norm[2,(n_iter/2):n_iter])
print(ess_nu)

}

# plot effective sample size
ess_df <- t(ess_nu)
colnames(ess_df) <- 1:2
ess_df <- melt(ess_df)
colnames(ess_df)[2] <- 'k'
ess_df[,2] <- as.factor(ess_df[,2])
print(
ggplot(ess_df) + geom_line(aes(x = Var1,y = value,color = k),size = 2) + 
  xlab(expression(nu[lambda])) + ylab("ESS") + 
  scale_x_continuous(breaks = 1:6,labels = c('100','1','0.01','1e-04','1e-06','1e-08')) + 
  theme(text = element_text(size = 24)) 
)
