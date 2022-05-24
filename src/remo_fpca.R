##### MCMC wrapper function for FPCA using ReMO processes #####

remo_fpca_mcmc <- function(y_common,obs_grid,time_grid,
                           K,nu_eta,nu_lambda,
                           prior_a,prior_b,prior_var,
                           n_iter,lag){
  
# description: MCMC wrapper function for FPCA using ReMO processes
# 
# arguments:
#   y_common - sparse and noisy functional observations
#   obs_grid - binary matrix indicating subject-specific observation grid
#   time_grid - common time grid
#   K - assumed number of FPCs
#   nu_eta - penaly parameter to enforce sum-to-zero constraint for FPC scores
#   nu_lambda - penaly parameter to enforce orthogonality constraint for FPCs
#   prior_a - length-scale hyperparameter for FPCs and mean process
#   prior_b - length-scale hyperparameter for FPCs and mean process
#   prior_var - scale hyperparameter for FPCs and mean process
#   n_iter - number of MCMC iterations to run 
#   lag - lag for saving MCMC samples
#
# output: MCMC samples 
#
 
# number of timepoints on common grid
M <- dim(y_common)[1] 
# number of observations
N <- dim(y_common)[2] 
one_vec_N <- rep(1,N)

# set up weights used riemann sum Riemann sum
w <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
W <- diag(w)


# initialize current parameter values for MCMC
lambda_cur <- matrix(rnorm(M*K),nrow = M,ncol = K)
eta_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
psi_cur <- rep(1,K)
sig_sq_cur <- var((y-lambda_cur%*%eta_cur)[obs_grid>0])
l_mu_cur <- .1
scale_mu_cur <- 1
l_param_cur <- rep(.1,K)
scale_param_cur <- rep(1,K)

# prior parameters for parameter expansion
px_b <- 10
px_a <- px_b + 1


# initialize transition kernel proposal parameters for Metropolis-Hastings (MH) steps
proposal_l_mu <- .0001
proposal_scale_mu <- .0001
proposal_l <- rep(.001,K)
proposal_scale <- rep(.001,K)
accept_l_mu_count <- 0
accept_scale_mu_count <- 0
accept_l_count <- rep(0,K)
accept_scale_count <- rep(0,K)

# pre-allocate memory for MCMC samples
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
  
  # Gibbs step for mean process
  y_st_cur <- compute_y_st(y_common,lambda_cur%*%eta_cur,obs_grid)
  mu_cur <- sample_mu_sparse(y_st_cur, sig_sq_cur, time_grid, l_mu_cur,scale_mu_cur,obs_grid)
  
  # MH step for length-scale hyperparameter of mean process
  out <- l_mu_MH(y_st_cur,time_grid,obs_grid,l_mu_cur,scale_mu_cur,mu_cur,sig_sq_cur,proposal_l_mu,1,iter,prior_a,prior_b) 
  l_mu_cur <- out[1]
  accept_l_mu_count <- accept_l_mu_count + out[2]
  proposal_l_mu <- out[3]
  
  # MH step for scale hyperparameter of mean process
  out <- scale_mu_MH(y_st_cur, time_grid,obs_grid,l_mu_cur, scale_mu_cur,mu_cur,sig_sq_cur,proposal_scale_mu,1,iter,prior_var)
  scale_mu_cur<- out[1]
  accept_scale_mu_count <- accept_scale_mu_count + out[2]
  proposal_scale_mu<- out[3]
  
  # Gibbs step for FPCs
  y_st_cur <- compute_y_st(y_common,mu_cur%*%t(one_vec_N),obs_grid)
  lambda_cur <- sample_lambda_sparse(y_st_cur,W,time_grid,l_param_cur,scale_param_cur,lambda_cur,eta_cur,sig_sq_cur,nu_lambda,obs_grid)
  
  # MH step for length-scale hyperparameter of mean process
  out <- l_param_MH(l_param_cur,proposal_l,scale_param_cur,time_grid,obs_grid,W,y_st_cur,lambda_cur,eta_cur,sig_sq_cur,nu_lambda,1,iter,prior_a,prior_b)
  l_param_cur <- out[,1]
  accept_l_count <- accept_l_count + out[,2]
  proposal_l <- out[,3]
  
  # MH step for scale hyperparameter of FPCs
  out <- scale_param_MH(scale_param_cur,proposal_scale,l_param_cur,time_grid,obs_grid,W,y_st_cur,lambda_cur,eta_cur,sig_sq_cur,nu_lambda,1,iter,prior_var)
  scale_param_cur <- out[,1]
  accept_scale_count <- accept_scale_count + out[,2]
  proposal_scale <- out[,3]
  
  # Gibbs step for FPC scores
  eta_cur <- sample_eta_sparse(y_st_cur,eta_cur,psi_cur,lambda_cur,sig_sq_cur,nu_eta,obs_grid)
  
  # Gibbs step for parameter expansion scale parameter 
  psi_cur <- sample_psi_sparse(eta_cur,nu_eta,px_a,px_b)
  
  # Gibbs step for error variance
  sig_sq_cur <- sample_sig_sq_sparse(y_st_cur,lambda_cur,eta_cur,obs_grid)
  
  # save current state of MCMC (depending on lag)
  if((iter %% lag) == 0){
    iter_save <- iter/lag
    mu_save[,iter_save] <- mu_cur
    scale_mu_save[iter_save] <- scale_mu_cur
    l_mu_save[iter_save] <- l_mu_cur
    lambda_save[,,iter_save] <- lambda_cur
    eta_save[,,iter_save] <- eta_cur
    psi_save[,iter_save] <- psi_cur
    l_param_save[,iter_save] <- l_param_cur
    scale_param_save[,iter_save] <- scale_param_cur
    sig_sq_save[iter_save] <- sig_sq_cur
  }
  
}

# store saved mcmc samples in a list
mcmc_output <- list()

mcmc_output$mu_save <- mu_save
mcmc_output$scale_mu_save <- scale_mu_save
mcmc_output$l_mu_save <- l_mu_save
mcmc_output$lambda_save <- lambda_save
mcmc_output$eta_save <- eta_save
mcmc_output$psi_save <- psi_save
mcmc_output$l_param_save <- l_param_save
mcmc_output$scale_param_save <- scale_param_save
mcmc_output$sig_sq_save <- sig_sq_save

# return saved mcmc samples
return(mcmc_output)
}

##### MCMC wrapper function for binary  FPCA using ReMO processes #####

remo_binary_fpca_mcmc <- function(z_common,obs_grid,time_grid,
                      K,nu_eta,nu_lambda,
                      prior_a,prior_b,prior_var,
                      n_iter,lag){
  
# description: MCMC wrapper function for binary FPCA using ReMO processes
# 
# arguments:
#   z_common - sparse and noisy binary functional observations
#   obs_grid - binary matrix indicating subject-specific observation grid
#   time_grid - common time grid
#   K - assumed number of FPCs
#   nu_eta - penaly parameter to enforce sum-to-zero constraint for FPC scores
#   nu_lambda - penaly parameter to enforce orthogonality constraint for FPCs
#   prior_a - length-scale hyperparameter for FPCs and mean process
#   prior_b - length-scale hyperparameter for FPCs and mean process
#   prior_var - scale hyperparameter for FPCs and mean process
#   n_iter - number of MCMC iterations to run 
#   lag - lag for saving MCMC samples
#
# output: MCMC samples 
#  
  
# number of time points on common grid
M <- length(time_grid)  
# number of observations
N <- dim(z_common)[2]
one_vec_N <- rep(1,N)
# compute weights for Riemann sum
w <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
W <- diag(w)
# parameter expansion hyperparameters
px_b <- 10
px_a <- px_b + 1

# initialize current parameter states for MCMC
mu_cur <- rep(0,M)
lambda_cur <- matrix(0,nrow = M,ncol = K)
eta_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
psi_cur <- rep(1,K)
sig_sq_cur <- 1
l_mu_cur <- 1
scale_mu_cur <- 1
l_param_cur <- rep(1,K)
scale_param_cur <- rep(1,K)

# initialize transition kernel proposal parameters for Metropolis-Hastings (MH) steps
proposal_l_mu <- .0001
proposal_scale_mu <- .0001
proposal_l <- rep(.001,K)
proposal_scale <- rep(.001,K)
accept_l_mu_count <- 0
accept_scale_mu_count <- 0
accept_l_count <- rep(0,K)
accept_scale_count <- rep(0,K)

# pre-allocate memory for MCMC samples
n_save <- n_iter/lag
mu_save <- matrix(0,nrow = M,ncol = n_save)
scale_mu_save <- rep(0,n_save)
l_mu_save <- rep(0,n_save)
lambda_save <- array(0,dim = c(M,K,n_save))
eta_save <- array(0,dim = c(K,N,n_save))
psi_save <- array(0,dim = c(K,n_save))
l_param_save <- array(0,dim = c(K,n_save))
scale_param_save <- array(0,dim = c(K,n_save))

for(iter in 1:n_iter){
  
  # Gibbs step for latent continuous process 
  z_latent_cur <- sample_z_latent(z_common,mu_cur%*%t(one_vec_N) +lambda_cur%*%eta_cur,obs_grid)
  
  # Gibbs step for mean process
  z_star_cur <- compute_y_st(z_latent_cur,lambda_cur%*%eta_cur,obs_grid)
  mu_cur <- sample_mu_sparse(z_star_cur, sig_sq_cur, time_grid, l_mu_cur,scale_mu_cur,obs_grid)
  
  # MH step for length-scale hyperparameter of mean process
  out <- l_mu_MH(z_star_cur,time_grid,obs_grid,l_mu_cur,scale_mu_cur,mu_cur,sig_sq_cur,proposal_l_mu,1,iter,prior_a,prior_b) 
  l_mu_cur <- out[1]
  accept_l_mu_count <- accept_l_mu_count + out[2]
  proposal_l_mu <- out[3]
  
  # MH step for scale hyperparameter of mean process
  out <- scale_mu_MH(z_star_cur, time_grid,obs_grid,l_mu_cur, scale_mu_cur,mu_cur,sig_sq_cur,proposal_scale_mu,1,iter,prior_var)
  scale_mu_cur<- out[1]
  accept_scale_mu_count <- accept_scale_mu_count + out[2]
  proposal_scale_mu<- out[3]
  
  # Gibbs step for FPCs
  z_star_cur <- compute_y_st(z_latent_cur,mu_cur%*%t(one_vec_N),obs_grid)
  lambda_cur <- sample_lambda_sparse(z_star_cur,W,time_grid,l_param_cur,scale_param_cur,lambda_cur,eta_cur,sig_sq_cur,nu_lambda,obs_grid)
  
  # MH step for length-scale hyperparameter of FPCs
  out <- l_param_MH(l_param_cur,proposal_l,scale_param_cur,time_grid,obs_grid,W,z_star_cur,lambda_cur,eta_cur,sig_sq_cur,nu_lambda,1,iter,prior_a,prior_b)
  l_param_cur <- out[,1]
  accept_l_count <- accept_l_count + out[,2]
  proposal_l <- out[,3]
  
  # MH step for scale hyperparameter of FPCs
  out <- scale_param_MH(scale_param_cur,proposal_scale,l_param_cur,time_grid,obs_grid,W,z_star_cur,lambda_cur,eta_cur,sig_sq_cur,nu_lambda,1,iter,prior_var)
  scale_param_cur <- out[,1]
  accept_scale_count <- accept_scale_count + out[,2]
  proposal_scale <- out[,3]
  
  # Gibbs step FPC scores
  eta_cur <- sample_eta_cholupdate(z_star_cur,eta_cur,psi_cur,lambda_cur,sig_sq_cur,nu_eta,obs_grid)
  
  # Gibbs step for parameter expansion scale parameter 
  psi_cur <- sample_psi_sparse(eta_cur,nu_eta,px_a,px_b)
  
  # save current state of MCMC (depending on lag)
  if(iter%%lag == 0){
    iter_save <- iter/lag
    mu_save[,iter_save] <- mu_cur
    scale_mu_save[iter_save] <- scale_mu_cur
    l_mu_save[iter_save] <- l_mu_cur
    lambda_save[,,iter_save] <- lambda_cur
    eta_save[,,iter_save] <- eta_cur
    psi_save[,iter_save] <- psi_cur
    l_param_save[,iter_save] <- l_param_cur
    scale_param_save[,iter_save] <- scale_param_cur
  }
}


# store saved mcmc samples in a list
mcmc_output <- list()

mcmc_output$mu_save <- mu_save
mcmc_output$scale_mu_save <- scale_mu_save
mcmc_output$l_mu_save <- l_mu_save
mcmc_output$lambda_save <- lambda_save
mcmc_output$eta_save <- eta_save
mcmc_output$psi_save <- psi_save
mcmc_output$l_param_save <- l_param_save
mcmc_output$scale_param_save <- scale_param_save

# return saved mcmc samples
return(mcmc_output)
}

##### organize MCMC samples #####
# resolve label and sign switching
organize_mcmc <- function(mcmc_output,time_grid,pivot_ind){
  
# description: resolve label and sign switching ambiguity in MCMC samples
# 
# arguments:
#   mcmc_output - MCMC samples (output from remo_fpca_mcmc/remo_binary_fpca_mcmc)
#   time_grid - common time grid
#   pivot_ind - MCMC index used as pivot for MSF algorithm
#
# output: MCMC samples with ambiguity resolved 
#  
  

  lambda_post <- mcmc_output$lambda_save
  eta_post <- mcmc_output$eta_save
  psi_post <- mcmc_output$psi_save
  M <- dim(lambda_post[,,1])[1]
  K <- dim(lambda_post[,,1])[2]
  N <- dim(eta_post[,,1])[2]
  n_iter <- dim(lambda_post)[3]
  one_M <- rep(1,M)
  one_N <- rep(1,N)
  # scale lamda and eta by psi
  for(iter in 1:n_iter){
    lambda_post[,,iter] <- lambda_post[,,iter]*(one_M%*%t(sqrt(psi_post[,iter])))
    eta_post[,,iter] <- eta_post[,,iter]/(sqrt(psi_post[,iter])%*%t(one_N))
  }
  
  # resolve label switching an sign ambiguity
  w <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
  W <- diag(w)
  lambda_save_processed <- array(0,dim = c(M,K,n_iter))
  eta_save_processed <- array(0,dim = c(K,N,n_iter))
  l_param_save_processed <- array(0,dim = c(K,n_iter))
  scale_param_save_processed <- array(0,dim = c(K,n_iter))
  for(iter in 1:n_iter){
    perm_out <- msfOUT_functional(lambda_post[,,iter],W,lambda_post[,,pivot_ind])
    lambda_save_processed[,,iter] <- aplr(lambda_post[,,iter],perm_out)
    eta_save_processed[,,iter] <- t(aplr(t(eta_post[,,iter]),perm_out))
    l_param_save_processed[,iter] <- mcmc_output$l_param_save[abs(perm_out),iter]
    scale_param_save_processed[,iter] <- mcmc_output$scale_param_save[abs(perm_out),iter]
  }
  
  # arrange parameters by magnitude of lambda
  lambda_post_mean <- apply(lambda_save_processed,c(1,2),mean)
  lambda_mean_norm <- sqrt(diag(t(lambda_post_mean)%*%W%*%lambda_post_mean))
  lambda_mean_order <- order(lambda_mean_norm,decreasing = T)
  
  lambda_save_processed <- lambda_save_processed[,lambda_mean_order,]
  eta_save_processed <- eta_save_processed[lambda_mean_order,,]
  l_param_save_processed <- l_param_save_processed[lambda_mean_order,]
  scale_param_save_processed <- scale_param_save_processed[lambda_mean_order,]

  mcmc_output$lambda_save <- lambda_save_processed
  mcmc_output$eta_save <- eta_save_processed
  mcmc_output$l_param_save <- l_param_save_processed
  mcmc_output$scale_param_save <- scale_param_save_processed
  
  return(mcmc_output)
}

# compute fitted functions based on MCMC samples of parameters
compute_fit_samples <- function(mcmc_output){
  
# description: computes functions underlying sparse and noisy observations
#              based on mu, lambda, and eta
# 
#   mcmc_output - MCMC samples (output from remo_fpca_mcmc/remo_binary_fpca_mcmc)
#
# output: MCMC samplesof underlying functions (MxNxn_iter dimensional matrix)
#  
  
  M <- dim(mcmc_output$lambda_save)[1]
  N <- dim(mcmc_output$eta_save)[2]
  n_iter <- dim(mcmc_output$lambda_save)[3]
  # sum the inferred mean and FPCs multiplied by FPC scores  
  fit_samples <- array(0,dim = c(M,N,n_iter))
  for(iter in 1:n_iter){
    fit_samples[,,iter] <- mcmc_output$mu_save[,iter]%*%t(rep(1,N)) + mcmc_output$lambda_save[,,iter]%*%mcmc_output$eta_save[,,iter]
  }
  return(fit_samples)
}

##### MCMC visualization functions #####

# plot sparse and noisy funcional observations (spagetti plot)
plot_data <- function(time_grid,y_common,obs_grid){
  
# description: plots sparse and noisy functional observations
# 
# arguments:
#   y_common - Sparse and noisy observations on common grid
#   time_grid - common time grid
#   obs_grid - accounts for subject-specific time grids
#
# output: ggplot object (spaghetti plot of functional observations)
#    

  N <- dim(y_common)[2]
  # format long array 
  all_id <- NULL
  all_time <- NULL
  all_value <- NULL
  for(i in 1:N){
    temp_ind <- which(obs_grid[,i] == 1)
    all_id <- c(all_id,rep(i,length(temp_ind)))
    all_time <- c(all_time,time_grid[temp_ind])
    all_value <- c(all_value,y_common[temp_ind,i])
  }
  # display functional data
  cebu_all = data.frame(id = as.factor(all_id), time = all_time, value = all_value)
  ggplot(cebu_all) + geom_line(aes(x = time,y = value,color = id)) + 
    ylab("y(t)") + xlab("t")+ theme(text = element_text(size = 24),legend.position = "none")
}

# plot functional parameters
plot_fun_parameter <- function(fun_param_samples,time_grid,y_label = "value",yl = NULL){
# description: plots thinned mcmc samples of functional paramaters
# 
# arguments:
#   fun_param_samples - (Mxn_iter-dimensional matrix) functional parameter MCMC samples
#   time_grid - common time grid
#   y_label - label for y-axis
#   yl - limits for y-axis
#
# output: ggplot object (spaghetti plot of functional observations)
#    
  
  # thin parameter samples for easier plot rendering
  n_iter <- dim(fun_param_samples)[2]
  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
  param_thin_df <- as.data.frame(cbind(time_grid,fun_param_samples[,thin_iter]))
  param_thin_df <- melt(param_thin_df,id = "time_grid")
  # plot thinned samples
  if(is.null(yl)){yl = c(min(param_thin_df$value),max(param_thin_df$value))}
  ggplot() + 
    geom_line(data = param_thin_df, aes(x = time_grid, y = value, group = variable),alpha = .05,size = 2) +  
    ylab(y_label) + xlab("t") + ylim(yl) + theme(text = element_text(size = 24))
}

# plot binary FPCs
plot_binary_fpcs <- function(mcmc_output,time_grid,k,title_label = NULL){
# description: visualizes binary FPCs thinned mcmc samples of functional paramaters
# 
# arguments:
#   mcmc_output -  MCMC samples (output from remo_fpca_mcmc/remo_binary_fpca_mcmc)
#   time_grid - common time grid
#   k - FPC index to visualize
#   title_label - title for ggplot object
#
# output: ggplot object that visualizes +/- 1 sd from mean 
#         in direction  of FPCs
#    
# make color scheme consistent with ggplot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_colors <- gg_color_hue(2)

M <-  dim(mcmc_output$mu_save)[1]
n_iter <- dim(mcmc_output$mu_save)[2]
thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
# perterb the mean in the direction of FPCs
mu_thin_df <- as.data.frame(cbind(time_grid,mcmc_output$mu_save[,thin_iter]))
mu_thin_df <- melt(mu_thin_df,id = "time_grid")
mean_perterb <- array(0,dim = c(M,length(thin_iter)))
for(iter in 1:length(thin_iter)){
  mean_perterb[,iter] <- mcmc_output$lambda_save[,k,thin_iter[iter]]*sd(mcmc_output$eta_save[k,,thin_iter[iter]])
}
mean_perterb_df <- as.data.frame(cbind(time_grid,mean_perterb))
mean_perterb_df <- melt(mean_perterb_df,id = "time_grid")

# plot perterbations from mean
yl <- c(0,1)
    ggplot() + 
      geom_line(aes(x = mu_thin_df$time_grid, y = pnorm(mu_thin_df$value), group = mu_thin_df$variable),alpha = .05,size = 2) + 
      geom_line(aes(x = mu_thin_df$time_grid, y = pnorm(mu_thin_df$value + mean_perterb_df$value), group = mu_thin_df$variable),alpha = .05,size = 2,color = gg_colors[1]) + 
      geom_line(aes(x = mu_thin_df$time_grid, y = pnorm(mu_thin_df$value - mean_perterb_df$value), group = mu_thin_df$variable),alpha = .05,size = 2,color = gg_colors[2]) + 
      ylab("probability") + xlab("t")+ ggtitle(title_label) +theme(text = element_text(size = 24))  + ylim(yl)
}

##### ancillary functions #####

# automatically determine hyperparameters for length-scale prior
length_scale_hyper <- function(time_grid){
  
# description: selects reasonable hyperparameter values for length-scale
#              parameter in the squared-exponential covariance function              
# 
# arguments:
#   time_grid - common time grid
#
# output: vector of reasonable hyperparameter values for length-scale
#              parameter in the squared-exponential covariance function       
#
  
  M <- length(time_grid)
  # determine cutoffs for length-scale parameter
  lower <- min(diff(time_grid))/3
  upper <- (max(time_grid) - min(time_grid))/3
  # define function to approximate tails of a inverse gamma density
  invgamma_tails_approx <- function(hyperparams){
    c(
      F1 = lower - hyperparams[2]/(hyperparams[1] - 1) + 3*hyperparams[2]/sqrt((hyperparams[1] - 1)^2*(hyperparams[1] - 2)),
      F2 = upper - hyperparams[2]/(hyperparams[1] - 1) - 3*hyperparams[2]/sqrt((hyperparams[1] - 1)^2*(hyperparams[1] - 2))
    )
  }
  # find hyperparameters that satisfy tail probability constraint
  root_out <- multiroot(f = invgamma_tails_approx, start = c(11,10),positive = T)
  return(root_out$root)
}
