library(MASS)
library(pracma)
library(ggplot2)
library(reshape2)
library(Rcpp)

# source c++ files
sourceCpp("../src/remo_fpca.cpp")

########## Simulate Function for Orthogonality ##########

# covariance for simulated observation
set.seed(4444)
M_grid <- 100
time_grid <- seq(0,1,length.out = M_grid)
w_grid <- diff(c(time_grid[1],(time_grid[2:M_grid]+time_grid[1:(M_grid-1)])/2,time_grid[M_grid]))
W_grid <- diag(w_grid)
lambda_l <- .1
lambda_scale <- 1
cov_mat <- make_cov(time_grid,lambda_l,lambda_scale) 

lambda_minus_k <- mvrnorm(mu = rep(0,length(time_grid)),Sigma = cov_mat)

temp_df = as.data.frame(cbind(time_grid,lambda_minus_k))
ggplot(temp_df) + geom_line(aes(x = time_grid,y = lambda_minus_k),size = 2)+ xlab("t") + 
  ylab(expression(lambda[1](t))) + theme(legend.position = "none") + 
  theme(text = element_text(size = 24)) 

########## Vary length-scale ##########

prior_scale <- 1
nu_param <- .0001

params <- c(.001,.01,.1,1)
 
legend_labs <- rep("",length(params))
for(ind in 1:length(legend_labs)){
  legend_labs[ind] <- paste0(params[ind])
}

# Look at some draws from the prior
num_draws <- 100
ortho_draws <- array(0,dim = c(M_grid,num_draws,length(params)))
for(ind in 1:length(params)){
  prior_l <- params[ind]
  prior_cov <- make_cov(time_grid,prior_l,prior_scale) 
  
  ortho_cov <- prior_cov - 
    prior_cov%*%W_grid%*%lambda_minus_k%*%
    inv(nu_param + t(lambda_minus_k)%*%W_grid%*%prior_cov%*%W_grid%*%lambda_minus_k)%*%
    t(lambda_minus_k)%*%W_grid%*%prior_cov
  
  ortho_draws[,,ind] <- t(mvrnorm(num_draws, mu = rep(0,length(time_grid)),Sigma = ortho_cov))
}

yl <- c(min(ortho_draws),max(ortho_draws))
count <- 1
for(ind in 1:length(params)){
  temp_df <- as.data.frame(cbind(time_grid,ortho_draws[,,ind]))
  temp_df <- melt(temp_df,id = "time_grid")
  print(
    ggplot() + 
      geom_line(data = temp_df, aes(x = time_grid, y = value, group = variable),alpha = .1,size = 2) +  
      ylab(expression(lambda[2](t))) + xlab("t")+theme(text = element_text(size = 24))  + ylim(yl)
  )
  count <- count + 1
}


# what happens to inner product
num_draws <- 1000
inner_prod <- array(0,dim = c(num_draws,length(params)))
ortho_draws <- array(0,dim = c(M_grid,num_draws,length(params)))
for(ind in 1:length(params)){
  prior_l <- params[ind]
  prior_cov <- make_cov(time_grid,prior_l,prior_scale) 
  
  ortho_cov <- prior_cov - 
    prior_cov%*%W_grid%*%lambda_minus_k%*%
    inv(nu_param + t(lambda_minus_k)%*%W_grid%*%prior_cov%*%W_grid%*%lambda_minus_k)%*%
    t(lambda_minus_k)%*%W_grid%*%prior_cov
  
  ortho_draws[,,ind] <- t(mvrnorm(num_draws, mu = rep(0,length(time_grid)),Sigma = ortho_cov))
  
  inner_prod[,ind] <- t(ortho_draws[,,ind])%*%W_grid%*%lambda_minus_k
}

temp_df <- as.data.frame(inner_prod)
colnames(temp_df) <- legend_labs
temp_df <- melt(temp_df)
ggplot(data = as.data.frame(temp_df)) + 
geom_boxplot(mapping = aes(x = variable,y = value))  + 
ylab("inner product") + xlab("")+ theme(text = element_text(size = 24))

########## Vary Scale #########

prior_l <- .1
nu_param <- .0001

params <- c(.25,.5,1,2)

legend_labs <- rep("",length(params))
for(ind in 1:length(legend_labs)){
  legend_labs[ind] <- paste0(params[ind])
}

# Look at some draws from the prior
num_draws <- 100
ortho_draws <- array(0,dim = c(M_grid,num_draws,length(params)))
for(ind in 1:length(params)){
  prior_scale <- params[ind]
  prior_cov <- make_cov(time_grid,prior_l,prior_scale) 
  
  ortho_cov <- prior_cov - 
    prior_cov%*%W_grid%*%lambda_minus_k%*%
    inv(nu_param + t(lambda_minus_k)%*%W_grid%*%prior_cov%*%W_grid%*%lambda_minus_k)%*%
    t(lambda_minus_k)%*%W_grid%*%prior_cov
  
  ortho_draws[,,ind] <- t(mvrnorm(num_draws, mu = rep(0,length(time_grid)),Sigma = ortho_cov))
}

yl <- c(min(ortho_draws),max(ortho_draws))
count <- 1
for(ind in 1:length(params)){
  temp_df <- as.data.frame(cbind(time_grid,ortho_draws[,,ind]))
  temp_df <- melt(temp_df,id = "time_grid")
  print(
    ggplot() + 
      geom_line(data = temp_df, aes(x = time_grid, y = value, group = variable),alpha = .1,size = 2) +  
      ylab(expression(lambda[2](t))) + xlab("t")+theme(text = element_text(size = 24))  + ylim(yl)
  )
  count <- count + 1
}

# what happens to inner product
num_draws <- 1000
inner_prod <- array(0,dim = c(num_draws,length(params)))
ortho_draws <- array(0,dim = c(M_grid,num_draws,length(params)))
for(ind in 1:length(params)){
  prior_scale <- params[ind]
  prior_cov <- make_cov(time_grid,prior_l,prior_scale) 
  
  ortho_cov <- prior_cov - 
    prior_cov%*%W_grid%*%lambda_minus_k%*%
    inv(nu_param + t(lambda_minus_k)%*%W_grid%*%prior_cov%*%W_grid%*%lambda_minus_k)%*%
    t(lambda_minus_k)%*%W_grid%*%prior_cov
  
  ortho_draws[,,ind] <- t(mvrnorm(num_draws, mu = rep(0,length(time_grid)),Sigma = ortho_cov))
  
  inner_prod[,ind] <- t(ortho_draws[,,ind])%*%W_grid%*%lambda_minus_k
}

temp_df <- as.data.frame(inner_prod)
colnames(temp_df) <- legend_labs
temp_df <- melt(temp_df)
ggplot(data = as.data.frame(temp_df)) + 
  geom_boxplot(mapping = aes(x = variable,y = value)) + 
  ylab("inner product") + xlab("")+ theme(text = element_text(size = 24))

########## Vary Nu #########

prior_l <- .1
prior_scale <- 1

params <- c(.0001,.01,1,100)

legend_labs <- rep("",length(params))
for(ind in 1:length(legend_labs)){
  legend_labs[ind] <- paste0(params[ind])
}

# Look at some draws from the prior
num_draws <- 100
ortho_draws <- array(0,dim = c(M_grid,num_draws,length(params)))
for(ind in 1:length(params)){
  prior_cov <- make_cov(time_grid,prior_l,prior_scale)
  
  nu_param <- params[ind]
  ortho_cov <- prior_cov - 
    prior_cov%*%W_grid%*%lambda_minus_k%*%
    inv(nu_param + t(lambda_minus_k)%*%W_grid%*%prior_cov%*%W_grid%*%lambda_minus_k)%*%
    t(lambda_minus_k)%*%W_grid%*%prior_cov
  
  ortho_draws[,,ind] <- t(mvrnorm(num_draws, mu = rep(0,length(time_grid)),Sigma = ortho_cov))
}

yl <- c(min(ortho_draws),max(ortho_draws))
count <- 1
for(ind in 1:length(params)){
  temp_df <- as.data.frame(cbind(time_grid,ortho_draws[,,ind]))
  temp_df <- melt(temp_df,id = "time_grid")
  print(
    ggplot() + 
      geom_line(data = temp_df, aes(x = time_grid, y = value, group = variable),alpha = .1,size = 2) +  
      ylab(expression(lambda[2](t))) + xlab("t")+theme(text = element_text(size = 24)) + ylim(yl)
  )
  count <- count + 1
}

# what happens to inner product
num_draws <- 1000
inner_prod <- array(0,dim = c(num_draws,length(params)))
ortho_draws <- array(0,dim = c(M_grid,num_draws,length(params)))
for(ind in 1:length(params)){
  
  prior_cov <- make_cov(time_grid,prior_l,prior_scale) 
  
  nu_param <- params[ind]
  ortho_cov <- prior_cov - 
    prior_cov%*%W_grid%*%lambda_minus_k%*%
    inv(nu_param + t(lambda_minus_k)%*%W_grid%*%prior_cov%*%W_grid%*%lambda_minus_k)%*%
    t(lambda_minus_k)%*%W_grid%*%prior_cov
  
  ortho_draws[,,ind] <- t(mvrnorm(num_draws, mu = rep(0,length(time_grid)),Sigma = ortho_cov))
  
  inner_prod[,ind] <- t(ortho_draws[,,ind])%*%W_grid%*%lambda_minus_k
}

temp_df <- as.data.frame(inner_prod)
colnames(temp_df) <- legend_labs
temp_df <- melt(temp_df)
ggplot(data = as.data.frame(temp_df)) + 
  geom_boxplot(mapping = aes(x = variable,y = value)) + 
  ylab("inner product")+ xlab("") + theme(text = element_text(size = 24))

