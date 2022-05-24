library(tidyverse)
library(lubridate)
library(reshape2)
library(haven)
library(Rcpp)
library(ggplot2)
library(splines)
library(boot)
library(fda)
library(fdapace)
library(SLFPCA)

load("./results/cebu_wb_mcmc.RData")
sourceCpp("../src/remo_fpca.cpp")
hd_figs <- "./figures/"
im_save = TRUE

# mean process
if(im_save){png(file=paste0(hd_figs,'wb_mu_y.png'), width=480, height=480)}
print(
  ggplot() + 
    geom_line(aes(x = age_grid, y = y_scale*Y_mean),size = 2) + 
    ylab("weight (grams)") + xlab("age (years)")+theme(text = element_text(size = 24)) 
  #+ geom_line(aes(x = time,y = lambda_mean[,k]),size = 2)
)
if(im_save){dev.off()}

# FPCs
fpc_temp <- y_scale*E*(rep(1,M)%*%t(sqrt(eigen(sm_Gt)$values[1:5])))
lambda_df <- melt(as.data.frame(cbind(age_grid,fpc_temp)),id = "age_grid")
lambda_df$variable<-as.factor(as.numeric(lambda_df$variable))
colnames(lambda_df)[2]<- "FPC"
if(im_save){png(file=paste0(hd_figs,'wb_lambda_y.png'), width=480, height=480)}
print(
  ggplot(lambda_df) + 
    geom_line(aes(x = age_grid, y = value, color = FPC),size = 2) + 
    ylab("weight (grams)") + xlab("age (years)")+theme(text = element_text(size = 24)) 
)
if(im_save){dev.off()}

library(GGally)
# FPC scores
freq_eta <- apply(wb_output$sims.list$xi,c(2,3),mean)
freq_eta_df <- as.data.frame(freq_eta)
colnames(freq_eta_df) <- 1:K
if(im_save){png(file=paste0(hd_figs,'wb_xi.png'), width=480, height=480)}
ggpairs(freq_eta_df,axisLabels = "internal")+ ylab(expression(xi)) + xlab(expression(xi))+theme(text = element_text(size = 24))
if(im_save){dev.off()}


# influence of covariates on fit
cov_titles <- c("mother's height","sex","strata","season of birth")
fit_all <- array(0,dim = c(M,4,n_iter))
for(cov_ind in 1:4){
  x_new <- as.numeric(1:4 == cov_ind)
  for(iter in 1:n_iter){
    fit_all[,cov_ind,iter] <- E%*%(wb_output$sims.list$beta_scalar[iter,,]%*%x_new)*y_scale
  }
}
thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
yl <- c(min(fit_all[,,thin_iter]),max(fit_all[,,thin_iter]))

for(cov_ind in 1:4){
  x_new <- as.numeric(1:4 == cov_ind)
  
  fit_save <- array(0,dim = c(M,n_iter))
  for(iter in 1:n_iter){
    fit_save[,iter] <-  E%*%(wb_output$sims.list$beta_scalar[iter,,]%*%x_new)*y_scale
  }
  
  thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
  fit_thin_df <- as.data.frame(cbind(age_grid,fit_save[,thin_iter]))
  fit_thin_df <- melt(fit_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'wb_fit_b_',cov_ind,'.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) + 
      ylim(yl) + ylab("weight (grams)") + xlab("age (years)")+ ggtitle(cov_titles[cov_ind]) + 
      theme(text = element_text(size = 24)) 
  )
  if(im_save){dev.off()}
}


# plot some overall fits

n_iter <- 1000
fit_save <- array(0,dim = c(M,N,n_iter))
Y_mean <- colMeans(t(y),na.rm=TRUE)

for(iter in 1:n_iter){
  # fpca regression
  fit_save[,,iter] <- Y_mean%*%t(rep(1,N)) + E%*%wb_output$sims.list$beta_scalar[iter,,]%*%x + 
    E%*%t(wb_output$sims.list$xi[iter,,])
  # historical regression
  fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_fun_reg_fit(historical_integrals,wb_output$sims.list$beta_bf[iter,]))
  # concurrent regression 
  fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_concurrent_fit(z_ill_complete,beta_ill_basis,wb_output$sims.list$beta_ill[iter,,]))
}

all_id <- NULL
all_age <- NULL
all_weight <- NULL

for(i in 1:N){
  temp_ind <- which(obs_grid_y[,i] == 1)
  all_id <- c(all_id,rep(i,length(temp_ind)))
  all_age <- c(all_age,age_grid[temp_ind])
  all_weight <- c(all_weight,y_common[temp_ind,i])
}
cebu_all = data.frame(id = as.factor(all_id), age = all_age, weight = y_scale*all_weight)

# What subjects do you want to look at?
# extreme magnitude - 1905
# extreme shape - 2601, 1926, 2205
# missing observations - 1626, 404, 1333
# typical observaions - 914,421
vis_subset <- c(404,421,914,1333,1626,1905,1926,2205,2601)
for(obs_ind in vis_subset){
  thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
  fit_thin_df <- as.data.frame(cbind(age_grid,y_scale*fit_save[,obs_ind,thin_iter]))
  fit_thin_df <- melt(fit_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'wb_fit_',obs_ind,'.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .025,size = 2) +  
      ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24)) +
      #ggtitle(obs_ind)+
      geom_point(aes(x = cebu_all[cebu_all$id == obs_ind,2],y = cebu_all[cebu_all$id == obs_ind,3]), size = 5) 
  )
  if(im_save){dev.off()}
}


# breastfeeding profiles
ones_mat <- matrix(1,nrow = M,ncol = M)
ones_mat[lower.tri(ones_mat)] <- 0
bf_profiles <-  ones_mat
n_profiles <- dim(bf_profiles)[2]
profile_integrals <- compute_historical_integrals(bf_profiles,basis_funs,grid_weight)
profile_fit <- array(0,dim = c(M-1,n_profiles,n_iter))
for(iter in 1:n_iter){
  profile_fit[,,iter] <- compute_fun_reg_fit(profile_integrals,wb_output$sims.list$beta_bf[iter,])
}
yl <- c(min(y_scale*profile_fit),max(y_scale*profile_fit))

# plot cross sections of covariance surface
for(age_ind in 1:n_profiles){
  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
  profile_thin_df <- as.data.frame(cbind(age_grid[2:M],y_scale*profile_fit[,age_ind,thin_iter]))
  colnames(profile_thin_df)[1] <- "age_grid"
  profile_thin_df <- melt(profile_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'wb_bf_profile_',age_ind,'_fit.png'), width=480, height=480)}
  #yl <- c(min(Y_atgpxr[,3]),max(Y_atgpxr[,3]))
  print(
    ggplot() + 
      geom_line(data = profile_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
      ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24)) + ylim(yl)
  )
  if(im_save){dev.off()}
}

# how much do the profiles change from zero?
profile_norm <- array(0,dim = c(n_profiles,n_iter))
for(iter in 1:n_iter){
  profile_norm[,iter] <- y_scale*sqrt(diag(t(profile_fit[,,iter])%*%grid_weight[2:M,2:M]%*%profile_fit[,,iter]))
}
profile_norm_df <- melt(profile_norm[,(n_iter/2):n_iter])
yl <- c(0,max(profile_norm_df[,3]))
if(im_save){png(file=paste0(hd_figs,'wb_bf_profile_norm_fit.png'), width=480, height=480)}
ggplot() + geom_boxplot(aes(group = as.factor(age_grid[profile_norm_df$Var1]),y = profile_norm_df$value)) +
  scale_x_continuous(breaks=c(-.3455,-.3455/2,0,.3455/2,.3455),
                     labels=c("0",".5", "1","1.5", "2")) +
  ylab("norm") + xlab("age (years)")+ ylim(yl)+ theme(text = element_text(size = 24)) 
if(im_save){dev.off()}


# illness inference
yl <- c(-150,50)
cov_titles <- c("diarrhea","fever","cough")
for(cov_ind in 1:num_covs){
  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
  mu_thin_df <- as.data.frame(cbind(age_grid[2:M],y_scale*beta_ill_basis%*%t(wb_output$sims.list$beta_ill[thin_iter,,cov_ind])))
  colnames(mu_thin_df)[1] <- "age_grid"
  mu_thin_df <- melt(mu_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'wb_beta_ill_',cov_ind,'_samples.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = mu_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ggtitle(cov_titles[cov_ind]) + 
      #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
      ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24)) + ylim(yl)
  )
  if(im_save){dev.off()}
}

##### compute WAIC #####

n_iter <- 1000
fit_save <- array(0,dim = c(M,N,n_iter))
Y_mean <- colMeans(t(y),na.rm=TRUE)

for(iter in 1:n_iter){
  # fpca regression
  fit_save[,,iter] <- Y_mean%*%t(rep(1,N)) + E%*%wb_output$sims.list$beta_scalar[iter,,]%*%x + 
    E%*%t(wb_output$sims.list$xi[iter,,])
  # historical regression
  fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_fun_reg_fit(historical_integrals,wb_output$sims.list$beta_bf[iter,]))
  # concurrent regression 
  fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_concurrent_fit(z_ill_complete,beta_ill_basis,wb_output$sims.list$beta_ill[iter,,]))
}

thin_iter <- round(seq(n_iter/2 +1,n_iter,length.out = 500))
waic_vec_wb <- NULL
lppd <- 0
p_waic <- 0
for(i in 1:N){
  temp_ind <- which(obs_grid_y[,i] == 1)
  lik_temp <- array(0,dim = c(length(temp_ind),length(thin_iter)))
  for(m in 1:length(temp_ind)){
    lik_temp[m,] <- dnorm(y_common[temp_ind[m],i],
                          mean = fit_save[temp_ind[m],i,thin_iter],
                          sd = sqrt(1/wb_output$sims.list$taueps[thin_iter]))
  }
  lppd <- lppd + sum(log(apply(lik_temp,1,mean)))
  
  p_waic <- p_waic + sum(apply(log(lik_temp),1,var))
  waic_vec_wb <- c(waic_vec_wb,log(apply(lik_temp,1,mean)) - apply(log(lik_temp),1,var))
}
waic_wb <- -2*(lppd - p_waic)


