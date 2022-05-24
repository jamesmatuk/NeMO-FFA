library(tidyverse)
library(lubridate)
library(reshape2)
library(haven)
library(Rcpp)
library(ggplot2)
library(splines)

##### Import data #####
sourceCpp("../src/remo_fpca.cpp")
cebu_long <- read_dta("data/mlong.dta")
cebu_mort <- read_dta("data/mmort.dta")
cebu_birth <- read_dta("data/mbirth2.dta")
  
num_subjects <-dim(unique( cebu_long[,c("basebrgy","basewman")]))[1]

##### subset relevant data #####
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
K <- 5
K_bf <- 3
Q <- dim(x)[1]
one_vec_N <- rep(1,N)
nu_eta <- 1e-4
nu_param <- 1e-4

PX <- FALSE
px_b <- 10
px_a <- px_b + 1

prior_a <- 15
prior_b <- 5
prior_var <- 1

# initilize ortho gp 
lambda_y_cur <- matrix(0,nrow = M,ncol = K)
xi_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
b_cur <- matrix(rnorm(K*Q,0,1),nrow = K,ncol = Q)
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
num_centers <- 4
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

##### MCMC preprocessing #####

im_save <- TRUE
hd_figs <- "./figures/"
mcmc_output <- readRDS("./results/cebu_mcmc.rds")

mu_save <- mcmc_output$mu_y_save 
scale_mu_save <- mcmc_output$scale_mu_y_save
l_mu_save <- mcmc_output$l_mu_y_save 
lambda_save <- mcmc_output$lambda_y_save
xi_save <- mcmc_output$xi_save
b_save <- mcmc_output$b_save
eta_save <- mcmc_output$eta_y_save 
psi_save <- mcmc_output$psi_y_save 
l_param_save <- mcmc_output$l_y_param_save
scale_param_save <- mcmc_output$scale_y_param_save
sig_sq_save <- mcmc_output$sig_sq_save 
n_iter <- length(sig_sq_save)


# Does there appear to be any label/sign switching
thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
yl <- c(y_scale*min(lambda_save[,,thin_iter]),y_scale*max(lambda_save[,,thin_iter]))
for(k in 1:dim(lambda_save)[2]){
  lambda_thin_df <- as.data.frame(cbind(age_grid,y_scale*lambda_save[,k,thin_iter]))
  lambda_thin_df <- melt(lambda_thin_df,id = "age_grid")
  print(
    ggplot() + 
      geom_line(data = lambda_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ylab("") + xlab("")+theme(text = element_text(size = 24))  + ylim(yl)
    #+ geom_line(aes(x = time,y = lambda_mean[,k]),size = 2)
  )
}
# it looks like there is some switching between the 1st and 2nd lambda

# post - process samples
lambda_post <- lambda_save
xi_post <- xi_save
b_post <- b_save
eta_post <- eta_save
M <- dim(lambda_save[,,1])[1]
K <- dim(lambda_save[,,1])[2]
N <- dim(eta_save[,,1])[2]
one_M <- rep(1,M)
one_N <- rep(1,N)
for(iter in 1:n_iter){
  lambda_post[,,iter] <- lambda_save[,,iter]*(one_M%*%t(sqrt(psi_save[,iter])))
  eta_post[,,iter] <- eta_save[,,iter]/(sqrt(psi_save[,iter])%*%t(one_N))
  b_post[,,iter] <- b_save[,,iter]/(sqrt(psi_save[,iter])%*%t(rep(1,Q)))
  xi_post[,,iter] <- xi_save[,,iter]/(sqrt(psi_save[,iter])%*%t(one_N))
}

# resolve label switching an sign ambiguity
sourceCpp("../src/msf.cpp")
post_burn_iter <- 1:n_iter
lambda_save_processed <- array(0,dim = c(M,K,length(post_burn_iter)))
eta_save_processed <- array(0,dim = c(K,N,length(post_burn_iter)))
xi_save_processed <- array(0,dim = dim(xi_save))
b_save_processed <- array(0,dim = dim(b_save))
l_param_save_processed <- array(0,dim = c(K,length(post_burn_iter)))
scale_param_save_processed <- array(0,dim = c(K,length(post_burn_iter)))
for(iter in 1:length(post_burn_iter)){
  perm_out <- msfOUT_functional(lambda_post[,,post_burn_iter[iter]],W,lambda_post[,,round(n_iter/(1.1))])
  #perm_out <- 1:5
  lambda_save_processed[,,iter] <- aplr(lambda_post[,,post_burn_iter[iter]],perm_out)
  eta_save_processed[,,iter] <- t(aplr(t(eta_post[,,post_burn_iter[iter]]),perm_out))
  xi_save_processed[,,iter] <- t(aplr(t(xi_post[,,post_burn_iter[iter]]),perm_out))
  b_save_processed[,,iter] <- t(aplr(t(b_post[,,post_burn_iter[iter]]),perm_out))
  l_param_save_processed[,iter] <- l_param_save[abs(perm_out),post_burn_iter[iter]]
  scale_param_save_processed[,iter] <- scale_param_save[abs(perm_out),post_burn_iter[iter]]
}

# arrange lambda related parameters by magnitude of lambda
lambda_post_mean <- apply(lambda_save_processed,c(1,2),mean)
lambda_mean_norm <- sqrt(diag(t(lambda_post_mean)%*%W%*%lambda_post_mean))
lambda_mean_order <- order(lambda_mean_norm,decreasing = T)

lambda_save_processed <- lambda_save_processed[,lambda_mean_order,]
eta_save_processed <- eta_save_processed[lambda_mean_order,,]
xi_save_processed <- xi_save_processed[lambda_mean_order,,]
b_save_processed <- b_save_processed[lambda_mean_order,,]
l_param_save_processed <- l_param_save_processed[lambda_mean_order,]
scale_param_save_processed <- scale_param_save_processed[lambda_mean_order,]

eta_post_mean <- apply(eta_save_processed[,,(n_iter/2):n_iter],c(1,2),mean)
extreme_obs <- order(diag(t(eta_post_mean)%*%eta_post_mean),decreasing = T)[1:10]

##### trace plots for regression #####
# error variance
if(im_save){png(file=paste0(hd_figs,'trace_error_variance.png'), width=480, height=480)}
qplot((.1*n_iter):n_iter,y_scale^2*sig_sq_save[(.1*n_iter):n_iter],geom = "path") + xlab("iteration") + ylab("value") + 
  ggtitle("") + theme(text = element_text(size = 28)) 
if(im_save){dev.off()}

# mu covariance scale 
if(im_save){png(file=paste0(hd_figs,'trace_mu_y_cov_scale.png'), width=480, height=480)}
qplot(post_burn_iter,y_scale*sqrt(scale_mu_save[post_burn_iter]),geom = "path") + xlab("iteration") + ylab("value") + 
  ggtitle("") +theme(text = element_text(size = 28)) 
if(im_save){dev.off()}

# mu covariance length
if(im_save){png(file=paste0(hd_figs,'trace_mu_y_cov_length.png'), width=480, height=480)}
qplot(post_burn_iter,l_mu_save[post_burn_iter],geom = "path") + xlab("iteration") + ylab("value") + 
  ggtitle("") +theme(text = element_text(size = 28))
if(im_save){dev.off()}

# length parameter for loadings
for(k in 1:K){
  if(im_save){png(file=paste0(hd_figs,'trace_lambda_y_cov_l_',k,'.png'), width=480, height=480)}
  print(
    qplot(post_burn_iter,l_param_save_processed[k,post_burn_iter],geom = "path") + 
      ggtitle("") + xlab("iteration") + ylab("value") +  theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
}

# scale parameter for loadings
for(k in 1:K){
  if(im_save){png(file=paste0(hd_figs,'trace_lambda_y_cov_scale_',k,'.png'), width=480, height=480)}
  print(
    qplot(post_burn_iter[(.1*n_iter):n_iter],y_scale*sqrt(scale_param_save_processed[k,post_burn_iter[(.1*n_iter):n_iter]]),geom = "path") + 
      ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
}

# one observation of eta parameter
obs_ind <- 421
for(k in 1:K){
  if(im_save){png(file=paste0(hd_figs,'trace_xi_',k,'_obs_',obs_ind,'.png'), width=480, height=480)}
  print(
    qplot(post_burn_iter,xi_save_processed[k,obs_ind,post_burn_iter],geom = "path") + 
      ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
}


# regression trace plots
Q <- dim(b_save_processed)[2]
for(k in 1:K){
  for(q in 1:Q){
    if(im_save){png(file=paste0(hd_figs,'trace_b_k_',k,'_q_',q,'.png'), width=480, height=480)}
    print(
      qplot(post_burn_iter,b_save_processed[k,q,post_burn_iter],geom = "path") + # ggtitle(paste0("k = ",k,", q = ",q)) +
        ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28)) 
      
    )
    if(im_save){dev.off()}
  }
}

# fit of regression coefficients
fit_all <- array(0,dim = c(M,4,n_iter))
for(cov_ind in 1:4){
  x_new <- as.numeric(1:4 == cov_ind)
  for(iter in 1:n_iter){
    fit_all[,cov_ind,iter] <- (lambda_save[,,iter]%*%b_save[,,iter]%*%x_new)*y_scale
  }
  
}

for(cov_ind in 1:4){
if(im_save){png(file=paste0(hd_figs,'trace_fit_b_q_',cov_ind,'_t1.png'), width=480, height=480)}
print(
qplot(post_burn_iter,fit_all[1,cov_ind,],geom = "path") + 
  ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28))
)
if(im_save){dev.off()}
if(im_save){png(file=paste0(hd_figs,'trace_fit_b_q_',cov_ind,'_t6.png'), width=480, height=480)}
print(
qplot(post_burn_iter,fit_all[6,cov_ind,],geom = "path") + 
  ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28)) 
)
if(im_save){dev.off()}
if(im_save){png(file=paste0(hd_figs,'trace_fit_b_q_',cov_ind,'_t13.png'), width=480, height=480)}
print(
qplot(post_burn_iter,fit_all[13,cov_ind,],geom = "path") + 
  ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28)) 
)
if(im_save){dev.off()}
}


##### inference plots for regression #####

# plot the actual data
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
if(im_save){png(file=paste0(hd_figs,'cebu_data.png'), width=480, height=480)}
ggplot(cebu_all) + geom_point(aes(x = age,y = weight),alpha = .5) + 
  ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24))
if(im_save){dev.off()}


# pattern of missingness for weight
prop_na <- apply(obs_grid_y == 0,1,mean)
temp_grid <- NULL
is_missing <- NULL
for(m in 1:length(prop_na)){
  temp_grid <- c(temp_grid,rep(age_grid[m],100))
  prop_na_m <- round(prop_na[m]*100)
  is_missing <- c(is_missing,c(rep("missing",prop_na_m),rep("recorded",100-prop_na_m)))
}
if(im_save){png(file=paste0(hd_figs,'proportion_weight_missing.png'), width=480, height=480)}
na_df <- data.frame(age_grid = temp_grid,is_missing = is_missing)
ggplot(na_df) + geom_bar(aes(x = age_grid,fill = is_missing),position="fill") + 
  ylab("missingness") + xlab("age (years)") + theme(text = element_text(size = 24),legend.title = element_blank())
if(im_save){dev.off()}

# plot thinned posterior samples of mu
thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
mu_thin_df <- as.data.frame(cbind(age_grid,y_scale*mu_save[,thin_iter]))
mu_thin_df <- melt(mu_thin_df,id = "age_grid")
if(im_save){png(file=paste0(hd_figs,'mu_y_samples.png'), width=480, height=480)}
#yl <- c(min(Y_atgpxr[,3]),max(Y_atgpxr[,3]))
print(
  ggplot() + 
    geom_line(data = mu_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
    #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
    ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24)) #+ ylim(yl)
)
if(im_save){dev.off()}

# plot thinned posterior samples of lambda
thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
yl <- c(y_scale*min(lambda_save_processed[,,thin_iter]),y_scale*max(lambda_save_processed[,,thin_iter]))
for(k in 1:K){
  lambda_thin_df <- as.data.frame(cbind(age_grid,y_scale*lambda_save_processed[,k,thin_iter]))
  lambda_thin_df <- melt(lambda_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'lambda_y_',k,'_samples.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = lambda_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ylab("weight (grams)") + xlab("age (years)")+theme(text = element_text(size = 24))  + ylim(yl)
    #+ geom_line(aes(x = time,y = lambda_mean[,k]),size = 2)
  )
  if(im_save){dev.off()}

}


# posterior means for FPC scores
library(GGally)
xi_post_mean <- apply(xi_save_processed,c(1,2),mean)
xi_df <- as.data.frame(t(xi_post_mean))
colnames(xi_df) <- 1:length(colnames(xi_df))
if(im_save){png(file=paste0(hd_figs,'xi_posterior_means.png'), width=480, height=480)}
ggpairs(xi_df,axisLabels = "internal")+ ylab(expression(xi)) + xlab(expression(xi))+theme(text = element_text(size = 24))
if(im_save){dev.off()}

# inference for regression coefficients
# boxplots of posterior samples
b_df <- melt(b_save_processed[,,(n_iter/2):n_iter])
colnames(b_df)[1:3] <- c("k","q","mcmc_iter")
covariate_names = c("mother's height","sex","strata","season of birth")
for(cov in 1:length(covariate_names)){
  xi_df_subset <- b_df[b_df$q == cov,]
  if(im_save){png(file=paste0(hd_figs,'boxplot_b_',cov,'.png'), width=480, height=480)}
  print(
    ggplot(xi_df_subset) + geom_boxplot(aes(x = as.factor(k),y = value,fill=as.factor(q)))+ 
      facet_wrap(~covariate_names[q])+ geom_hline(yintercept = 0,color = "black",linetype = "dashed",size = 2) + 
      ggtitle("") + xlab("k") + ylab("value") + theme(text = element_text(size = 24),legend.position = "none") 
  )
  if(im_save){dev.off()}
}

# influence of covariates on fit
cov_titles <- c("mother's height","sex","strata","season of birth")
fit_all <- array(0,dim = c(M,4,n_iter))
for(cov_ind in 1:4){
  x_new <- as.numeric(1:4 == cov_ind)
  for(iter in 1:n_iter){
    fit_all[,cov_ind,iter] <- (lambda_save[,,iter]%*%b_save[,,iter]%*%x_new)*y_scale
  }
}
thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
yl <- c(min(fit_all[,,thin_iter]),max(fit_all[,,thin_iter]))

for(cov_ind in 1:4){
  x_new <- as.numeric(1:4 == cov_ind)
  
  fit_save <- array(0,dim = c(M,n_iter))
  for(iter in 1:n_iter){
    fit_save[,iter] <- lambda_save[,,iter]%*%b_save[,,iter]%*%x_new
  }
  
  thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
  fit_thin_df <- as.data.frame(cbind(age_grid,y_scale*fit_save[,thin_iter]))
  fit_thin_df <- melt(fit_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'fit_b_',cov_ind,'.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) + 
      ylim(yl) + ylab("weight (grams)") + xlab("age (years)")+ ggtitle(cov_titles[cov_ind]) + 
      theme(text = element_text(size = 24)) 
  )
  if(im_save){dev.off()}
}

##### binary FPCA for breastfeeding #####

mu_save <- mcmc_output$mu_bf_save 
scale_mu_save <- mcmc_output$scale_mu_bf_save
l_mu_save <- mcmc_output$l_mu_bf_save 
lambda_save <- mcmc_output$lambda_bf_save
eta_save <- mcmc_output$eta_bf_save 
psi_save <- mcmc_output$psi_bf_save 
l_param_save <- mcmc_output$l_bf_param_save
scale_param_save <- mcmc_output$scale_bf_param_save

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
  perm_out <- msfOUT_functional(lambda_post[,,post_burn_iter[iter]],W,lambda_post[,,n_iter/(1)])
  lambda_save_processed[,,iter] <- aplr(lambda_post[,,post_burn_iter[iter]],perm_out)
  eta_save_processed[,,iter] <- t(aplr(t(eta_post[,,post_burn_iter[iter]]),perm_out))
  l_param_save_processed[,iter] <- l_param_save[abs(perm_out),post_burn_iter[iter]]
  scale_param_save_processed[,iter] <- scale_param_save[abs(perm_out),post_burn_iter[iter]]
}


lambda_post_mean <- apply(lambda_save_processed,c(1,2),mean)
lambda_mean_norm <- sqrt(diag(t(lambda_post_mean)%*%W%*%lambda_post_mean))
lambda_order <- order(lambda_mean_norm,decreasing = T)

lambda_save_processed <- lambda_save_processed[,lambda_order,] 
eta_save_processed <- eta_save_processed[lambda_order,,] 
l_param_save_processed <- l_param_save_processed[lambda_order,] 
scale_param_save_processed <- scale_param_save_processed[lambda_order,] 

##### Interpretting parameters on scale of [0,1] #####

# pattern of missingness for breastfeeding
prop_na <- apply(obs_grid_bf == 0,1,mean)
temp_grid <- NULL
is_missing <- NULL
for(m in 1:length(prop_na)){
  temp_grid <- c(temp_grid,rep(age_grid[m],100))
  prop_na_m <- round(prop_na[m]*100)
  is_missing <- c(is_missing,c(rep("missing",prop_na_m),rep("recorded",100-prop_na_m)))
}
if(im_save){png(file=paste0(hd_figs,'proportion_bf_missing.png'), width=480, height=480)}
na_df <- data.frame(age_grid = temp_grid,is_missing = is_missing)
ggplot(na_df) + geom_bar(aes(x = age_grid,fill = is_missing),position="fill") +
  ggtitle("breastfeeding") + 
  ylab("missingness") + xlab("age (years)") + theme(text = element_text(size = 24),legend.title = element_blank())
if(im_save){dev.off()}

# plot thinned posterior samples of mu
thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
mu_thin_df <- as.data.frame(cbind(age_grid,pnorm(mu_save[,thin_iter])))
mu_thin_df <- melt(mu_thin_df,id = "age_grid")
if(im_save){png(file=paste0(hd_figs,'mu_bf_probit_samples.png'), width=480, height=480)}
#yl <- c(min(Y_atgpxr[,3]),max(Y_atgpxr[,3]))
yl <- c(0,1)
print(
  ggplot() + 
    geom_line(data = mu_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
    #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
    ylab("probability") + xlab("age (years)")+ theme(text = element_text(size = 24)) + ylim(yl)
)
if(im_save){dev.off()}


# interpretting lambdas
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_colors <- gg_color_hue(2)

thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
yl <- c(0,1)
mu_thin_df <- as.data.frame(cbind(age_grid,mu_save[,thin_iter]))
mu_thin_df <- melt(mu_thin_df,id = "age_grid")
for(k in 1:K){
  mean_perterb <- array(0,dim = c(M,length(thin_iter)))
  for(iter in 1:length(thin_iter)){
    mean_perterb[,iter] <- lambda_save_processed[,k,thin_iter[iter]]*sd(eta_save_processed[k,,thin_iter[iter]])
  }
  mean_perterb_df <- as.data.frame(cbind(age_grid,mean_perterb))
  mean_perterb_df <- melt(mean_perterb_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'lambda_bf_',k,'_probit_samples.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(aes(x = mu_thin_df$age_grid, y = pnorm(mu_thin_df$value), group = mu_thin_df$variable),alpha = .05,size = 2) + 
      geom_line(aes(x = mu_thin_df$age_grid, y = pnorm(mu_thin_df$value + mean_perterb_df$value), group = mu_thin_df$variable),alpha = .05,size = 2,color = gg_colors[1]) + 
      geom_line(aes(x = mu_thin_df$age_grid, y = pnorm(mu_thin_df$value - mean_perterb_df$value), group = mu_thin_df$variable),alpha = .05,size = 2,color = gg_colors[2]) + 
      ylab("probability") + xlab("age (years)")+theme(text = element_text(size = 24))  + ylim(yl)
    #+ geom_line(aes(x = time,y = lambda_mean[,k]),size = 2)
  )
  if(im_save){dev.off()}
}

# look at some fits 

# plot fitted observations
fit_save <- array(0,dim = c(M,N,n_iter))
for(s in 1:N){
  for(iter in 1:n_iter){
    fit_save[,s,iter] <- mu_save[,iter] + lambda_save[,,iter]%*%eta_save[,s,iter]
  }
}

# all ones: 412
# all zeros: 2162
# changed from bf to non bf: 1188,1969
# some missing: 2575,763
#vis_subset <- c(412,2162,1188,1969,2575,763)

vis_subset <- c(350,1079,2019,2087)
#vis_subset <- sample(which(apply(is.na(z_bf),2,sum)>0),20)
for(obs_ind in vis_subset){
  thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
  fit_thin_df <- as.data.frame(cbind(age_grid,pnorm(fit_save[,obs_ind,thin_iter])))
  fit_thin_df <- melt(fit_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'fit_bf_',obs_ind,'.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ylim(c(0,1)) + ylab("probability") + xlab("age (years)")+ theme(text = element_text(size = 24)) +
      ggtitle("")+
      geom_point(aes(x = age_grid,y = as.numeric(z_bf[,obs_ind])),size = 10) 
  )
  if(im_save){dev.off()}
}

##### Trace plots #####

# scale parameter for loadings
if(im_save){png(file=paste0(hd_figs,'trace_fit_bf_t1.png'), width=480, height=480)}
  qplot(post_burn_iter,pnorm(fit_save[1,2087,]),geom = "path") + 
    ggtitle("") + xlab("iteration") + ylab("probability") + theme(text = element_text(size = 28)) + ylim(c(0,1))
if(im_save){dev.off()}
  if(im_save){png(file=paste0(hd_figs,'trace_fit_bf_t6.png'), width=480, height=480)}
  qplot(post_burn_iter,pnorm(fit_save[6,2087,]),geom = "path") + 
    ggtitle("") + xlab("iteration") + ylab("probability") + theme(text = element_text(size = 28)) + ylim(c(0,1))
if(im_save){dev.off()}
if(im_save){png(file=paste0(hd_figs,'trace_fit_bf_t13.png'), width=480, height=480)}
qplot(post_burn_iter,pnorm(fit_save[13,2087,]),geom = "path") + 
    ggtitle("") + xlab("iteration") + ylab("probability") + theme(text = element_text(size = 28)) + ylim(c(0,1))
if(im_save){dev.off()}

# Look at all eta posterior means 
eta_post_mean <- apply(eta_save_processed[,,(n_iter/2):n_iter],c(1,2),mean)
library(GGally)
eta_df <- as.data.frame(t(eta_post_mean))
colnames(eta_df) <- 1:length(colnames(eta_df))
if(im_save){png(file=paste0(hd_figs,'eta_bf_posterior_means.png'), width=480, height=480)}
ggpairs(eta_df,axisLabels = "internal")+ ylab(expression(eta)) + xlab(expression(eta))+theme(text = element_text(size = 24))
if(im_save){dev.off()}

# what is the crazy structure in the etas? 

vis_subset <-  order(abs(eta_post_mean[3,]),decreasing = T)[1:3]
#vis_subset <- sample(N,10,replace = F)
for(obs_ind in vis_subset){
  thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
  fit_thin_df <- as.data.frame(cbind(age_grid,pnorm(fit_save[,obs_ind,thin_iter])))
  fit_thin_df <- melt(fit_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'fit_bf_',obs_ind,'_eta_3_outliers.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ylim(c(0,1)) + ylab("probability") + xlab("age (years)")+ theme(text = element_text(size = 24)) +
      #ggtitle(obs_ind)+
      geom_point(aes(x = age_grid,y = as.numeric(z_bf[,obs_ind])),size = 10) 
  )
  if(im_save){dev.off()}
}

vis_subset <-  order(abs(eta_post_mean[2,]),decreasing = T)[1:3]
#vis_subset <- sample(N,10,replace = F)
for(obs_ind in vis_subset){
  thin_iter <- seq(n_iter/2,n_iter,length.out = 100)
  fit_thin_df <- as.data.frame(cbind(age_grid,pnorm(fit_save[,obs_ind,thin_iter])))
  fit_thin_df <- melt(fit_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'fit_bf_',obs_ind,'_eta_2_outliers.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ylim(c(0,1)) + ylab("probability") + xlab("age (years)")+ theme(text = element_text(size = 24)) +
      #ggtitle(obs_ind)+
      geom_point(aes(x = age_grid,y = as.numeric(z_bf[,obs_ind])),size = 10) 
  )
  if(im_save){dev.off()}
}

##### historical regression - breastfeeding #####

beta_save <- mcmc_output$beta_bf_save

# pattern of missingness for breastfeeding
cov_titles <- c("diarrhea","fever","cough")
for(cov_ind in 1:3){
prop_na <- apply(obs_grid_ill[,,cov_ind] == 0,1,mean)
temp_grid <- NULL
is_missing <- NULL
for(m in 2:length(prop_na)){
  temp_grid <- c(temp_grid,rep(age_grid[m],100))
  prop_na_m <- round(prop_na[m]*100)
  is_missing <- c(is_missing,c(rep("missing",prop_na_m),rep("recorded",100-prop_na_m)))
}
if(im_save){png(file=paste0(hd_figs,'proportion_ill_',cov_ind,'_missing.png'), width=480, height=480)}
na_df <- data.frame(age_grid = temp_grid,is_missing = is_missing)
print(
ggplot(na_df) + geom_bar(aes(x = age_grid,fill = is_missing),position="fill") + 
  ylab("missingness") + xlab("age (years)")  + ggtitle(cov_titles[cov_ind])+
  theme(text = element_text(size = 24),legend.title = element_blank())
)
if(im_save){dev.off()}
}

# mean coefficient surface
coefficient_surface_samples <- array(0,dim = c(M-1,M,n_iter))
for(p in 1:dim(beta_save)[1]){
  for(iter in 1:n_iter){
    coefficient_surface_samples[,,iter] <- coefficient_surface_samples[,,iter] + basis_funs[,,p]*beta_save[p,iter]
  }
}
coefficient_surface_mean <- apply(coefficient_surface_samples[,,(n_iter/2):n_iter],c(1,2),mean)
coefficient_surface_mean[coefficient_surface_mean == 0] <- NA

if(im_save){png(file=paste0(hd_figs,'mean_beta_surface_bf.png'), width=480, height=480)}
surf3D(x = mesh(age_grid[2:M],age_grid)$x,y= mesh(age_grid[2:M],age_grid)$y,
       z = coefficient_surface_mean, bty = "g",d = 2,ticktype = "detailed",
       xlab = 't',ylab = 's',zlab = '',shade = 1,colkey = F,phi = 30,theta = 40,
       lighting = T,cex.axis = 1, cex.lab = 2)
if(im_save){dev.off()}



## plot cross sections of covariance surface
#for(age_ind in 1:12){
#  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
#  beta_thin_df <- as.data.frame(cbind(age_grid,y_scale*coefficient_surface_samples[age_ind,,thin_iter]))
#  beta_thin_df <- melt(beta_thin_df,id = "age_grid")
#  if(im_save){png(file=paste0(hd_figs,'beta_samples_age_',age_ind,'.png'), width=480, height=480)}
#  #yl <- c(min(Y_atgpxr[,3]),max(Y_atgpxr[,3]))
#  print(
#    ggplot() + 
#      geom_line(data = beta_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
#      #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
#      ylab("") + xlab("")+ theme(text = element_text(size = 24)) #+ ylim(yl)
#  )
#  if(im_save){dev.off()}
#}

# profiles of breastfeading and weight outcomes
ones_mat <- matrix(1,nrow = M,ncol = M)
ones_mat[lower.tri(ones_mat)] <- 0
bf_profiles <-  ones_mat
n_profiles <- dim(bf_profiles)[2]
profile_integrals <- compute_historical_integrals(bf_profiles,basis_funs,W)
profile_fit <- array(0,dim = c(M-1,n_profiles,n_iter))
for(iter in 1:n_iter){
  profile_fit[,,iter] <- compute_fun_reg_fit(profile_integrals,beta_save[,iter])
}
yl <- c(min(y_scale*profile_fit),max(y_scale*profile_fit))

# plot cross sections of covariance surface
for(age_ind in 1:n_profiles){
  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
  profile_thin_df <- as.data.frame(cbind(age_grid[2:M],y_scale*profile_fit[,age_ind,thin_iter]))
  colnames(profile_thin_df)[1] <- "age_grid"
  profile_thin_df <- melt(profile_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'bf_profile_',age_ind,'_fit.png'), width=480, height=480)}
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
  profile_norm[,iter] <- y_scale*sqrt(diag(t(profile_fit[,,iter])%*%W[2:M,2:M]%*%profile_fit[,,iter]))
}
profile_norm_df <- melt(profile_norm[,(n_iter/2):n_iter])
yl <- c(0,max(profile_norm_df[,3]))
if(im_save){png(file=paste0(hd_figs,'bf_profile_norm_fit.png'), width=480, height=480)}
ggplot() + geom_boxplot(aes(group = as.factor(age_grid[profile_norm_df$Var1]),y = profile_norm_df$value)) +
  scale_x_continuous(breaks=c(-.3455,-.3455/2,0,.3455/2,.3455),
                   labels=c("0",".5", "1","1.5", "2")) +
  ylab("norm") + xlab("age (years)")+ ylim(yl)+ theme(text = element_text(size = 24)) 
if(im_save){dev.off()}

##### trace plots #####

# one observation of beta parameter
for(p in 1:dim(beta_save)[1]){
  if(im_save){png(file=paste0(hd_figs,'beta_bf_',p,'.png'), width=480, height=480)}
  print(
    qplot(1:n_iter,beta_save[p,1:n_iter],geom = "path") + 
      ggtitle("") + xlab("iteration") + ylab("value") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
}

##### concurrent regression  #####

mu_z_save <- mcmc_output$mu_ill_save
l_mu_z_save <- mcmc_output$l_mu_ill_save
scale_mu_z_save <- mcmc_output$scale_mu_ill_save
beta_z_coefs_save <- mcmc_output$beta_ill_coefs_save  

# plot thinned posterior samples of mu_z
yl <- c(0,1)
cov_titles <- c("diarrhea","fever","cough")
for(cov_ind in 1:num_covs){
  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
  mu_thin_df <- as.data.frame(cbind(age_grid[2:M],pnorm(mu_z_save[,cov_ind,thin_iter])))
  colnames(mu_thin_df)[1] <- "age_grid"
  mu_thin_df <- melt(mu_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'mu_ill_',cov_ind,'_samples.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = mu_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ggtitle(cov_titles[cov_ind]) + 
      #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
      ylab("probability") + xlab("age (years)")+ theme(text = element_text(size = 24)) + ylim(yl)
  )
  if(im_save){dev.off()}
}

# plot thind poster samples of beta
yl <- c(-150,50)
cov_titles <- c("diarrhea","fever","cough")
for(cov_ind in 1:num_covs){
  thin_iter <- round(seq((n_iter*.5),n_iter,length.out = 100))
  mu_thin_df <- as.data.frame(cbind(age_grid[2:M],y_scale*beta_ill_basis%*%beta_z_coefs_save[,cov_ind,thin_iter]))
  colnames(mu_thin_df)[1] <- "age_grid"
  mu_thin_df <- melt(mu_thin_df,id = "age_grid")
  if(im_save){png(file=paste0(hd_figs,'beta_ill_',cov_ind,'_samples.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = mu_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .05,size = 2) +  
      ggtitle(cov_titles[cov_ind]) + 
      #geom_point(aes(x = cebu$visit_num/6,y = cebu$weight)) + 
      ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24)) #+ ylim(yl)
  )
  if(im_save){dev.off()}
}


# length-scale of illness means
for(cov_ind in 1:num_covs){
if(im_save){png(file=paste0(hd_figs,'trace_ill_l_mu_',cov_ind,'.png'), width=480, height=480)}
print(
qplot(1:n_iter,l_mu_z_save[cov_ind,],geom = "path") + xlab("iteration") + ylab("value") + 
  ggtitle("") + theme(text = element_text(size = 28)) 
)
if(im_save){dev.off()}
}

# scale of illness means
for(cov_ind in 1:num_covs){
  if(im_save){png(file=paste0(hd_figs,'trace_ill_scale_mu_',cov_ind,'.png'), width=480, height=480)}
  print(
    qplot(1:n_iter,scale_mu_z_save[cov_ind,],geom = "path") + xlab("iteration") + ylab("value") + 
      ggtitle("") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
}

# illness means themselves 
for(cov_ind in 1:num_covs){
  if(im_save){png(file=paste0(hd_figs,'trace_ill_mu_',cov_ind,'_t1.png'), width=480, height=480)}
  print(
    qplot(1:n_iter,pnorm(mu_z_save[1,cov_ind,]),geom = "path") + xlab("iteration") + ylab("value") + 
      ggtitle("") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
  if(im_save){png(file=paste0(hd_figs,'trace_ill_mu_',cov_ind,'_t6.png'), width=480, height=480)}
  print(
    qplot(1:n_iter,pnorm(mu_z_save[6,cov_ind,]),geom = "path") + xlab("iteration") + ylab("value") + 
      ggtitle("") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
  if(im_save){png(file=paste0(hd_figs,'trace_ill_mu_',cov_ind,'_t12.png'), width=480, height=480)}
  print(
    qplot(1:n_iter,pnorm(mu_z_save[12,cov_ind,]),geom = "path") + xlab("iteration") + ylab("value") + 
      ggtitle("") + theme(text = element_text(size = 28)) 
  )
  if(im_save){dev.off()}
}

# trace plots for regression coefficients
for(cov_ind in 1:num_covs){
  for(p in 1:beta_ill_basis_size){
    if(im_save){png(file=paste0(hd_figs,'trace_ill_coefs_cov_',cov_ind,'_p_',p,'.png'), width=480, height=480)}
    print(
    qplot(1:n_iter,beta_z_coefs_save[p,cov_ind,],geom = "path") + xlab("iteration") + ylab("value") + 
      ggtitle("") + theme(text = element_text(size = 28)) 
    )
    if(im_save){dev.off()}
  }
}

##### investigating overall mcmc fit #####

# need to add all aspects together for overall fit 
# plot fitted observations

fit_save <- array(0,dim = c(M,N,n_iter))
for(iter in 1:n_iter){
    # fpca regression
    fit_save[,,iter] <- mcmc_output$mu_y_save[,iter]%*%t(rep(1,N)) + mcmc_output$lambda_y_save[,,iter]%*%mcmc_output$eta_y_save[,,iter] 
    # historical regression
    z_bf_complete <- sample_z_complete(z_bf, mcmc_output$mu_bf_save[,iter]%*%t(rep(1,N)) + mcmc_output$lambda_bf_save[,,iter]%*%mcmc_output$eta_bf_save[,,iter],obs_grid_bf)
    historical_integrals <- compute_historical_integrals(z_bf_complete,hist_basis,W)
    fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_fun_reg_fit(historical_integrals,mcmc_output$beta_bf_save[,iter]))
    # concurrent regression 
    for(cov_ind in 1:num_covs){
      z_ill_complete[,,cov_ind] <- sample_z_complete(z_ill[,,cov_ind], mcmc_output$mu_ill_save[,cov_ind,iter]%*%t(rep(1,N)),obs_grid_ill[,,cov_ind])
    }
    fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_concurrent_fit(z_ill_complete,beta_ill_basis,mcmc_output$beta_ill_coefs_save[,,iter]))
}


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
  if(im_save){png(file=paste0(hd_figs,'fit_',obs_ind,'.png'), width=480, height=480)}
  print(
    ggplot() + 
      geom_line(data = fit_thin_df, aes(x = age_grid, y = value, group = variable),alpha = .025,size = 2) +  
      ylab("weight (grams)") + xlab("age (years)")+ theme(text = element_text(size = 24)) +
      #ggtitle(obs_ind)+
      geom_point(aes(x = cebu_all[cebu_all$id == obs_ind,2],y = cebu_all[cebu_all$id == obs_ind,3]), size = 5) 
  )
  if(im_save){dev.off()}
}

##### residual analysis #####

im_save <- FALSE

fit_post_mean <- apply(fit_save[,,(n_iter/2):n_iter],c(1,2),mean)

residual_age <- NULL
residual_weight <- NULL
predicted_weight <- NULL

for(i in 1:N){
  temp_ind <- which(obs_grid_y[,i] == 1)
  residual_age <- c(residual_age,age_grid[temp_ind])
  residual_weight <- c(residual_weight,y_scale*(y_common[temp_ind,i] - fit_post_mean[temp_ind,i]))
  predicted_weight <- c(predicted_weight,y_scale*fit_post_mean[temp_ind,i])
}

y <- quantile(residual_weight, c(0.1, 0.9))
x <- qnorm(c(0.1, 0.9))
slope <- diff(y)/diff(x)
int <- y[1L] - slope * x[1L]

if(im_save){png(file=paste0(hd_figs,'residuals_by_pred_weight.png'), width=480, height=480)}
ggplot() + geom_boxplot(mapping = aes(group = residual_age,y = residual_weight),alpha = .05) + 
  ylab("weight (grams)") + xlab("age (years)")+theme(text = element_text(size = 24))  + 
  scale_x_continuous(breaks=c(-.3455,-.3455/2,0,.3455/2,.3455),labels=c("0",".5", "1","1.5", "2"))
if(im_save){dev.off()}


##### subject-specific posterior-predictive checks #####

# norms of individual functios
y_pred_norm <- array(0,dim = c(N,n_iter))
y_norm <- array(0,dim = c(N))
for(i in 1:N){
  temp_ind <- which(obs_grid_y[,i] == 1)
  y_temp <- y_common[temp_ind,i]
  y_pred <- array(0,dim = c(length(y_temp),n_iter))
  for(m in 1:length(temp_ind)){
    y_pred[m,] <- rnorm(1000,mean = fit_save[temp_ind[m],i,],sd = sqrt(sig_sq_save))
  }
  y_pred_norm[i,] <- sqrt(diag(t(y_pred)%*%y_pred))
  y_norm[i] <- sqrt(t(y_temp)%*%y_temp)
}

post_burn_iter <- ((n_iter/2)+1):n_iter
# worst case scenario
norm_order <- order(abs(apply(y_pred_norm[,post_burn_iter],1,mean) - y_norm),decreasing = TRUE)
norm_order[1:3]
for(temp_ind in 1:3){
  i <- norm_order[temp_ind]
  if(im_save){png(file=paste0(hd_figs,'obs_norm_worst_',i,'.png'), width=480, height=480)}
  print(
    ggplot() + geom_histogram(aes(x = y_scale*y_pred_norm[i,post_burn_iter],y = ..density..),
                              color="black", fill="gray") + xlab('norm')+
      ggtitle(paste0("i = ",i)) + theme(text = element_text(size = 24)) +
      geom_vline(xintercept = y_scale*y_norm[i],color = 'red',linetype = 'dashed',size = 2)
  )
  if(im_save){dev.off()}
}

# random subjects
rand_subjects <- c(716, 1943,  219)
for(temp_ind in 1:3){
  i <- rand_subjects[temp_ind]
  if(im_save){png(file=paste0(hd_figs,'obs_norm_rand_',i,'.png'), width=480, height=480)}
  print(
    ggplot() + geom_histogram(aes(x = y_scale*y_pred_norm[i,post_burn_iter],y = ..density..),
                              color="black", fill="gray") + xlab('norm')+
      ggtitle(paste0("i = ",i)) + theme(text = element_text(size = 24)) +
      geom_vline(xintercept = y_scale*y_norm[i],color = 'red',linetype = 'dashed',size = 2)
  )
  if(im_save){dev.off()}
}

norm_diff <- abs(apply(y_pred_norm[,post_burn_iter],1,mean) - y_norm)
if(im_save){png(file=paste0(hd_figs,'obs_norm_means.png'), width=480, height=480)}
ggplot() + geom_histogram(aes(x = y_scale*norm_diff,y = ..density..),
                          color="black", fill="gray") + xlab('norm')+
  theme(text = element_text(size = 24)) 
if(im_save){dev.off()}

##### population level quantities #####

# population level mean
y_pred_mean <- array(0,dim = c(M,n_iter))
for(m in 1:M){
  temp_ind <- which(obs_grid_y[m,] == 1)
  y_temp <- y_common[m,temp_ind]
  y_pred <- array(0,dim = c(length(y_temp),n_iter))
  for(i in 1:length(temp_ind)){
    y_pred[i,] <- rnorm(1000,mean = fit_save[m,temp_ind[i],],sd = sqrt(sig_sq_save))
  }
  y_pred_mean[m,] <- apply(y_pred,2,mean)
}

y_common_NA <- y_common
y_common_NA[y_common_NA == 0] <- NA

y_obs_mean <- apply(y_common_NA,1,mean,na.rm = T)

for(m in 1:M){
  if(im_save){png(file=paste0(hd_figs,'pop_mean_time_',m,'.png'), width=480, height=480)}
  print(
    ggplot() + geom_histogram(aes(x = y_scale*y_pred_mean[m,post_burn_iter],y = ..density..),
                              color="black", fill="gray") + xlab(expression(mu(t)))+
      theme(text = element_text(size = 24)) +
      geom_vline(xintercept = y_scale*y_obs_mean[m],color = 'red',linetype = 'dashed',size = 2)
  )
  if(im_save){dev.off()}
}

# population covariance
y_pred_eig <- array(0,dim = c(5,n_iter))
for(iter in 1:n_iter){
  y_pred <- array(NA,dim = c(M,N))
  for(i in 1:N){
    for(m in 1:M){
      if(obs_grid_y[m,i]>0){
        y_pred[m,i] <- rnorm(1,mean = fit_save[m,i,iter],sd = sqrt(sig_sq_save[iter]))
      }
    }
  }
  y_pred_eig[,iter] <- eig(cov(t(y_pred),use = "pairwise.complete.obs"))[1:5]
}
y_eig <- eig(cov(t(y_common_NA),use = "pairwise.complete.obs"))[1:5]

for(e in 1:5){
  if(im_save){png(file=paste0(hd_figs,'cov_eig_',e,'.png'), width=480, height=480)}
  print(
    ggplot() + geom_histogram(aes(x = y_scale^2*y_pred_eig[e,post_burn_iter],y = ..density..),
                              color="black", fill="gray") + xlab('eigenvalue')+
      theme(text = element_text(size = 24)) + scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
      theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
      geom_vline(xintercept = y_scale^2*y_eig[e],color = 'red',linetype = 'dashed',size = 2)
  )
  if(im_save){dev.off()}
}

##### compute WAIC for model #####

# for ReMO output
n_iter <- 1000
fit_save <- array(0,dim = c(M,N,n_iter))
grid_weight <- diff(c(age_grid[1],(age_grid[2:M]+age_grid[1:(M-1)])/2,age_grid[M]))
grid_weight <- diag(grid_weight)
for(iter in 1:n_iter){
  # fpca regression
  fit_save[,,iter] <- mcmc_output$mu_y_save[,iter]%*%t(rep(1,N)) + mcmc_output$lambda_y_save[,,iter]%*%mcmc_output$eta_y_save[,,iter] 
  # historical regression
  z_bf_complete <- sample_z_complete(z_bf, mcmc_output$mu_bf_save[,iter]%*%t(rep(1,N)) + mcmc_output$lambda_bf_save[,,iter]%*%mcmc_output$eta_bf_save[,,iter],obs_grid_bf)
  historical_integrals <- compute_historical_integrals(z_bf_complete,hist_basis,grid_weight)
  fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_fun_reg_fit(historical_integrals,mcmc_output$beta_bf_save[,iter]))
  # concurrent regression 
  for(cov_ind in 1:num_covs){
    z_ill_complete[,,cov_ind] <- sample_z_complete(z_ill[,,cov_ind], mcmc_output$mu_ill_save[,cov_ind,iter]%*%t(rep(1,N)),obs_grid_ill[,,cov_ind])
  }
  fit_save[,,iter] <- fit_save[,,iter] + rbind(0*one_vec_N,compute_concurrent_fit(z_ill_complete,beta_ill_basis,mcmc_output$beta_ill_coefs_save[,,iter]))
}


thin_iter <- round(seq(n_iter/2 +1,n_iter,length.out = 500))
waic_vec_remo <- NULL
lppd <- 0
p_waic <- 0
for(i in 1:N){
  temp_ind <- which(obs_grid_y[,i] == 1)
  lik_temp <- array(0,dim = c(length(temp_ind),length(thin_iter)))
  for(m in 1:length(temp_ind)){
    lik_temp[m,] <- dnorm(y_common[temp_ind[m],i],
                          mean = fit_save[temp_ind[m],i,thin_iter],
                          sd = sqrt(mcmc_output$sig_sq_save[thin_iter]))
  }
  lppd <- lppd + sum(log(apply(lik_temp,1,mean)))
  p_waic <- p_waic + sum(apply(log(lik_temp),1,var))
  waic_vec_remo <- c(waic_vec_remo,log(apply(lik_temp,1,mean)) - apply(log(lik_temp),1,var))
}
waic_remo <- -2*(lppd - p_waic)
