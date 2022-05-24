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

# initialize historical regression
num_centers <- 4
source("create_tent_basis.R")

hist_basis <- basis_funs
P <- dim(hist_basis)[3]
beta_bf_cur <- rnorm(P,0,1)
z_bf_complete <- array(0,dim = dim(z_bf))
historical_integrals <- compute_historical_integrals(z_bf_complete,hist_basis,W)
historical_fit <- rbind(0*one_vec_N,compute_fun_reg_fit(historical_integrals,beta_bf_cur))

beta_ill_basis_size <- 4
beta_ill_basis <- bs(age_grid[2:M],df = beta_ill_basis_size,intercept = T)
beta_ill_basis <- matrix(beta_ill_basis,nrow = M-1)
beta_ill_coefs_cur <- array(0,dim = c(beta_ill_basis_size,num_covs))

z_ill_complete <- array(0,dim = dim(z_ill))
concurrent_fit <- rbind(0*one_vec_N,compute_concurrent_fit(z_ill_complete,beta_ill_basis,beta_ill_coefs_cur))


##### Impute longitudinal covaraites #####

z_bf <- z_bf[,1:N]

# impute bf indicator 
Ly <- list()
Lt <- list()
num_obs <- rep(0,N)
for(i in 1:N){
  inds_obs <- which(obs_grid_bf[,i] > 0)
  Ly[[i]] <- z_bf[inds_obs,i]
  Lt[[i]] <- age_grid[inds_obs]
  num_obs[i] <- length(inds_obs)
}

Ly <- Ly[num_obs>1]
Lt <- Lt[num_obs>1]

K_bf <- 3
slfpca_results_bf <- SLFPCA(Ly,Lt,interval = c(min(age_grid),max(age_grid)),npc = K_bf, L_list = c(5),
                            norder = 4, kappa_theta = 0.2, sparse_pen = 0,
                            nRegGrid = length(age_grid), stepmu = 0.005)

slfpca_mu <- eval.fd(age_grid,slfpca_results_bf$mufd)
slfpca_eig_vec <- array(0,dim = c(M,K_bf)) 
for(k in 1:K_bf){
  slfpca_eig_vec[,k] <- eval.fd(age_grid,slfpca_results_bf$eigfd_list[[k]])
}
slfpca_score <- slfpca_results_bf$score
fit_slfpca <- inv.logit(slfpca_mu%*%t(rep(1,length(Lt))) + slfpca_eig_vec%*%t(slfpca_score))

fit_slfpca_all <- array(NA,dim = c(length(age_grid),N))
fit_slfpca_all[,num_obs<=1] <- inv.logit(slfpca_mu)
fit_slfpca_all[,num_obs>1] <- fit_slfpca

z_bf_complete <- z_bf
# randomly impute
z_bf_complete[is.na(z_bf_complete)] <- runif(sum(is.na(z_bf_complete)))<fit_slfpca_all[is.na(z_bf_complete)]

# impute illness indicators
z_ill <- z_ill[,1:N,]
z_ill_complete <- z_ill
slfpca_results_illness_all <- list()

for(cov_ind in 1:num_covs){
  Ly <- list()
  Lt <- list()
  num_obs <- rep(0,N)
  for(i in 1:N){
    inds_obs <- which(obs_grid_ill[,i,cov_ind] > 0)
    Ly[[i]] <- z_ill[inds_obs,i,cov_ind]
    Lt[[i]] <- age_grid[inds_obs]
    num_obs[i] <- length(inds_obs)
  }
  Ly <- Ly[num_obs>1]
  Lt <- Lt[num_obs>1]
  K_ill <- 1
  slfpca_results_ill <- SLFPCA(Ly,Lt,interval = c(min(age_grid[2:M]),max(age_grid[2:M])),npc = K_ill, L_list = c(5),
                               norder = 4, kappa_theta = 0.2, sparse_pen = 0,
                               nRegGrid = length(age_grid[2:M]), stepmu = 0.005)
  slfpca_results_illness_all[[cov_ind]] <- slfpca_results_ill
  
  
  slfpca_mu <- eval.fd(age_grid[2:M],slfpca_results_ill$mufd)
  fit_slfpca_all <- repmat(inv.logit(slfpca_mu),1,N)
  
  z_ill_complete[,,cov_ind] <- z_ill[,,cov_ind]
  z_ill_complete[,,cov_ind][is.na(z_ill_complete[,,cov_ind])] <- runif(sum(is.na(z_ill_complete[,,cov_ind])))<fit_slfpca_all[is.na(z_ill_complete[,,cov_ind])]
}

##### after imputation, optimize for FPCA while accounting for longitudinal covariates #####

y_common <- y_common[,1:N]
obs_grid_y <- obs_grid_y[,1:N]

# initial estimate of mean based on population mean 
y_common_NA <- y_common
y_common_NA[y_common_NA == 0] <- NA 

mu_obs_grid <- apply(y_common_NA,1,mean,na.rm = T)
fit_pace <- mu_obs_grid%*%t(rep(1,N))

# initialize current optimal values
historical_integrals <- compute_historical_integrals(z_bf_complete,hist_basis,W)
historical_fit <- rbind(0*rep(1,N),compute_fun_reg_fit(historical_integrals,beta_bf_cur))
concurrent_integrals <- array(0,dim = c(M-1,beta_ill_basis_size,N,num_covs))
for(cov_ind in 1:num_covs){
  concurrent_integrals[,,,cov_ind] <- compute_fun_prod(z_ill_complete[,,cov_ind],beta_ill_basis)
}
concurrent_fit <- rbind(0*rep(1,N),compute_concurrent_fit(z_ill_complete,beta_ill_basis,beta_ill_coefs_cur))

##### Everything is set for WinBUGS #####
grid_weight <- diff(c(age_grid[1],(age_grid[2:M]+age_grid[1:(M-1)])/2,age_grid[M]))
grid_weight <- diag(grid_weight)
# set variables for WinBUGS
N_subj=dim(y_common)[2]
N_obs=dim(y_common)[1]
dim.space=5

# W is y - mean(y) 
y <- y_common
y[y == 0] = NA

W=t(y)
W<-W-matrix(rep(colMeans(W,na.rm=TRUE),nrow(W)),nrow=nrow(W),byrow=TRUE)

scalar_covs <- t(x)

# Estimate covariance function through bivaraite smoothing
Gt=cov(W,use="pairwise.complete.obs")
gw.temp <- Gt
diag(gw.temp) <- rep(NA, length(diag(Gt)))

data.gb <- data.frame(gw=as.vector(gw.temp),
                      x1=rep(seq(0,1,length=length(diag(Gt))), length(diag(Gt))), x2=rep(seq(0,1,length=length(diag(Gt))), each=length(diag(Gt))))
x1 <- data.gb$x1
x2 <- data.gb$x2
gw <- data.gb$gw

library(SemiPar)
myknots <- data.frame(x1=rep(seq(0,1,length=10),each=10), x2=rep(seq(0,1,length=10),10)) 
fit<- spm(gw ~ f(x1, x2,knots=myknots),omit.missing=T)        


newdata <- data.frame(x1=x1,x2=x2)
pred <- predict(fit,newdata)

var.noise <- mean( diag(Gt) - diag(matrix(pred,length(diag(Gt)))))

s.gt <-matrix(pred,length(diag(Gt))) 
sm_Gt <- (s.gt + t(s.gt) )/2 

#Obtain the eigenvectors of Gt up to the dim.space
E=(eigen(sm_Gt)$vectors[1:N_obs,1:dim.space])
E_norm <- sqrt(diag(t(E)%*%grid_weight%*%E))
E <- E/rep(1,M)%*%t(E_norm)

E_vals <- 1/(eigen(sm_Gt)$values[1:5]*(E_norm^2))

##### MCMC Using WinBUGS #####
P_scalar <- dim(x)[1]
P_bf <- dim(historical_integrals)[2]
P_ill <- dim(concurrent_integrals)[2]

data<-list("E","E_vals","W","N_subj","N_obs","dim.space",
           "P_scalar","scalar_covs",
           "P_bf","historical_integrals",
           "P_ill","concurrent_integrals")

#This is the name of the program
program.file.name="cebu_wb_model.txt"

#Define the initial values & parameters to record 
inits.W=matrix(rep(NA,N_subj*N_obs),ncol=N_obs)
inits.W[is.na(W)]=mean(mean(W,na.rm=TRUE))
inits<-function(){list(xi=matrix(rep(0,N_subj*dim.space),ncol=dim.space),
                       beta_scalar=array(0,dim = c(dim.space,4)),
                       beta_bf=rep(0,P_bf),
                       beta_ill=array(0,dim = c(4,3)),
                       taueps=0.01,W=inits.W)}

parameters=list("beta_ill","beta_bf","beta_scalar","xi","taueps")


# MCMC
n.thin=1
n.iter=1500
n.burnin=500

library(R2WinBUGS)
wb_output <- bugs(data, inits, parameters, model.file = program.file.name,
                 n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin,debug = FALSE, DIC = FALSE, digits = 5, codaPkg = FALSE,
                 bugs.directory = "c:/Program Files/WinBUGS14/")

save.image('./results/cebu_wb_mcmc.RData')


