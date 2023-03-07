library(MASS)
library(pracma)
library(fdapace)
library(fpca)
library(Rcpp)
library(SemiPar)
library(fda)
library(rootSolve)
library(reshape2)
library(ggplot2)
library(R2WinBUGS)
library(ggpubr)
library(grid)

sourceCpp("../src/nemo_ffa.cpp")
sourceCpp("../src/msf.cpp")
source("../src/nemo_ffa.R")
source("./supp_funs/fpca_funs.R")

##### FPCA of Growth curves From PACE #####

growth_data <- growth$hgtm
time <- growth$age
N <- dim(growth_data)[2]

growth_df <- cbind(melt(as.data.frame(growth_data)),rep(time,N))
colnames(growth_df)[3] <- "age"

ggplot() + geom_line(growth_df,mapping = aes(x = age,y = value,color = variable),size = 1) + 
  ylab("height (cm)")+ xlab("age (years)") + theme(text = element_text(size = 24),legend.position = 'none')

Ly <- list()
Lt <- list()
for(i in 1:N){
  Ly[[i]] <- growth_data[,i]
  Lt[[i]] <- time
}

pace_options <- list()
pace_options$dataType <- 'Dense'
pace_options$plot <- TRUE
pace_options$FVEthreshold <- .95

pace_results <- FPCA(Ly,Lt,optns = pace_options)

mu_true <- ConvertSupport(
  pace_results$workGrid,
  pace_results$obsGrid,
  phi = pace_results$mu
)

lambda_true <- ConvertSupport(
  pace_results$workGrid,
  pace_results$obsGrid,
  phi = pace_results$phi
)

eta_var_true <- pace_results$lambda

sig_true <- sqrt(pace_results$sigma2)

set.seed(1234)


n_rep <- 100
M <- length(time)
one_M <- rep(1,M)
one_N <- rep(1,N)
y_long_save <- list()
y_long_count <- 1

mu_MISE <- array(NA,dim = c(4,n_rep))
lambda_1_MISE <- array(NA,dim = c(4,n_rep))
lambda_2_MISE <- array(NA,dim = c(4,n_rep))
lambda_3_MISE <- array(NA,dim = c(4,n_rep))
est_MISE <- array(NA,dim = c(4,n_rep))
K_selected <- array(NA,dim = c(4,n_rep))


for(replicate in 1:n_rep){

##### generate underlying observations #####
f <- mu_true%*%t(rep(1,N)) + lambda_true%*%t(mvrnorm(N,mu = rep(0,length(eta_var_true)),Sigma = diag(eta_var_true)))
noise <- matrix(rnorm(N*length(time),0,sig_true),nrow = length(time),ncol = N)
y <- f + noise

obs_grid <- matrix(rep(1,M*N),nrow = M,ncol = N)
y <- f + noise
y_common <- y*obs_grid

y_long <- NULL
for(i in 1:N){
  id_temp <- rep(i,sum(obs_grid[,i]))
  t_temp <- time[which(obs_grid[,i]>0)]
  y_temp <- y_common[which(obs_grid[,i]>0),i]
  y_long <- rbind(y_long,cbind(id_temp,t_temp,y_temp))
}

##### run nemo ffa #####
y_sc <- sqrt(norm(y_common))
y_common <- y_common/y_sc
w <- diff(c(time[1],(time[2:M]+time[1:(M-1)])/2,time[M]))
W <- diag(w)
one_vec_N <- rep(1,N)

##### run nemo ffa #####
K_max <- 5
nu_eta <- 1e-4
nu_lambda <- 1e-4

inv_g_hyper <- length_scale_hyper(time)
prior_a <- inv_g_hyper[1]
prior_b <- inv_g_hyper[2]
prior_var <- 1
n_iter <- 1000
init_mcmc_params <- init_ffa_mcmc(y_common,obs_grid,time,
                                  K_max,nu_eta,nu_lambda,
                                  prior_a,prior_b,prior_var,
                                  n_iter)
n_iter <- 1000
lag <- 1
ffa_mcmc <- run_ffa_mcmc(y_common,obs_grid,time,
                         nu_eta,nu_lambda,
                         init_mcmc_params,
                         prior_a,prior_b,prior_var,
                         n_iter,lag)

K_selected[1,replicate] <- dim(init_mcmc_params$eta_cur)[1]

# compute mu MISE
mu_post_mean <- apply(ffa_mcmc$mu_save,1,mean)*y_sc
mu_MISE_nemo <- diag(t(mu_true - mu_post_mean)%*%W%*%(mu_true - mu_post_mean))

# compute lambda MISE
# resolve label switching an sign ambiguity
lambda_save_processed <- array(0,dim = c(M,3,n_iter/lag))
for(iter in 1:dim(lambda_save_processed)[3]){
  lambda_post <- ffa_mcmc$lambda_save[,,iter]*(one_M%*%t(sqrt(ffa_mcmc$psi_save[,iter])))
  lambda_temp_norm <- sqrt(diag(t(lambda_post)%*%W%*%lambda_post))
  lambda_order <- order(lambda_temp_norm,decreasing = T)
  lambda_post <- lambda_post[,lambda_order[1:3]]/(rep(1,M)%*%t(lambda_temp_norm[lambda_order[1:3]]))
  perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
  lambda_save_processed[,,iter] <- aplr(lambda_post,perm_out)
}
lambda_post_mean <- apply(lambda_save_processed,c(1,2),mean)
if(dim(lambda_post_mean)[2]<3){
  lambda_post_mean <- cbind(lambda_post_mean,
                            array(NA,dim = c(M,3 - dim(lambda_post_mean)[2])))
}

lambda_MISE_nemo <- diag(t(lambda_true - lambda_post_mean)%*%W%*%(lambda_true - lambda_post_mean))

# overall fit
n_iter <- n_iter/lag
fit_save <- array(0,dim = c(M,N,n_iter/2))
for(iter in (n_iter/2 + 1):n_iter){
  if(K_selected[1,replicate] == 1){
    fit_save[,,iter - n_iter/2] <- ffa_mcmc$mu_save[,iter]%*%t(one_N) + ffa_mcmc$lambda_save[,,iter]%*%t(ffa_mcmc$eta_save[,,iter])
  }else{
    fit_save[,,iter - n_iter/2] <- ffa_mcmc$mu_save[,iter]%*%t(one_N) + ffa_mcmc$lambda_save[,,iter]%*%ffa_mcmc$eta_save[,,iter]
  }
  
}
fit_mean <- apply(fit_save,c(1,2),mean)*y_sc
MISE <- diag(t(fit_mean - f)%*%W%*%(fit_mean - f))

# store estimated quantities

mu_MISE[1,replicate] <- mu_MISE_nemo
lambda_1_MISE[1,replicate] <- lambda_MISE_nemo[1]
lambda_2_MISE[1,replicate] <- lambda_MISE_nemo[2]
lambda_3_MISE[1,replicate] <- lambda_MISE_nemo[3]
est_MISE[1,replicate] <- mean(MISE)

y_common <- y_common*y_sc

##### Run PACE #####

Ly <- list()
Lt <- list()
for(i in 1:N){
  inds_obs <- which(obs_grid[,i] > 0)
  Ly[[i]] <-y_common[inds_obs,i]
  Lt[[i]] <- time[inds_obs]
}

pace_options <- list()
pace_options$dataType <- 'Dense'
pace_options$methodSelectK <- 3

pace_results <- FPCA(Ly,Lt,optns = pace_options)

mu_obs_grid <- ConvertSupport(
  pace_results$workGrid,
  pace_results$obsGrid,
  phi = pace_results$mu
)

phi_obs_grid <- ConvertSupport(
  pace_results$workGrid,
  pace_results$obsGrid,
  phi = pace_results$phi
)

# compute mu MISE
mu_MISE_pace <- diag(t(mu_true - mu_obs_grid)%*%W%*%(mu_true - mu_obs_grid))

# compute lambda MISE
lambda_post <- phi_obs_grid
if(pace_results$selectK > 3){
  lambda_post <- lambda_post[,1:3]
}
perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
lambda_post <- aplr(lambda_post,perm_out)
if(pace_results$selectK<3){
  lambda_post <- cbind(lambda_post,
                           array(NA,dim = c(M,3 -pace_results$selectK)))
}
lambda_MISE_pace <- diag(t(lambda_true - lambda_post)%*%W%*%(lambda_true - lambda_post))

# overall fit
fit_pace <- mu_obs_grid%*%t(one_N) + phi_obs_grid%*%t(pace_results$xiEst)

MISE_pace <- diag(t(fit_pace - f)%*%W%*%(fit_pace - f))

# store estimated quantities
K_selected[2,replicate] <- pace_results$selectK
mu_MISE[2,replicate] <- mu_MISE_pace
lambda_1_MISE[2,replicate] <- lambda_MISE_pace[1]
lambda_2_MISE[2,replicate] <- lambda_MISE_pace[2]
lambda_3_MISE[2,replicate] <- lambda_MISE_pace[3]
est_MISE[2,replicate] <- mean(MISE_pace)

##### Peng and Paul FPCA ##### commented out because of poor performance
#fpca_data <- cbind(y_long[,1],y_long[,3],y_long[,2]/max(y_long[,2]))
#fpca_results <- fpca.mle(fpca_data,c(5,10,15,20),c(2,3,4,5))
#
#tryCatch({ # There  sometimes when this code throws an error
#  fpca.score2 <- function(data.m, grids.u, muhat, eigenvals, eigenfuncs, sig2hat,K) 
#  {
#    temp <- table(data.m[, 1])
#    n <- length(temp)
#    m.l <- as.vector(temp)
#    result <- matrix(0, n, K)
#    N <- length(grids.u)
#    evalmat <- diag(eigenvals[1:K])
#    current <- 0
#    eigenfuncs.u <- t(eigenfuncs)
#    data.u <- matrix(as.numeric(as.vector(data.m[, -1])), nrow = nrow(data.m[, 
#                                                                             -1]), ncol = ncol(data.m[, -1]))
#    for (i in 1:n) {
#      Y <- as.vector(data.u[(current + 1):(current + m.l[i]), 
#                            1])
#      meastime <- data.u[(current + 1):(current + m.l[i]), 
#                         2]
#      gridtime <- ceiling(N * meastime) 
#      gridtime[gridtime == 0] <- 1 # there is an off by 1 error in the R package
#      muy <- muhat[gridtime]
#      Phiy <- matrix(eigenfuncs.u[gridtime, 1:K], ncol = K)
#      Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
#      temp.y <- matrix(Y - muy)
#      result[i, ] <- evalmat %*% t(Phiy) %*% solve(Sigy, temp.y)
#      current <- current + m.l[i]
#    }
#    return(result)
#  }
#  
#  fpca_scores <- fpca.score2(fpca_data,fpca_results$grid,fpca_results$fitted_mean,
#                             fpca_results$eigenvalues,fpca_results$eigenfunctions,fpca_results$error_var,fpca_results$selected_model[2])
#  mu_dense <- fpca_results$fitted_mean
#  lambda_dense <- fpca_results$eigenfunctions
  
#  mu_og_grid <- approx(x = fpca_results$grid,mu_dense,xout =time/max(y_long[,2]))$y
#  lambda_og_grid <- array(0,dim = c(length(mu_og_grid),fpca_results$selected_model[2]))
#  for(k in 1:fpca_results$selected_model[2]){
#    lambda_og_grid[,k] <- approx(x = fpca_results$grid,lambda_dense[k,],xout =time/max(y_long[,2]))$y
#  }
  
#  fit_fpca_dense <- fpca_results$fitted_mean%*%t(rep(1,N)) + t(fpca_results$eigenfunctions)%*%t(fpca_scores)
#  
#  fit_fpca_og <- array(0,dim = c(M,N))
#  for(i in 1:N){
#    fit_fpca_og[,i] <- approx(x = fpca_results$grid,fit_fpca_dense[,i],xout = time/max(y_long[,2]))$y
#  }
  
#  MISE_fpca <- diag(t(fit_fpca_og - f)%*%W%*%(fit_fpca_og - f))
  
  # compute mu MISE
#  mu_MISE_fpca <- diag(t(mu_true - mu_og_grid)%*%W%*%(mu_true - mu_og_grid))
  
  # compute lambda MISE
#  lambda_temp_norm <- sqrt(diag(t(lambda_og_grid)%*%W%*%lambda_og_grid))
#  lambda_post <- lambda_og_grid/(one_M%*%t(lambda_temp_norm))
#  if(fpca_results$selected_model[2] > 3){
#    lambda_post <- lambda_post[,1:3]
#  }
#  perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
#  lambda_post <- aplr(lambda_post,perm_out)
#  if(fpca_results$selected_model[2]<3){
#    lambda_post <- cbind(lambda_post,
#                         array(NA,dim = c(M,3 -pace_results$selectK)))
#  }
#  lambda_MISE_fpca <- diag(t(lambda_true - lambda_post)%*%W%*%(lambda_true - lambda_post))
  
#  # store estimated quantities
#  K_selected[3,replicate] <- fpca_results$selected_model[2]
#  mu_MISE[3,replicate] <- mu_MISE_fpca
#  lambda_1_MISE[3,replicate] <- lambda_MISE_fpca[1]
#  lambda_2_MISE[3,replicate] <- lambda_MISE_fpca[2]
#  lambda_3_MISE[3,replicate] <- lambda_MISE_fpca[3]
#  est_MISE[3,replicate] <- mean(MISE_fpca)
#},error = function(cond){
#  K_selected[3,replicate] <- NA
#  mu_MISE[3,replicate] <- NA
#  lambda_1_MISE[3,replicate] <- NA
#  lambda_2_MISE[3,replicate] <- NA
#  lambda_3_MISE[3,replicate] <- NA
#  est_MISE[3,replicate] <- NA
#})

##### Run GC model using WinBUGS #####

time_grid<- time
grid_weight <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
grid_weight <- diag(grid_weight)

y <- as.matrix(dcast(as.data.frame(y_long),id_temp~t_temp,value.var = "y_temp"))
y <- y[,2:(M+1)]

# set variables for WinBUGS
N_subj=dim(y)[1]
N_obs=dim(y)[2]
dim.space=3

# W is y - mean(y) 
W=y
W<-W-matrix(rep(colMeans(W,na.rm=TRUE),nrow(W)),nrow=nrow(W),byrow=TRUE)

# Estimate covariance function through bivaraite smoothing
Gt=cov(W,use="pairwise.complete.obs")
gw.temp <- Gt
diag(gw.temp) <- rep(NA, length(diag(Gt)))

data.gb <- data.frame(gw=as.vector(gw.temp),
                      x1=rep(seq(0,1,length=length(diag(Gt))), length(diag(Gt))), x2=rep(seq(0,1,length=length(diag(Gt))), each=length(diag(Gt))))

x1 <- data.gb$x1
x2 <- data.gb$x2
gw <- data.gb$gw


myknots <- data.frame(x1=rep(seq(0,1,length=10),each=10), x2=rep(seq(0,1,length=10),10)) 
fit<- spm(gw ~ f(x1, x2,knots=myknots),omit.missing=T)        

newdata <- data.frame(x1=x1,x2=x2)
pred <- predict(fit,newdata)


var.noise <- mean( diag(Gt) - diag(matrix(pred,length(diag(Gt)))))

s.gt <-matrix(pred,length(diag(Gt))) 
sm_Gt <- (s.gt + t(s.gt) )/2 

#Obtain the eigenvectors of Gt up to the dim.space
E=(eigen(sm_Gt)$vectors[1:N_obs,1:dim.space])

##### MCMC Using WinBUGS #####
data<-list("E","W","N_subj","N_obs","dim.space")

#This is the name of the program
program.file.name="wb_growth_model.txt"

#Define the initial values & parameters to record 
inits.W=matrix(rep(NA,N_subj*N_obs),ncol=N_obs)
inits.ll=rep(0.01,dim.space)
inits.W[is.na(W)]=mean(mean(W,na.rm=TRUE))
inits<-function(){list(xi=matrix(rep(0,N_subj*dim.space),ncol=dim.space),
                       taueps=0.01,ll=inits.ll,W=inits.W)}

parameters=list("lambda","xi")

# MCMC
n.thin=1
n.iter=1500
n.burnin=500



Bayes.fit<- bugs(data, inits, parameters, model.file = program.file.name,
                 n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
                 n.thin = n.thin,debug = FALSE, DIC = FALSE, digits = 5, codaPkg = FALSE,
                 bugs.directory = "c:/Program Files/WinBUGS14/")


y_mean=colMeans(y[1:N_subj,1:N_obs],na.rm=TRUE)
wb_fit <- t(Bayes.fit$mean$xi%*%t(E) + rep(1,dim(Bayes.fit$mean$xi)[1])%*%t(y_mean))

wb_MISE <- diag(t(wb_fit - f)%*%grid_weight%*%(wb_fit - f))

# compute mu MISE
y_mean=colMeans(y[1:N_subj,1:N_obs],na.rm=TRUE)
mu_MISE_cg <- diag(t(mu_true - y_mean)%*%grid_weight%*%(mu_true - y_mean))

# compute lambda MISE
lambda_temp_norm <- sqrt(diag(t(E)%*%grid_weight%*%E))
lambda_post <- E/(one_M%*%t(lambda_temp_norm))
perm_out <- msfOUT_functional(lambda_post,grid_weight,lambda_true)
lambda_post <- aplr(lambda_post,perm_out)
lambda_MISE_cg <- diag(t(lambda_true - lambda_post)%*%grid_weight%*%(lambda_true - lambda_post))

#overall fit
wb_fit <- t(Bayes.fit$mean$xi%*%t(E) + rep(1,dim(Bayes.fit$mean$xi)[1])%*%t(y_mean))
wb_MISE <- diag(t(wb_fit - f)%*%grid_weight%*%(wb_fit - f))

# store estimated quantities
K_selected[4,replicate] <- 3
mu_MISE[4,replicate] <- mu_MISE_cg
lambda_1_MISE[4,replicate] <- lambda_MISE_cg[1]
lambda_2_MISE[4,replicate] <- lambda_MISE_cg[2]
lambda_3_MISE[4,replicate] <- lambda_MISE_cg[3]
est_MISE[4,replicate] <- mean(wb_MISE)

print(replicate)
print(K_selected[,replicate])
print(mu_MISE[,replicate])
print(lambda_1_MISE[,replicate])
print(est_MISE[,replicate])

}

##### Visualize results #####

# accidentally saved incorrectly, so resolving manuyally
cube_diff <- function(cube_array){
  diff_array <- array(0,dim= c(3,2,100))
  for(method in 1:2){
    diff_array[,method,] <- cube_array[,1+method] - cube_array[,1]
  }
  return(diff_array)
}

reorg_results <- function(cube_array){
  s_levels <- c(.25,.5,.75)
  m_levels <- c("YMW2005","CG2010")
  long_rep <- NULL
  long_s <- NULL
  long_m <- NULL
  long_value <- NULL
  for(s_ind in 1){
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

temp_df <- reorg_results(cube_diff(t(mu_MISE[c(1,2,4),])))
scale_breaks <- round(seq(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1)),length.out = 4),2)
scale_breaks <- sort(c(scale_breaks,0))
scale_labels <- as.character(scale_breaks)  
scale_labels[scale_labels == "0"] <- "NeMO"

p1 <- ggplot(temp_df) + geom_boxplot(aes(x = method,y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + 
  theme(text = element_text(size = 24),axis.text.x = element_blank()) + ggtitle(expression(mu))

temp_df <- reorg_results(cube_diff(t(lambda_1_MISE[c(1,2,4),])))
temp_df <- na.omit(temp_df)
scale_breaks <- round(seq(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1)),length.out = 4),2)
scale_breaks <- sort(c(scale_breaks,0))
scale_labels <- as.character(scale_breaks)  
scale_labels[scale_labels == "0"] <- "NeMO"

p2 <- ggplot(temp_df) + geom_boxplot(aes(x = method,y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + 
  theme(text = element_text(size = 24),axis.text.x = element_blank())  + ggtitle(expression(lambda[1]))

temp_df <- reorg_results(cube_diff(t(lambda_2_MISE[c(1,2,4),])))
temp_df <- na.omit(temp_df)
scale_breaks <- round(seq(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1)),length.out = 4),2)
scale_breaks <- sort(c(scale_breaks,0))
scale_labels <- as.character(scale_breaks)  
scale_labels[scale_labels == "0"] <- "NeMO"

p3 <- ggplot(temp_df) + geom_boxplot(aes(x = method,y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + 
  theme(text = element_text(size = 24),axis.text.x = element_blank())  + ggtitle(expression(lambda[2]))

temp_df <- reorg_results(cube_diff(t(lambda_3_MISE[c(1,2,4),])))
temp_df <- na.omit(temp_df)
scale_breaks <- round(seq(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1)),length.out = 4),2)
scale_breaks <- sort(c(scale_breaks,0))
scale_labels <- as.character(scale_breaks)  
scale_labels[scale_labels == "0"] <- "NeMO"

p4 <- ggplot(temp_df) + geom_boxplot(aes(x = method,y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + 
  theme(text = element_text(size = 24),axis.text.x = element_blank())  + ggtitle(expression(lambda[3]))

temp_df <- reorg_results(cube_diff(t(est_MISE[c(1,2,4),])))
scale_breaks <- round(seq(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1)),length.out = 4),2)
scale_breaks <- sort(c(scale_breaks,0))
scale_labels <- as.character(scale_breaks)  
scale_labels[scale_labels == "0"] <- "NeMO"

p5 <- ggplot(temp_df) + geom_boxplot(aes(x = as.factor(sparsity),y = value,color = method)) +
  scale_y_continuous(limits = c(quantile(temp_df$value, c(0)),quantile(temp_df$value, c(1))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + 
  theme(text = element_text(size = 24),axis.text.x = element_blank())  + ggtitle("f")

common_plot <- ggarrange(p1,p2,p3,p4,p5,common.legend = TRUE,ncol = 5, legend = "right")
  annotate_figure(common_plot, left = textGrob("MISE", rot = 90, vjust = 1, gp = gpar(cex = 2)))
