library(MASS)
library(pracma)
library(fdapace)
library(fpca)
library(Rcpp)
library(SemiPar)
library(rootSolve)
library(reshape2)
library(R2WinBUGS)
library(ggpubr)
library(grid)

sourceCpp("../src/nemo_ffa.cpp")
sourceCpp("../src/msf.cpp")
source("../src/nemo_ffa.R")
source("./supp_funs/fpca_funs.R")

# noise setting - either "low" or "high"
noise_setting <- "high"

set.seed(1234)

n_rep <- 100
M <- 30
N <- 100
one_M <- rep(1,M)
one_N <- rep(1,N)

y_long_save <- list()
y_long_count <- 1

mu_MISE <- array(NA,dim = c(3,4,n_rep))
lambda_1_MISE <- array(NA,dim = c(3,4,n_rep))
lambda_2_MISE <- array(NA,dim = c(3,4,n_rep))
est_MISE <- array(NA,dim = c(3,4,n_rep))
K_selected <- array(NA,dim = c(3,4,n_rep))

sparsity_levels <- c(.25,.5,.75)
for(replicate in 1:n_rep){
  ##### Parameters that may need to be varied for simulation experiments #####
  K_true <- 2
  lambda_l <- .4*rep(1,K_true)
  lambda_scale <- 1*rep(1,K_true)
  mu_l <- .4
  mu_scale <- 1
  if(noise_setting == "low"){sig_true <- .5}
  if(noise_setting == "high"){sig_true <- 1}
  
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
  
  nu_eta <- 1e-4
  eta_cov <- diag(one_N) - one_N%*%t(one_N)/(nu_eta + N)
  eta_true <- mvrnorm(K_true,mu = rep(0,N),Sigma = eta_cov)
  eta_sd <- sqrt(apply(eta_true,1,var))
  
  lambda_true <- lambda_true*(one_M%*%t(eta_sd))
  eta_true <- eta_true/(eta_sd%*%t(one_N))
  
  lambda_norm <- sqrt(diag(t(lambda_true)%*%W%*%lambda_true))
  lambda_true <- lambda_true*(one_M%*%t(c(1.5,.75)/lambda_norm))
  
  lambda_norm_order <- order(lambda_norm,decreasing = T)
  lambda_true <- lambda_true[,lambda_norm_order]
  eta_true <- eta_true[lambda_norm_order,]
  
  f <- mu_true%*%t(rep(1,N)) + lambda_true%*%eta_true
  noise <- matrix(rnorm(N*M,0,sig_true),nrow = M,ncol = N)
  y <- f + noise
  
  
  for(sparsity_ind in 1:3){
    sparsity_level <- sparsity_levels[sparsity_ind]
    
    obs_grid <- matrix(sample(c(0,1),M*N,replace = T,prob = c(sparsity_level,1-sparsity_level)),nrow = M,ncol = N)
    while(prod(!apply(obs_grid,2,sum) <= 1) == 0){
      obs_grid <- matrix(sample(c(0,1),M*N,replace = T,prob = c(sparsity_level,1-sparsity_level)),nrow = M,ncol = N)
    }
    
    y <- f + noise
    y_common <- y*obs_grid
    
    y_long <- NULL
    for(i in 1:N){
      id_temp <- rep(i,sum(obs_grid[,i]))
      t_temp <- time[which(obs_grid[,i]>0)]
      y_temp <- y_common[which(obs_grid[,i]>0),i]
      y_long <- rbind(y_long,cbind(id_temp,t_temp,y_temp))
    }
    
    y_long_save[[y_long_count]] <- y_long
    y_long_count <- y_long_count + 1
    
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
    lag <- 1
    
    ffa_mcmc <- run_ffa_mcmc(y_common,obs_grid,time,
                             nu_eta,nu_lambda,
                             init_mcmc_params,
                             prior_a,prior_b,prior_var,
                             n_iter,lag)
    
    K_selected[sparsity_ind,1,replicate] <- dim(init_mcmc_params$eta_cur)[1]
    
    # compute mu MISE
    mu_post_mean <- apply(ffa_mcmc$mu_save,1,mean)
    mu_MISE_nemo <- diag(t(mu_true - mu_post_mean)%*%W%*%(mu_true - mu_post_mean))
    
    # compute lambda MISE
    # resolve label switching an sign ambiguity
    lambda_save_processed <- ffa_mcmc$lambda_save
    for(iter in 1:dim(lambda_save_processed)[3]){
      lambda_post <- ffa_mcmc$lambda_save[,,iter]*(one_M%*%t(sqrt(ffa_mcmc$psi_save[,iter])))
      eta_post <- ffa_mcmc$eta_save[,,iter]/(sqrt(ffa_mcmc$psi_save[,iter])%*%t(one_N))
      eta_sd <- sqrt(apply(eta_post,1,var))
      
      lambda_post <- lambda_post*(one_M%*%t(eta_sd))
      eta_post <- eta_post/(eta_sd%*%t(one_N))
      
      perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
      lambda_save_processed[,,iter] <- aplr(lambda_post,perm_out)
    }
    lambda_post_mean <- apply(lambda_save_processed,c(1,2),mean)
    if(dim(lambda_post_mean)[2]<K_true){
      lambda_post_mean <- cbind(lambda_post_mean,
                                t(repmat(rep(0,M),K_true - dim(lambda_post_mean)[2]) ))
    }
    
    lambda_MISE_nemo <- diag(t(lambda_true - lambda_post_mean)%*%W%*%(lambda_true - lambda_post_mean))
    
    # overall fit
    fit_save <- array(0,dim = c(M,N,n_iter/2))
    for(iter in (n_iter/2 + 1):n_iter){
      if(K_selected[sparsity_ind,1,replicate] == 1){
        fit_save[,,iter - n_iter/2] <- ffa_mcmc$mu_save[,iter]%*%t(one_N) + ffa_mcmc$lambda_save[,,iter]%*%t(ffa_mcmc$eta_save[,,iter])
      }else{
        fit_save[,,iter - n_iter/2] <- ffa_mcmc$mu_save[,iter]%*%t(one_N) + ffa_mcmc$lambda_save[,,iter]%*%ffa_mcmc$eta_save[,,iter]
      }
      
    }
    fit_mean <- apply(fit_save,c(1,2),mean)
    MISE <- diag(t(fit_mean - f)%*%W%*%(fit_mean - f))
    
    # store estimated quantities

    mu_MISE[sparsity_ind,1,replicate] <- mu_MISE_nemo
    lambda_1_MISE[sparsity_ind,1,replicate] <- lambda_MISE_nemo[1]
    lambda_2_MISE[sparsity_ind,1,replicate] <- lambda_MISE_nemo[2]
    
    
    
    ##### Run PACE #####
    
    Ly <- list()
    Lt <- list()
    for(i in 1:N){
      inds_obs <- which(obs_grid[,i] > 0)
      Ly[[i]] <-y_common[inds_obs,i]
      Lt[[i]] <- time[inds_obs]
    }
    
    pace_options <- list()
    pace_options$dataType <- 'Sparse'
    pace_options$maxK <- 5
    
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
    eta_sd <- sqrt(apply(pace_results$xiEst,2,var))
    lambda_post <- phi_obs_grid*(one_M%*%t(eta_sd))
    if(pace_results$selectK > K_true){
      lambda_post <- lambda_post[,1:K_true]
    }
    perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
    lambda_post <- aplr(lambda_post,perm_out)
    lambda_MISE_pace <- diag(t(lambda_true - lambda_post)%*%W%*%(lambda_true - lambda_post))
    
    # overall fit
    fit_pace <- mu_obs_grid%*%t(one_N) + phi_obs_grid%*%t(pace_results$xiEst)
    
    MISE_pace <- diag(t(fit_pace - f)%*%W%*%(fit_pace - f))
    
    # store estimated quantities
    K_selected[sparsity_ind,2,replicate] <- pace_results$selectK
    mu_MISE[sparsity_ind,2,replicate] <- mu_MISE_pace
    lambda_1_MISE[sparsity_ind,2,replicate] <- lambda_MISE_pace[1]
    lambda_2_MISE[sparsity_ind,2,replicate] <- lambda_MISE_pace[2]
    
    ##### Peng and Paul FPCA #####
    fpca_data <- cbind(y_long[,1],y_long[,3],y_long[,2]/2)
    fpca_results <- fpca.mle(fpca_data,c(5,10,15,20),2:5)
    fpca_scores <- fpca.score2(fpca_data,fpca_results$grid,fpca_results$fitted_mean,
                               fpca_results$eigenvalues,fpca_results$eigenfunctions,fpca_results$error_var,fpca_results$selected_model[2])
    mu_dense <- fpca_results$fitted_mean
    lambda_dense <- fpca_results$eigenfunctions
    
    mu_og_grid <- approx(x = fpca_results$grid,mu_dense,xout = time/2)$y
    lambda_og_grid <- array(0,dim = c(length(mu_og_grid),fpca_results$selected_model[2]))
    for(k in 1:fpca_results$selected_model[2]){
      lambda_og_grid[,k] <- approx(x = fpca_results$grid,lambda_dense[k,],xout = time/2)$y
    }
    
    fit_fpca_dense <- fpca_results$fitted_mean%*%t(rep(1,N)) + t(fpca_results$eigenfunctions)%*%t(fpca_scores)
    
    fit_fpca_og <- array(0,dim = c(M,N))
    for(i in 1:N){
      fit_fpca_og[,i] <- approx(x = fpca_results$grid,fit_fpca_dense[,i],xout = time/2)$y
    }
    
    MISE_fpca <- diag(t(fit_fpca_og - f)%*%W%*%(fit_fpca_og - f))

    # compute mu MISE
    mu_MISE_fpca <- diag(t(mu_true - mu_og_grid)%*%W%*%(mu_true - mu_og_grid))
    
    # compute lambda MISE
    eta_sd <- sqrt(apply(fpca_scores,2,var))
    lambda_post <- lambda_og_grid*(one_M%*%t(eta_sd))
    if(fpca_results$selected_model[2] > K_true){
      lambda_post <- lambda_post[,1:K_true]
    }
    perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
    lambda_post <- aplr(lambda_post,perm_out)
    lambda_MISE_fpca <- diag(t(lambda_true - lambda_post)%*%W%*%(lambda_true - lambda_post))
    
    # store estimated quantities
    K_selected[sparsity_ind,3,replicate] <- fpca_results$selected_model[2]
    mu_MISE[sparsity_ind,3,replicate] <- mu_MISE_fpca
    lambda_1_MISE[sparsity_ind,3,replicate] <- lambda_MISE_fpca[1]
    lambda_2_MISE[sparsity_ind,3,replicate] <- lambda_MISE_fpca[2]
    
    ##### Run GC model using WinBUGS #####
    
    y <- as.matrix(dcast(as.data.frame(y_long),id_temp~t_temp,value.var = "y_temp"))
    y <- y[,2:31]
    
    # set variables for WinBUGS
    N_subj=dim(y)[1]
    N_obs=dim(y)[2]
    dim.space=2
    
    y_centered=y
    y_centered<-y_centered-matrix(rep(colMeans(y_centered,na.rm=TRUE),nrow(y_centered)),nrow=nrow(y_centered),byrow=TRUE)
    
    # Estimate covariance function through bivaraite smoothing
    Gt=cov(y_centered,use="pairwise.complete.obs")
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
    
    perm_out <- msfOUT_functional(E,W,lambda_true)
    E <- aplr(E,perm_out)
    
    ##### MCMC Using WinBUGS #####
    data<-list("E","y_centered","N_subj","N_obs","dim.space")
    
    #This is the name of the program
    program.file.name="wb_model.txt"
    
    #Define the initial values & parameters to record 
    inits.y_centered=matrix(rep(NA,N_subj*N_obs),ncol=N_obs)
    inits.ll=rep(0.01,dim.space)
    inits.y_centered[is.na(y_centered)]=mean(mean(y_centered,na.rm=TRUE))
    inits<-function(){list(xi=matrix(rep(0,N_subj*dim.space),ncol=dim.space),
                           taueps=0.01,ll=inits.ll,y_centered=inits.y_centered)}
    
    parameters=list("lambda","xi")
    
    # MCMC
    n.thin=1
    n.iter=1500
    n.burnin=500

    Bayes.fit<- bugs(data, inits, parameters, model.file = program.file.name,
                     n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
                     n.thin = n.thin,debug = FALSE, DIC = FALSE, digits = 5, codaPkg = FALSE,
                     bugs.directory = "c:/Program Files/WinBUGS14/")
    
    # compute mu MISE
    y_mean=colMeans(y[1:N_subj,1:N_obs],na.rm=TRUE)
    mu_MISE_cg <- diag(t(mu_true - y_mean)%*%W%*%(mu_true - y_mean))
    
    # compute lambda MISE
    eta_post_mean <- Bayes.fit$mean$xi
    eta_sd <- sqrt(apply(eta_post_mean,2,var))
    lambda_post <- E*(one_M%*%t(eta_sd))
    perm_out <- msfOUT_functional(lambda_post,W,lambda_true)
    lambda_post <- aplr(lambda_post,perm_out)
    lambda_MISE_cg <- diag(t(lambda_true - lambda_post)%*%W%*%(lambda_true - lambda_post))
    
    # store estimated quantities
    K_selected[sparsity_ind,4,replicate] <- 2
    mu_MISE[sparsity_ind,4,replicate] <- mu_MISE_cg
    lambda_1_MISE[sparsity_ind,4,replicate] <- lambda_MISE_cg[1]
    lambda_2_MISE[sparsity_ind,4,replicate] <- lambda_MISE_cg[2]
    
    #overall fit
    wb_fit <- t(Bayes.fit$mean$xi%*%t(E) + rep(1,dim(Bayes.fit$mean$xi)[1])%*%t(y_mean))
    wb_MISE <- diag(t(wb_fit - f)%*%W%*%(wb_fit - f))
    
    est_MISE[sparsity_ind,1,replicate] <- mean(MISE)
    est_MISE[sparsity_ind,2,replicate] <- mean(MISE_pace)
    est_MISE[sparsity_ind,3,replicate] <- mean(MISE_fpca)
    est_MISE[sparsity_ind,4,replicate] <- mean(wb_MISE)
    
  }
  
  print(K_selected[,,replicate])
  print(mu_MISE[,,replicate])
  print(lambda_1_MISE[,,replicate])
  print(lambda_2_MISE[,,replicate])
  print(est_MISE[,,replicate])
}



##### Visualize results #####


# accidentally saved incorrectly, so resolving manuyally
pace_MISE <- est_MISE[sparsity_ind,3,replicate]
fpca_MISE <- est_MISE[sparsity_ind,2,replicate]

est_MISE[sparsity_ind,3,replicate] <- fpca_MISE
est_MISE[sparsity_ind,3,replicate] <- pace_MISE

cube_diff <- function(cube_array){
diff_array <- array(0,dim= c(3,3,100))
for(method in 1:3){
  diff_array[,method,] <- cube_array[,1+method,] - cube_array[,1,]
}
return(diff_array)
}

reorg_results <- function(cube_array){
s_levels <- c(.25,.5,.75)
m_levels <- c("YMW2005","PP2009","CG2010")
long_rep <- NULL
long_s <- NULL
long_m <- NULL
long_value <- NULL
for(s_ind in 1:3){
  for(m_ind in 1:3){
    long_rep <- c(long_rep,1:100)
    long_s <- c(long_s,rep(s_levels[s_ind],100))
    long_m <- c(long_m,rep(m_levels[m_ind],100))
    long_value <- c(long_value,cube_array[s_ind,m_ind,])
  }
}
long_df <- data.frame(replicate = long_rep,sparsity = long_s,method = long_m,value = long_value)
return(long_df)
}

temp_df <- reorg_results(cube_diff(mu_MISE))
scale_breaks <- round(seq(0,quantile(temp_df$value, c(0.95)),length.out = 4),2)
scale_labels <- as.character(scale_breaks)  
scale_labels[1] <- "NeMO"

p1 <- ggplot(temp_df) + geom_boxplot(aes(x = as.factor(sparsity),y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(0,quantile(temp_df$value, c(0.95))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + theme(text = element_text(size = 24)) + ggtitle(expression(mu))

temp_df <- reorg_results(cube_diff(lambda_1_MISE))
scale_breaks <- round(seq(0,quantile(temp_df$value, c(0.95)),length.out = 4),2)
scale_labels <- as.character(scale_breaks)  
scale_labels[1] <- "NeMO"

p2 <- ggplot(temp_df) + geom_boxplot(aes(x = as.factor(sparsity),y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(0,quantile(temp_df$value, c(0.95))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + theme(text = element_text(size = 24)) + ggtitle(expression(lambda[1]))

temp_df <- reorg_results(cube_diff(lambda_2_MISE))
scale_breaks <- round(seq(0,quantile(temp_df$value, c(0.95)),length.out = 4),2)
scale_labels <- as.character(scale_breaks)  
scale_labels[1] <- "NeMO"

p3 <- ggplot(temp_df) + geom_boxplot(aes(x = as.factor(sparsity),y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(0,quantile(temp_df$value, c(0.95))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + theme(text = element_text(size = 24)) + ggtitle(expression(lambda[2]))

temp_df <- reorg_results(cube_diff(est_MISE))
scale_breaks <- round(seq(0,quantile(temp_df$value, c(0.95)),length.out = 4),2)
scale_labels <- as.character(scale_breaks)  
scale_labels[1] <- "NeMO"

p4 <- ggplot(temp_df) + geom_boxplot(aes(x = as.factor(sparsity),y = value,color = method),outlier.shape = NA) +
  scale_y_continuous(limits = c(0,quantile(temp_df$value, c(0.95))),breaks = scale_breaks,labels = scale_labels) +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1)+ ylab("") + xlab("") + theme(text = element_text(size = 24)) + ggtitle("f")

common_plot <- ggarrange(p1,p2,p3,p4,common.legend = TRUE,ncol = 4, legend = "right")
annotate_figure(common_plot, left = textGrob("MISE", rot = 90, vjust = 1, gp = gpar(cex = 2)),
                bottom = textGrob("sparsity level", gp = gpar(cex = 2)))