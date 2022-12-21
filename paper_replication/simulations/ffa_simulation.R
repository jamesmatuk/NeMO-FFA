library(MASS)
library(pracma)
library(fdapace)
library(fpca)
library(Rcpp)
library(SemiPar)


sourceCpp("../src/nemo.cpp")

set.seed(1234)

n_rep <- 5
M <- 30
y_long_save <- list()
y_long_count <- 1

MISE_mat <- array(NA,dim = c(3,4,n_rep))
sparsity_levels <- c(.25,.5,.75)
for(replicate in 1:n_rep){
  ##### Parameters that may need to be varied for simulation experiments #####
  
  N <- 100
  K_true <- 2
  lambda_l <- .4*rep(1,K_true)
  lambda_scale <- 1*rep(1,K_true)
  mu_l <- .4
  mu_scale <- 1
  sig_true <- .5 # low noise setting
  # sig_true <- 1 # high noise setting
  
  
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
    w <- diff(c(time[1],(time[2:M]+time[1:(M-1)])/2,time[M]))
    W <- diag(w)
    one_vec_N <- rep(1,N)
    K <- 5
    nu_eta <- 1e-4
    nu_param <- 1e-4
    
    lambda_cur <- matrix(rnorm(M*K),nrow = M,ncol = K)
    eta_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
    psi_cur <- rep(1,K)
    sig_sq_cur <- var((y-lambda_cur%*%eta_cur)[obs_grid>0])
    l_mu_current <- .1
    scale_mu_current <- 1
    l_param_cur <- rep(.1,K)
    scale_param_cur <- rep(1,K)
    
    prior_a <- 12
    prior_b <- 4
    prior_var <- 1
    
    px_b <- 1
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
    n_iter <- 1000
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
    
    fit_save <- array(0,dim = c(M,N,n_iter/2))
    for(iter in (n_iter/2 + 1):n_iter){
      fit_save[,,iter - n_iter/2] <- mu_save[,iter]%*%t(one_vec_N) + lambda_save[,,iter]%*%eta_save[,,iter]
    }
    fit_mean <- apply(fit_save,c(1,2),mean)
    MISE <- diag(t(fit_mean - f)%*%W%*%(fit_mean - f))
    
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
    pace_options$plot <- TRUE
    pace_options$methodSelectK <- 5
    
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
    
    fit_pace <- mu_obs_grid%*%t(one_vec_N) + phi_obs_grid%*%t(pace_results$xiEst)
    
    MISE_pace <- diag(t(fit_pace - f)%*%W%*%(fit_pace - f))
    
    
    ##### Peng and Paul FPCA #####
    tryCatch({ # There  sometimes when this code throws an error
      fpca_data <- cbind(y_long[,1],y_long[,3],y_long[,2]/2)
      fpca_results <- fpca.mle(fpca_data,c(5,10,15,20),5)
      
      fpca.score2 <- function(data.m, grids.u, muhat, eigenvals, eigenfuncs, sig2hat,K) 
      {
        temp <- table(data.m[, 1])
        n <- length(temp)
        m.l <- as.vector(temp)
        result <- matrix(0, n, K)
        N <- length(grids.u)
        evalmat <- diag(eigenvals[1:K])
        current <- 0
        eigenfuncs.u <- t(eigenfuncs)
        data.u <- matrix(as.numeric(as.vector(data.m[, -1])), nrow = nrow(data.m[, 
                                                                                 -1]), ncol = ncol(data.m[, -1]))
        for (i in 1:n) {
          Y <- as.vector(data.u[(current + 1):(current + m.l[i]), 
                                1])
          meastime <- data.u[(current + 1):(current + m.l[i]), 
                             2]
          gridtime <- ceiling(N * meastime) 
          gridtime[gridtime == 0] <- 1 # there is an off by 1 error in the R package
          muy <- muhat[gridtime]
          Phiy <- matrix(eigenfuncs.u[gridtime, 1:K], ncol = K)
          Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
          temp.y <- matrix(Y - muy)
          result[i, ] <- evalmat %*% t(Phiy) %*% solve(Sigy, temp.y)
          current <- current + m.l[i]
        }
        return(result)
      }
      
      fpca_scores <- fpca.score2(fpca_data,fpca_results$grid,fpca_results$fitted_mean,
                                 fpca_results$eigenvalues,fpca_results$eigenfunctions,fpca_results$error_var,5)
      fit_fpca_dense <- fpca_results$fitted_mean%*%t(rep(1,N)) + t(fpca_results$eigenfunctions)%*%t(fpca_scores)
      
      fit_fpca_og <- array(0,dim = c(M,N))
      for(i in 1:N){
        fit_fpca_og[,i] <- approx(x = fpca_results$grid,fit_fpca_dense[,i],xout = time/2)$y
      }
      
      MISE_fpca <- diag(t(fit_fpca_og - f)%*%W%*%(fit_fpca_og - f))
      MISE_mat[sparsity_ind,2,replicate] <- mean(MISE_fpca)
    },error = function(cond){
      MISE_mat[sparsity_ind,2,replicate] <- NA
    })
    
    ##### Run GC model using WinBUGS #####
    
    time_grid<- time
    grid_weight <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
    grid_weight <- diag(grid_weight)
    
    y <- as.matrix(dcast(as.data.frame(y_long),id_temp~t_temp,value.var = "y_temp"))
    y <- y[,2:31]
    
    # set variables for WinBUGS
    N_subj=dim(y)[1]
    N_obs=dim(y)[2]
    dim.space=5
    
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
    program.file.name="wb_model.txt"
    
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
    
    
    library(R2WinBUGS)
    Bayes.fit<- bugs(data, inits, parameters, model.file = program.file.name,
                     n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
                     n.thin = n.thin,debug = FALSE, DIC = FALSE, digits = 5, codaPkg = FALSE,
                     bugs.directory = "c:/Program Files/WinBUGS14/")
    
    
    y_mean=colMeans(y[1:N_subj,1:N_obs],na.rm=TRUE)
    wb_fit <- t(Bayes.fit$mean$xi%*%t(E) + rep(1,dim(Bayes.fit$mean$xi)[1])%*%t(y_mean))
    
    wb_MISE <- diag(t(wb_fit - f)%*%grid_weight%*%(wb_fit - f))
    
    MISE_mat[sparsity_ind,1,replicate] <- mean(MISE)
    MISE_mat[sparsity_ind,3,replicate] <- mean(MISE_pace)
    MISE_mat[sparsity_ind,4,replicate] <- mean(wb_MISE)
    
  }
  
  print(replicate)
  print(MISE_mat[,,replicate])
  
}

##### Visualize results #####

MISE_mean <- apply(MISE_mat,c(1,2),mean)
MISE_sd <- apply(MISE_mat,c(1,2),sd)
MISE_lower<- MISE_mean - 1.96*MISE_sd/sqrt(dim(MISE_mat)[3])
MISE_upper<- MISE_mean + 1.96*MISE_sd/sqrt(dim(MISE_mat)[3])

method_fact <- factor(c(rep("ReMO",3),rep("Newton",3),rep("PACE",3),rep("CG",3)),levels = c("ReMO","Newton","PACE","CG"))
mean_df = data.frame(method = method_fact,value = c(MISE_mean))
lower_df = data.frame(method = method_fact,value = c(MISE_lower))
upper_df = data.frame(method = method_fact,value = c(MISE_upper))

sparsity_levels <- c(.25,.5,.75)
yl <- c(min(MISE_lower),max(MISE_upper))
ggplot(mean_df) + geom_line(aes(x = rep(sparsity_levels,4),y = value,color = method),size = 2) + 
  geom_ribbon(mapping = aes(x = rep(sparsity_levels,4),ymin = lower_df$value,ymax = upper_df$value ,fill = method),alpha = .2) + 
  ylab("MISE") + xlab("sparsity level") + ylim(yl) + theme(text = element_text(size = 24)) +
  scale_x_continuous(breaks = sparsity_levels)