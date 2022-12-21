library(MASS)
library(pracma)
library(fdapace)
library(fpca)
library(Rcpp)
library(R2WinBUGS)

sourceCpp("../src/nemo_ffa.cpp")

set.seed(1234)

n_rep <- 100
M <- 30
N <- 100
f_save <- array(0,dim = c(M,N,n_rep))
y_long_save <- list()
y_long_count <- 1
lambda_true_save <- array(0,dim = c(M,n_rep))
x_save <- array(0,dim = c(N,n_rep))
b_true_save <- array(0,dim = c(n_rep))


covered_mat <- array(NA,dim = c(3,4,n_rep))
sparsity_levels <- c(.25,.5,.75)
for(replicate in 1:n_rep){
  ##### Parameters that may need to be varied for simulation experiments #####
  
  N <- 100
  K_true <- 1
  lambda_l <- .4*rep(1,K_true)
  lambda_scale <- 1*rep(1,K_true)
  mu_l <- .4
  mu_scale <- 1
  sig_true <- .5 # low noise
  #sig_true <- 1 # high noise
  
  ##### generate underlying observations #####
  
  M <- 30 
  time <- seq(0,2,length.out = M)
  w <- diff(c(time[1],(time[2:M]+time[1:(M-1)])/2,time[M]))
  W <- diag(w)
  
  mu_cov <- make_cov(time,mu_l,mu_scale) 
  mu_true <-  mvrnorm(mu = rep(0,length(time)),Sigma = mu_cov)
  
  prior_cov <- make_cov(time,mu_l,mu_scale) 
  lambda_true <- as.matrix(mvrnorm(1,mu = rep(0,length(time)),Sigma = prior_cov))
  lambda_norm <- sqrt(t(lambda_true)%*%W%*%lambda_true)
  lambda_true <- lambda_true/as.numeric(lambda_norm)
  
  xi_true <- matrix(rnorm(K_true*N),nrow = K_true,ncol = N)
  x <- rnorm(N)
  b_true <- rnorm(1)
  eta_true <- t(x)*b_true + xi_true
  
  f <- mu_true%*%t(rep(1,N)) + lambda_true%*%eta_true
  noise <- matrix(rnorm(N*M,0,sig_true),nrow = M,ncol = N)
  y <- f + noise
  
  f_save[,,replicate] <- f
  lambda_true_save[,replicate] <- lambda_true
  x_save[,replicate] <- x
  b_true_save[replicate] <- b_true
  
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
    
    ##### Run MCMC #####

    w <- diff(c(time[1],(time[2:M]+time[1:(M-1)])/2,time[M]))
    W <- diag(w)
    
    N <- dim(y_common)[2]
    K <- 2
    Q <- 1
    one_vec_N <- rep(1,N)
    nu_eta <- 1e-4
    nu_param <- 1e-4
    
    prior_a <- 12
    prior_b <- 4
    prior_var <- 1
    
    age_grid <- time
    M <- length(age_grid)
    w <- diff(c(age_grid[1],(age_grid[2:M]+age_grid[1:(M-1)])/2,age_grid[M]))
    W <- diag(w)
    
    #lambda_cur <- matrix(rnorm(M*K),nrow = M,ncol = K)
    lambda_cur <- matrix(0,nrow = M,ncol = K)
    xi_cur <- matrix(rnorm(K*N,0,1),nrow = K,ncol = N)
    b_cur <- matrix(rnorm(K*Q,0,1),nrow = K,ncol = Q)
    eta_cur <- b_cur%*%x + xi_cur
    psi_cur <- rep(1,K)
    sig_sq_cur <- var((y_common-lambda_cur%*%eta_cur)[obs_grid>0])
    l_mu_current <- 1
    scale_mu_current <- 1
    l_param_cur <- rep(1,K)
    scale_param_cur <- rep(1,K)
    
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
    n_iter <- 20000
    lag <- 20
    n_save <- n_iter/lag
    mu_save <- matrix(0,nrow = M,ncol = n_save)
    scale_mu_save <- rep(0,n_save)
    l_mu_save <- rep(0,n_save)
    lambda_save <- array(0,dim = c(M,K,n_save))
    xi_save <- array(0,dim = c(K,N,n_save))
    b_save <- array(0,dim = c(K,Q,n_save))
    eta_save <- array(0,dim = c(K,N,n_save))
    psi_save <- array(0,dim = c(K,n_save))
    l_param_save <- array(0,dim = c(K,n_save))
    scale_param_save <- array(0,dim = c(K,n_save))
    sig_sq_save <- rep(0,n_save)
    
    for(iter in 1:n_iter){
      
      y_st_cur <- compute_y_st(y_common,lambda_cur%*%eta_cur,obs_grid)
      mu_cur <- sample_mu_sparse(y_st_cur, sig_sq_cur, age_grid, l_mu_current,scale_mu_current,obs_grid)
      
      out <- l_mu_MH(y_st_cur,age_grid,obs_grid,l_mu_current,scale_mu_current,mu_cur,sig_sq_cur,proposal_l_mu,1,iter,prior_a,prior_b) 
      l_mu_current <- out[1]
      accept_l_mu_count <- accept_l_mu_count + out[2]
      proposal_l_mu <- out[3]
      
      out <- scale_mu_MH(y_st_cur, age_grid,obs_grid,l_mu_current, scale_mu_current,mu_cur,sig_sq_cur,proposal_scale_mu,1,iter,prior_var)
      scale_mu_current<- out[1]
      accept_scale_mu_count <- accept_scale_mu_count + out[2]
      proposal_scale_mu<- out[3]
      
      y_st_cur <- compute_y_st(y_common,mu_cur%*%t(one_vec_N),obs_grid)
      lambda_cur <- sample_lambda_sparse(y_st_cur,W,age_grid,l_param_cur,scale_param_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,obs_grid)
      
      out <- l_param_MH(l_param_cur,proposal_l,scale_param_cur,age_grid,obs_grid,W,y_st_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,1,iter,prior_a,prior_b)
      l_param_cur <- out[,1]
      accept_l_count <- accept_l_count + out[,2]
      proposal_l <- out[,3]
      
      out <- scale_param_MH(scale_param_cur,proposal_scale,l_param_cur,age_grid,obs_grid,W,y_st_cur,lambda_cur,eta_cur,sig_sq_cur,nu_param,1,iter,prior_var)
      scale_param_cur <- out[,1]
      accept_scale_count <- accept_scale_count + out[,2]
      proposal_scale <- out[,3]
      
      y_st_cur <- compute_y_st(y_common,mu_cur%*%t(one_vec_N) + lambda_cur%*%b_cur%*%x,obs_grid)
      xi_cur <- sample_eta_cholupdate(y_st_cur,xi_cur,psi_cur,lambda_cur,sig_sq_cur,nu_eta,obs_grid)
      
      y_st_cur <- compute_y_st(y_common,mu_cur%*%t(one_vec_N) + lambda_cur%*%xi_cur,obs_grid)  
      b_cur <- sample_b(y_st_cur,t(as.matrix(x)),b_cur,lambda_cur,psi_cur,sig_sq_cur,obs_grid)
      
      eta_cur <- b_cur%*%x + xi_cur
      
      psi_cur <- sample_psi_regression(xi_cur,b_cur,nu_eta,11,10)
      
      y_st_cur <- compute_y_st(y_common,mu_cur%*%t(one_vec_N),obs_grid)
      sig_sq_cur <- sample_sig_sq_sparse(y_st_cur,lambda_cur,eta_cur,obs_grid)
      
      if(iter%%lag == 0){
        iter_save <- iter/lag
        mu_save[,iter_save] <- mu_cur
        scale_mu_save[iter_save] <- scale_mu_current
        l_mu_save[iter_save] <- l_mu_current
        lambda_save[,,iter_save] <- lambda_cur
        xi_save[,,iter_save] <- xi_cur
        b_save[,,iter_save] <- b_cur
        eta_save[,,iter_save] <- eta_cur
        psi_save[,iter_save] <- psi_cur
        l_param_save[,iter_save] <- l_param_cur
        scale_param_save[,iter_save] <- scale_param_cur
        sig_sq_save[iter_save] <- sig_sq_cur  
      }
      
    }
    
    # post-process mcmc samples
    # need to resolve order, sign, and scale ambiguity in each mcmc sample
    b_processed <- rep(0,n_save)
    for(iter in 1:n_save){
      lambda_norm <- sqrt(diag(t(lambda_save[,,iter])%*%W%*%lambda_save[,,iter]))
      
      lambda_temp <- lambda_save[,order(lambda_norm,decreasing = T)[1],iter]
      lambda_temp_norm <- sqrt(t(lambda_temp)%*%W%*%lambda_temp)
      lambda_temp <- lambda_temp/as.numeric(lambda_temp_norm)
      
      pos_neg <- c(1,-1)
      sign_est <- pos_neg[order(c(sum((lambda_true - lambda_temp)^2),sum((lambda_true + lambda_temp)^2)))[1]]
      lambda_temp <- lambda_temp*sign_est
      
      b_processed[iter]  <- sign_est*as.numeric(lambda_temp_norm)*b_save[order(lambda_norm,decreasing = T)[1],1,iter]
    }
    ortho_ci <- as.numeric(quantile(b_processed[(n_save/2):n_save],probs = c(.025,.975)))
    
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
    pace_options$plot <- FALSE
    pace_options$methodSelectK <- 2
    
    pace_results <- FPCA(Ly,Lt,optns = pace_options)
    
    phi_obs_grid <- ConvertSupport(
      pace_results$workGrid,
      time,
      phi = pace_results$phi
    )
    
    phi_norm <- sqrt(t(phi_obs_grid[,1])%*%W%*%phi_obs_grid[,1])
    phi_est <- phi_obs_grid[,1]/as.numeric(phi_norm)
    
    pos_neg <- c(1,-1)
    sign_est <- pos_neg[order(c(sum((lambda_true - phi_est)^2),sum((lambda_true + phi_est)^2)))[1]]
    phi_est <- phi_est*sign_est
    
    xi_est <- sign_est*as.numeric(phi_norm)*pace_results$xiEst[,1]
    
    pace_lm <- lm(xi_est ~ x + 0)
    pace_ci <- as.numeric(confint(pace_lm))
    
    #plot(xi_est,eta_true)
    #print(b_true)
    #print(pace_ci)
    
    ##### Peng and Paul FPCA #####
    tryCatch({ # There are sometimes when this code throws an error
      fpca_data <- cbind(y_long[,1],y_long[,3],y_long[,2]/2)
      fpca_results <- fpca.mle(fpca_data,c(5,10,15,20),2)
      
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
                                 fpca_results$eigenvalues,fpca_results$eigenfunctions,fpca_results$error_var,2)
      
      fpca_efun_og <- approx(x = fpca_results$grid,t(fpca_results$eigenfunctions)[,1],xout = time/2)$y
      fpca_efun_norm <- sqrt(t(fpca_efun_og)%*%W%*%fpca_efun_og)
      fpca_efun <- fpca_efun_og/as.numeric(fpca_efun_norm)
      
      pos_neg <- c(1,-1)
      sign_est <- pos_neg[order(c(sum((lambda_true - fpca_efun)^2),sum((lambda_true + fpca_efun)^2)))[1]]
      fpca_efun <- fpca_efun*sign_est
      
      fpca_score <- fpca_scores[,1]*sign_est*as.numeric(fpca_efun_norm)
      
      fpca_lm <- lm(fpca_score ~ x + 0)
      fpca_ci <- as.numeric(confint(fpca_lm))
      covered_mat[sparsity_ind,2,replicate] <- b_true > fpca_ci[1] & b_true < fpca_ci[2]
    },error = function(cond){covered_mat[sparsity_ind,2,replicate] <- NA})
    
    
    time_grid <- sort(unique(y_long[,2]))
    M <- length(time_grid)
    grid_weight <- diff(c(time_grid[1],(time_grid[2:M]+time_grid[1:(M-1)])/2,time_grid[M]))
    grid_weight <- diag(grid_weight)
    
    y <- as.matrix(dcast(as.data.frame(y_long),id_temp~t_temp,value.var = "y_temp"))
    y <- y[,2:31]
    
    # set variables for WinBUGS
    N_subj=dim(y)[1]
    N_obs=dim(y)[2]
    dim.space=2
    
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
    
    pos_neg <- c(1,-1)
    sign_est <- pos_neg[order(c(sum((lambda_true - E[,1])^2),sum((lambda_true + E[,1])^2)))[1]]
    E[,1] <- E[,1]*sign_est
    
    ##### MCMC Using WinBUGS #####
    data<-list("E","W","x","N_subj","N_obs","dim.space")
    
    #This is the name of the program
    program.file.name="wb_reg_model.txt"
    
    #Define the initial values & parameters to record 
    inits.W=matrix(rep(NA,N_subj*N_obs),ncol=N_obs)
    inits.ll_b=rep(0.01,dim.space)
    inits.ll_xi=rep(0.01,dim.space)
    inits.W[is.na(W)]=mean(mean(W,na.rm=TRUE))
    inits<-function(){list(xi=matrix(rep(0,N_subj*dim.space),ncol=dim.space),
                           beta=rep(0,dim.space),
                           taueps=0.01,ll_b=inits.ll_b,ll_xi=inits.ll_xi,W=inits.W)}
    
    parameters=list("beta","xi")
    
    # MCMC
    n.thin=1
    n.iter=1500
    n.burnin=500
    
    Bayes.fit<- bugs(data, inits, parameters, model.file = program.file.name,
                     n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
                     n.thin = n.thin,debug = FALSE, DIC = FALSE, digits = 5, codaPkg = FALSE,
                     bugs.directory = "c:/Program Files/WinBUGS14/")
    
    
    wb_ci <- quantile(Bayes.fit$sims.list$beta[,1],probs = c(.025,.975))
    
    

    covered_mat[sparsity_ind,1,replicate] <- b_true > ortho_ci[1] & b_true < ortho_ci[2]
    covered_mat[sparsity_ind,3,replicate] <- b_true > pace_ci[1] & b_true < pace_ci[2]
    covered_mat[sparsity_ind,4,replicate] <- b_true > wb_ci[1] & b_true <wb_ci[2]
      
  }
  
  print(replicate)
  print(apply(covered_mat[,,1:replicate],c(1,2),mean))
  
}

##### visualize simulation results #####

covered_proportion <- apply(covered_mat,c(1,2),mean)

method_fact <- factor(c(rep("NeMO",3),rep("Newton",3),rep("PACE",3),rep("CG",3)),levels = c("NeMO","Newton","PACE","CG"))
mean_df = data.frame(method = method_fact,value = c(covered_proportion))

sparsity_levels <- c(.25,.5,.75)
yl <- c(min(mean_df$value),1)
ggplot(mean_df) + geom_line(aes(x = rep(sparsity_levels,4),y = value,color = method),size = 2) + 
  geom_hline(yintercept = .95,linetype = 'dashed') + 
  ylab("coverage") + xlab("sparsity level") + ylim(yl) + theme(text = element_text(size = 24))+
  scale_x_continuous(breaks = sparsity_levels)
