library(ggplot2)
library(tidyr)
library(magrittr)
library(splines)
library(dplyr)
library(gridExtra)
library(Rcpp)

data_clean <- function (data, family = "binomial") 
{
  value = index = NULL
  if (!all(c("id", "index", "value") %in% 
           names(data))) {
    stop("Input dataset must have variables 'id', 'index', and 'value'.")
  }
  if (!identical(data, data %>% arrange(id, index))) {
    message("Data have been sorted by id and index; all output will be in this format")
    data = data %>% arrange(id, index)
  }
  data_rows = data %>% select(id, index, value) %>% mutate(row = row_number()) %>% 
    group_by(id) %>% filter(row_number() == 1 | row_number() == 
                              n()) %>% mutate(index = c("first", "last")) %>% 
    select(-value) %>% spread(index, row) %>% ungroup() %>% 
    mutate(subject = row_number())
  data = data %>% group_by(id) %>% mutate(index_scaled = (index - 
                                                            min(index))/max(index)) %>% ungroup()
  colnames(data_rows) <- c("id", "first_row", "last_row", 
                           "subject")
  I = dim(data_rows)[1]
  return(list(Y = data, I = I, Y_rows = data_rows))
}

bfpca2 <- function(Y, npc = 1, Kt = 8, maxiter = 50, t_min = NULL, t_max = NULL, 
                   print.iter = FALSE, row_obj = NULL, seed = 1988, ...) 
{
  curr_iter = 1
  error = rep(NA, maxiter)
  error[1] = 100
  if (is.null(row_obj)) {
    data = data_clean(Y)
    Y = data$Y
    rows = data$Y_rows
    I = data$I
  }
  else {
    rows = row_obj
    I = dim(rows)[1]
  }
  if (Kt < 3) {
    stop("Kt must be greater than or equal to 3.")
  }
  time = Y$index
  if (any(!(Y$value %in% c(0, 1)))) {
    stop("'binomial' family requires data with binary values of 0 or 1")
  }
  if (is.null(t_min)) {
    t_min = min(time)
  }
  if (is.null(t_max)) {
    t_max = max(time)
  }
  knots = quantile(time, probs = seq(0, 1, length = Kt - 2))[-c(1, 
                                                                Kt - 2)]
  Theta_phi = bs(c(t_min, t_max, time), knots = knots, intercept = TRUE)[-(1:2), 
                                                                         ]
  #set.seed(seed)
  xi = matrix(rnorm(dim(Y)[1]), ncol = 1) * 0.5
  alpha_coefs = matrix(coef(glm(Y$value ~ 0 + Theta_phi, family = "binomial")), 
                       Kt, 1)
  psi_coefs = matrix(rnorm(Kt * npc), Kt, npc) * 0.5
  temp_alpha_coefs = alpha_coefs
  temp_psi_coefs = psi_coefs
  phi_a = list(NA, I)
  phi_b = matrix(0, nrow = Kt * (npc + 1), ncol = I)
  scores = matrix(NA, I, npc)
  while (curr_iter < maxiter && error[curr_iter] > 0.001) {
    if (print.iter) {
      message("current iteration: ", curr_iter)
      message("current error: ", error[curr_iter])
    }
    for (i in 1:I) {
      subject_rows = rows$first_row[i]:rows$last_row[i]
      Yi = Y$value[subject_rows]
      Di = length(Yi)
      Theta_i = Theta_phi[subject_rows, ]
      xi_i = xi[subject_rows, ]
      Theta_i_quad = squareTheta(xi_i, Theta_i)
      mlist = expectedScores(Yi, temp_alpha_coefs, temp_psi_coefs, 
                             Theta_i, Theta_i_quad)
      Ci = mlist$Ci
      mi = mlist$mi
      mm = Ci + tcrossprod(mi)
      xi[subject_rows, 1] = expectedXi(Theta_i, temp_alpha_coefs, 
                                       mi, temp_psi_coefs, Ci)
      xi_i = xi[subject_rows, ]
      Theta_i_quad = squareTheta(xi_i, Theta_i)
      si = rbind(mi, 1)
      ss = cbind(rbind(mm, t(mi)), si)
      phi_a[[i]] = 2 * kronecker(Theta_i_quad, ss)
      phi_b[, i] = t((Yi - 0.5) %*% kronecker(Theta_i, 
                                              t(si)))
      scores[i, ] = mi
    }
    phi_a_sum = Reduce("+", phi_a)
    phi_vec = -solve(phi_a_sum) %*% rowSums(phi_b)
    phi_mat = matrix(phi_vec, nrow = Kt, ncol = npc + 1, 
                     byrow = TRUE)
    alpha_coefs = phi_mat[, npc + 1]
    psi_coefs = phi_mat[, 1:npc]
    if (npc == 1) {
      psi_coefs = matrix(psi_coefs, ncol = 1)
    }
    curr_iter = curr_iter + 1
    error[curr_iter] = sum((psi_coefs - temp_psi_coefs)^2) + 
      sum((alpha_coefs - temp_alpha_coefs)^2)
    temp_psi_coefs = psi_coefs
    temp_alpha_coefs = alpha_coefs
  }
  fits = rep(NA, dim(Y)[1])
  subject_coef = alpha_coefs + tcrossprod(psi_coefs, scores)
  for (i in 1:I) {
    subject_rows = rows$first_row[i]:rows$last_row[i]
    fits[subject_rows] = Theta_phi[subject_rows, ] %*% subject_coef[, 
                                                                    i]
  }
  fittedVals = data.frame(id = Y$id, index = Y$index, value = fits)
  Theta_phi_mean = bs(seq(t_min, t_max, length.out = M), knots = knots, 
                      intercept = TRUE)
  psi_svd = svd(Theta_phi_mean %*% psi_coefs)
  efunctions = psi_svd$u
  evalues = (psi_svd$d)^2
  scores = scores %*% psi_svd$v
  ret = list(knots = knots, alpha = Theta_phi_mean %*% alpha_coefs,
             mu = Theta_phi_mean %*% alpha_coefs, efunctions = efunctions, 
             evalues = evalues, npc = npc, scores = scores, subject_coefs = subject_coef, 
             Yhat = fittedVals, Y = Y, family = "binomial", 
             error = error[!is.na(error)])
  class(ret) = "fpca"
  return(ret)
}