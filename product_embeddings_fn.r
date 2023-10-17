

max_loglik_mixture_similarity <- function(log_prior_similarities, log_estimated_similarities)
{
  loglik_mixture_similarity <- function(
    log_prior_similarities, log_estimated_similarities,
    p_max_related,
    mu_unrelated, sigma_unrelated, 
    mu_related, sigma_related)
  {
    p_related <- p_max_related * exp(log_prior_similarities)
    log_p_mix <- log( p_related )
    log_pm1_mix <- log( 1-p_related )
    
    loglik_unrelated <- dnorm(log_estimated_similarities, mean = mu_unrelated, sd = sigma_unrelated, log = T)
    loglik_related <- dnorm(log_estimated_similarities, mean = mu_related, sd = sigma_related, log = T)
    
    loglik_rowwise <- cbind(log_p_mix+loglik_related, log_pm1_mix+loglik_unrelated)
    matrixStats::rowLogSumExps(loglik_rowwise) %>% sum()
  }
  
  transform_par <- function(par) {
    par[["p_max_related"]] %<>% plogis()
    par[["sigma_unrelated"]] %<>% exp()
    par[["sigma_related"]] %<>% exp()
    par
  }
  
  loglik_mixture_similarity_optim <- function(par, 
                                              log_prior_similarities, 
                                              log_estimated_similarities)
  {
    par %<>% transform_par()
    
    loglik <-  
      loglik_mixture_similarity(log_prior_similarities = log_prior_similarities,
                                log_estimated_similarities = log_estimated_similarities,
                                p_max_related = par[["p_max_related"]],
                                mu_unrelated = par[["mu_unrelated"]], 
                                sigma_unrelated = par[["sigma_unrelated"]],
                                mu_related = par[["mu_related"]], 
                                sigma_related = par[["sigma_related"]])
    
    #print(loglik)
    loglik
  }
  
  sample_median = median(log_estimated_similarities)
  sample_sd = sd(log_estimated_similarities)
  par_start <- c(p_max_related = qlogis(.1),
                 mu_unrelated = sample_median, sigma_unrelated = log(sample_sd), 
                 mu_related = sample_median/2, sigma_related = log(sample_sd))
  
  res <-
    optim(par = par_start, loglik_mixture_similarity_optim, 
          log_prior_similarities = log_prior_similarities, 
          log_estimated_similarities = log_estimated_similarities, 
          control = list(fnscale=-1, maxit = 10000) )
  
  res$par %<>% transform_par()
  
  res
}

corr_posterior <- function(y1, y2, n_samples = 1000)
{
  # samples from the posterior of the poisson rates (lambda) associated with the vector x
  lambda_posterior_sample_fn <- function(x) {
    # set prior for all elements of x based on the average of x
    beta = 0.05
    alpha = mean(x)*beta
    
    # determine parameters of conjugate posteriors for all elements of x
    alpha_prime = x + alpha
    beta_prime = beta + 1
    
    # return function to sample from posterior 
    function(n) rgamma(n=length(x), shape = alpha_prime, scale = 1/beta_prime)
  }
  
  # initialize functions to sample from posteriors associated with y1 and y2
  y1_lambda_posterior_sample <- lambda_posterior_sample_fn(y1)
  y2_lambda_posterior_sample <- lambda_posterior_sample_fn(y2)
  
  # obtain n samples from the y1 and y2 posteriors, and compute correlation 
  corr_posterior_samples <- sapply(1:n_samples, function(i) {
    cor(y1_lambda_posterior_sample(), y2_lambda_posterior_sample()) 
  })
  
  # return 95% CI
  quantile(corr_posterior_samples, c(0.025, 0.975))
}


corr_posterior_mat <- function(seasonality_mat)
{
  n_products <- nrow(seasonality_mat)
  products <- rownames(seasonality_mat)

  corr_mat <- map(1:n_products, function(i)
  {
      product_1 <- products[i] %>% as.character()
      map_dbl(1:n_products, function(j)
      {
          product_2 <- products[j] %>% as.character()
          if (i < j) { NA } 
          else if (i == j) { 1 } 
          else { 
              series_1 <- seasonality_mat[as.character(product_1),]
              series_2 <- seasonality_mat[as.character(product_2),]
              corr_posterior(series_1, series_2)[2]
          }
      })
  }, .progress = TRUE)
  
  corr_mat <- data.frame(corr_mat) %>% t() %>% as.matrix()
  dimnames(corr_mat) <- NULL
  
  corr_mat <- ifelse(is.na(corr_mat), 0, corr_mat)
  corr_mat <- corr_mat + t(corr_mat) - diag(diag(corr_mat))
  
  dimnames(corr_mat) <- list(products, products)
  
  corr_mat
}


# replace half of the calibration products by alternative labels
sample_bool_approx_equal_n <- function(n) { 
  pool_samples = rep(c(T,F), each=ceiling(n()/2)); 
  sample(pool_samples, size = n, replace = F) 
}

