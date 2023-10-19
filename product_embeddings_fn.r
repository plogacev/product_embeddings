
get_word2vec <- function(x, dim, window, iter)
{
  fname <- sprintf("models/w2v_%d_%d_%d.model", dim, window, iter)
  if (file.exists(fname)) {
    model <- word2vec::read.word2vec(file = fname)
    
  } else {
    set.seed(123456789)
    
    tic()
    model <- word2vec(x = x, type = "cbow", dim = dim, window = window, iter = iter, threads = 8L)
    word2vec::write.word2vec(model, file = fname)
    toc()
  }
  model
}

evaluate_model_by_similarity_proximity <- function(calibration_products, dim, window, iter = 10)
{
  # using: log_list_prices, transaction_baskets_df
  
  model <- get_word2vec(x = transaction_baskets_df$basket, dim = dim, window = window, iter = iter)
  embeddings <- as.matrix(model)
  
  vocabulary <- rownames(embeddings)
  selected_calibration_products <- calibration_products[calibration_products %in% vocabulary]
  
  word2vec_similarities <- word2vec::word2vec_similarity(
    embeddings[selected_calibration_products,],
    embeddings[selected_calibration_products,], 
    type = "cosine")/2 + 0.5
  log_prior_similarities <- log_prior_similarity_mat[selected_calibration_products,
                                                     selected_calibration_products]
  prior_similarities <- exp(log_prior_similarities)

  proximity_metric <- function(prior_similarities, word2vec_similarities)
  {
    word2vec_similarities_lower_tri <- word2vec_similarities %>% .[lower.tri(.)]
    prior_similarities_lower_tri <- prior_similarities %>% .[lower.tri(.)]
    
    similarity    <- -1 * prior_similarities_lower_tri * log(word2vec_similarities_lower_tri)
    dissimilarity <- -1 * (1-prior_similarities_lower_tri) * log(1-word2vec_similarities_lower_tri)
    
    -sum(log(similarity + dissimilarity))
  }
  
  proximity_metric(prior_similarities, word2vec_similarities)
}


evaluate_model_by_dept_match <- function(selected_products, products_dept_match_mat, dim, window, iter = 10)
{
  # using: log_list_prices, transaction_baskets_df
  
  model <- get_word2vec(x = transaction_baskets_df$basket, dim = dim, window = window, iter = iter)
  embeddings <- as.matrix(model)
  
  vocabulary <- rownames(embeddings)
  selected_products <- selected_products[selected_products %in% vocabulary] %>% as.character()
  
  word2vec_similarities <- word2vec::word2vec_similarity(embeddings[selected_products,], embeddings[selected_products,], type = "cosine")/2 + 0.5
  products_dept_match_mat <- products_dept_match_mat[selected_products, selected_products]

  logloss_dept_match <- function(products_dept_match_mat, word2vec_similarities)
  {
    word2vec_similarities_lower_tri <- word2vec_similarities %>% .[lower.tri(.)]
    products_dept_match_lower_tri <- products_dept_match_mat %>% .[lower.tri(.)]
    
    similarity    <- -1 * products_dept_match_lower_tri * log(word2vec_similarities_lower_tri)
    dissimilarity <- -1 * (1-products_dept_match_lower_tri) * log(1-word2vec_similarities_lower_tri)
    
    -sum(log(similarity + dissimilarity))
  }
  
  logloss_dept_match(products_dept_match_mat, word2vec_similarities)
}

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

