
library(Rcpp)

source_cpp <- '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::export]]
double calculate_distances(Rcpp::NumericVector x, Rcpp::NumericVector y)
{
  // Calculate dot product
  double similarities = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

  // Calculate Euclidean norms
  double x_scale = std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));
  double y_scale = std::sqrt(std::inner_product(y.begin(), y.end(), y.begin(), 0.0));

  // Normalize and return
  return similarities / (x_scale * y_scale);
}

// [[Rcpp::export]]
Rcpp::NumericVector extract_distances(Rcpp::List embeddings, int idx_product_1, int idx_product_2)
{
  int n_models = embeddings.size();
  Rcpp::NumericVector result(n_models);

  for (int i_model = 0; i_model < n_models; ++i_model) {
    Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(embeddings[i_model]);
    Rcpp::NumericVector x = mat(idx_product_1-1, _);
    Rcpp::NumericVector y = mat(idx_product_2-1, _);

    result[i_model] = calculate_distances(x, y);
  }

  return result;
}

// [[Rcpp::export]]
NumericVector summary_similarities_ij(Rcpp::List embeddings, int idx_product_1, int idx_product_2)
{
  NumericVector distances = extract_distances(embeddings, idx_product_1, idx_product_2);
  arma::vec similarities = distances/2 + 0.5;

  arma::vec target_quantiles = {0.025, 0.975};
  arma::vec quantiles = arma::quantile(similarities, target_quantiles );
  
  return NumericVector::create(
    mean(similarities),
    stddev(similarities),
    quantiles[0],
    quantiles[1]
  );
}
'

# Save the C++ code to a file and load using sourceCpp
cppFile <- tempfile(fileext = ".cpp")
writeLines(source_cpp, cppFile)
sourceCpp(cppFile)


summary_similarities_i <- function(embeddings, idx_products, idx_target_product)
{
  x <- purrr::map_dfc(idx_products, \(j) {
    summary_similarities_ij(embeddings, idx_product_1 = idx_target_product, idx_product_2 = j)
  })
  similarities_df <- matrix( unlist(x), ncol = 4, byrow = TRUE) %>% as.data.frame()
  colnames(similarities_df) <- c("est", "se", "lower", "upper")
  similarities_df$product_2 <- names(idx_products)
  similarities_df$product_1 <- names(idx_target_product)
  similarities_df %>% select(product_1, product_2, est, se, lower, upper)
}

compute_similarities <- function(embeddings, idx_products, idx_target_products)
{
  similarities_lst <- 
    purrr::map(seq_along(idx_target_products), \(i) {
      similarities_df <- summary_similarities_i(embeddings, idx_products, idx_target_products[i])
      similarities_df %<>% arrange(desc(est))
      similarities_df
    }, .progress = TRUE)
  names(similarities_lst) <- names(idx_target_products)
  similarities_lst
}


embeddings_ensemble_word2vec <- function(baskets_lst, dim, window, iter, n_models) {
  set.seed(12345)
  fname <- sprintf("models/w2v_ensemble_%d_%d_%d_%d.rds", dim, window, iter, n_models)
  if (file.exists(fname)) {
    embeddings <- readRDS(file = fname)
    
  } else {
    models <- purrr::map(1:n_models, \(i) {
      baskets <- baskets_lst %>% sapply( function(x) paste(sample(x), collapse = " "))
      word2vec(x = sample(baskets), type = "cbow", dim = dim, window = window, iter = iter, threads = 8L)
    }, .progress = TRUE)
    embeddings <- lapply(models, as.matrix)
    saveRDS(embeddings, fname)
  }
  embeddings
}



get_word2vec <- function(x, dim, window, iter, n_resampled_baskets = 1)
{
  set.seed(12345678)
  
  cat(".")
  fname <- sprintf("models/w2v_%d_%d_%d_%d.model", dim, window, iter, n_resampled_baskets)
  if (file.exists(fname)) {
    model <- word2vec::read.word2vec(file = fname)
    
  } else {
    
    tic()
    if (n_resampled_baskets > 1) {
      x_original <- x
      for (i in 2:n_resampled_baskets) {
        x_resampled <- sample(x_original, replace = T)
        x %<>% c( x_resampled )
      }
    }
    
    baskets <- x %>% sapply( function(x) paste(sample(x), collapse = " "))
    model <- word2vec(x = baskets, type = "cbow", dim = dim, window = window, iter = iter, threads = 8L)
    word2vec::write.word2vec(model, file = fname)
    toc()
  }
  
  model
}

evaluate_model_loglik <- function(baskets_lst, calibration_products, dim, window, iter, n_resampled_baskets = 1)
{
  model <- get_word2vec(x = baskets_lst, dim = dim, window = window, iter = iter, n_resampled_baskets = n_resampled_baskets)
  if (is.null(model)) {
    return(NA)
  }
  
  embeddings <- as.matrix(model)
  vocabulary <- rownames(embeddings) %>% intersect( rownames(match_mat) )
  selected_products <- calibration_products[calibration_products %in% vocabulary]
  remaining_products <- vocabulary %>% setdiff(selected_products)
  
  distances <- word2vec::word2vec_similarity( embeddings[selected_products,], embeddings, type = "cosine")/2 + 0.5
  similarities <- distances/2 + 0.5
  
  similarities %<>% t() %>% .[c(selected_products, remaining_products), selected_products]
  match_mat %<>% .[c(selected_products, remaining_products), selected_products]
  mismatch_mat %<>% .[c(selected_products, remaining_products), selected_products]

  diag(similarities[selected_products, selected_products]) <- NA
  similarities[selected_products, selected_products][upper.tri(matrix(nrow = 200, ncol = 200))] <- NA
  
  match_mat_lower_tri <- match_mat %>% .[lower.tri(.)]
  mismatch_mat_lower_tri <- mismatch_mat %>% .[lower.tri(.)]
  similarities_lower_tri <- similarities %>% .[lower.tri(.)]

  loglik <- match_mat_lower_tri * log(similarities_lower_tri) + mismatch_mat_lower_tri * log(1-similarities_lower_tri)
  
  sum(loglik)
}

# replace half of the calibration products by alternative labels
sample_bool_approx_equal_n <- function(n) { 
  pool_samples = rep(c(T,F), each=ceiling(n()/2)); 
  sample(pool_samples, size = n, replace = F) 
}

