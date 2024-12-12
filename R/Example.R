#' @importFrom mvsusieR mvsusie_rss create_mixture_prior 
#' @importFrom susieR susie_rss
#' @importFrom CorBin cBern
#' @importFrom stats cor lm rnorm
generate_example_data <- function(){
  set.seed(1130)
  N = 500
  M = 1000
  C = 3
  
  X = cbind(cBern(N, rep(0.3, M/4), rho = 0.95, type = "DCP"),
            cBern(N, rep(0.3, M/2), rho = 0.85, type = "DCP"),
            cBern(N, rep(0.3, M/4), rho = 0.98, type = "DCP"))
  
  B = matrix(0, M, C)
  idx1 = 10
  idx2 = 950
  B[idx1, ] = c(sqrt(0.3), sqrt(0.2), sqrt(0.1))
  B[idx2, ] = c(sqrt(0.2), sqrt(0.2), 0)
  e = matrix(rnorm(N * C), N, C)
  Y = X %*% B + e
  
  
  ### heritability
  # apply(X %*% B, 2, var) / apply(Y, 2, var)
  
  R = cor(X)
  
  Z = matrix(0, M, C)
  for (i in 1:M){
    m1 = summary(lm(Y ~ X[, i]))
    Z[i, ] = sapply(m1, function(x) x$coefficients[2, 3])
  }
  Z[c(idx1, idx2), ]
  B[c(idx1, idx2), ]
  
  fit <- CASE(Z = Z, R = R, N = rep(N, C))
  # fit$sets
  
  f1 = susie_rss(z = Z[, 1], R = R, n = N)
  # f1$sets
  
  prior <- create_mixture_prior(R = C)
  f2 = mvsusie_rss(Z = Z, R = R, N = N, prior_variance = prior)
  # f2$sets
  # f2$single_effect_lfsr
  
  # mvsusie_plot(f2)
}



