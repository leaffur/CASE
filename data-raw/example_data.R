## code to prepare `example_data` dataset goes here
library(CorBin)

set.seed(1024)
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
B[idx2, ] = c(sqrt(0.25), sqrt(0.15), 0)

e = matrix(rnorm(N * C), N, C)
Y = X %*% B + e

example_data = list(Y = Y, X = X, B = B)

usethis::use_data(example_data, overwrite = TRUE)
