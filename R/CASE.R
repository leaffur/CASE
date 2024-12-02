#' @title CASE Model
#' 
#' @description Perform Multi-trait Fine-mapping
#' @author Chen Lin, Hongyu Zhao
#' @details TBD
#' @references TBD
#' @param Z M * C matrix of z scores.
#' @param R M * M matrix of LD.
#' @param hatB M * C matrix of the estimated effects. Alternative summary data (together with hatS) to be provided instead of Z.
#' @param hatS M * C matrix of standard errors of the estimated effects. Alternative summary data (together with hatB) to be provided instead of Z.
#' @param N either C vector of the sample size, or C * C matrix of the sample size (diagonal) and ovelaps  (off-diagonal). If provided with a vector, CASE assumes that each pair of traits overlaps with their minimal sample size.
#' @param V (optional) C * C covariance (correlation) matrix for the noise between traits. If not provided, the default is an identity matrix.
#' @param cs logical, whether to get credible sets.
#' @param ... additional arguments.
#' @return A \code{"CASE"} object with the following elements:
#' \item{pi:}{L-vector, the prior probabilities of sharing patterns.}
#' \item{U:}{L-list of C * C matrix, the prior covariances of sharing patterns.}
#' \item{V:}{C * C matrix, the sample-adjusted phenotypical variance.}
#' @examples
#' ## A single-trait example.
#' set.seed(3)
#' N = 500
#' M = 1000
#' X = matrix(sample(0:2, size = N * M, replace = TRUE), N, M)
#' B = rep(0, M)
#' B[1:4] = runif(4, 1, 2)
#' e = rnorm(N)
#' Y = X %*% B + e
#' 
#' X = scale(X)
#' Y = scale(Y)
#' R = cor(X)
#' 
#' Z = hatB = hatS = rep(0, M)
#' hatB <- as.vector(t(X) %*% Y / (N - 1))
#' hatS <- sqrt((sum(Y^2)  - hatB^2 * (N - 1)) / N / (N - 2))
#' Z <- hatB / hatS
#' 
#' fit <- CASE(Z = Z, R = R, N = N)
#' # print(fit$sets)
#' 
#' 
#' ## A multi-trait example.
#' set.seed(3)
#' N = 500
#' M = 1000
#' C = 3
#' X = matrix(sample(0:2, size = N * M, replace = TRUE), N, M)
#' B = matrix(0, M, C)
#' B[1, ] = runif(C, 1, 2)
#' B[2, 1] = runif(1, 2, 3)
#' e = matrix(rnorm(N * C), N, C)
#' Y = X %*% B + e
#' 
#' X = scale(X)
#' Y = scale(Y)
#' R = cor(X)
#' 
#' Z = hatB = hatS = matrix(0, M, 3)
#' for (c in 1:3){
#'   hatB[, c] <- t(X) %*% Y[, c] / (N - 1)
#'   hatS[, c] <- sqrt((sum(Y[, c]^2)  - hatB[, c]^2 * (N - 1)) / N / (N - 2))
#'   Z[, c] <- hatB[, c] / hatS[, c]
#' }
#' 
#' fit <- CASE(Z = Z, R = R, N = rep(N, 3))
#' # print(fit$sets)
#' @export
CASE <- function(Z = NULL, R, hatB = NULL, hatS = NULL, N, V = NULL, cs = TRUE, ...){
    t1 = Sys.time()
    if (is.null(Z)){
        Z = hatB / hatS
    }
    Z = as.matrix(Z)
    hatBS = transform_Z(Z, N)
    hatB = hatBS$hatB
    hatS = hatBS$hatS
    
    # V = estimate_null_correlation_simple(mash_set_data(do.call(rbind, raw.data$hatB), do.call(rbind, raw.data$hatS)))
    
    m1 <- CASE_train(hatB = hatB, hatS = hatS, V = V, R = R, N = N, Z = NULL, list(...))
    
    res <- CASE_test(hatB = hatB, hatS = hatS, R = R, N = N, CASE_training = m1, Z = NULL, list(...))
    t2 = Sys.time()
    res$time = difftime(t2, t1, units = "secs")
    if (cs){
        res$sets <- get_credible_sets(res$pvalue, R = R)
    }
    
    return(res)
}
