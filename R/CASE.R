#' @title CASE Model
#' 
#' @description Perform Multi-trait Fine-mapping
#' @author Chen Lin, Hongyu Zhao
#' @details TBD
#' @references TBD
#' @param Z M * C matrix of z scores.
#' @param R M * M matrix of LD.
#' @param hatB M * C matrix of the estimated effects. Alternative summary data (together with hatS) to be provided instead of Z.
#' @param hatB M * C matrix of standard errors of the estimated effects. Alternative summary data (together with hatB) to be provided instead of Z.
#' @return A \code{"CASE"} object with the following elements:
#' \item{pi:}{L-vector, the prior probabilities of sharing patterns.}
#' \item{U:}{L-list of C * C matrix, the prior covariances of sharing patterns.}
#' \item{V:}{C * C matrix, the sample-adjusted phenotypical variance.}
#' @examples 
#' TBD
#' @export
CASE <- function(Z = NULL, R, N, hatB = NULL, hatS = NULL, 
                     V = NULL, h = NULL,
                     cs = TRUE,
                     pi.init = NULL, U.orig = NULL,  
                     V.fix = TRUE, pi.fix = FALSE, 
                     n.iter = 45, MC.max = 125, tol = 1e-2, MC.seq = "exp", 
                     MC.sim = 41, MC.sample = 58,
                     significant_thres = 0.1, noise_threshold = 0.1){
    t1 = Sys.time()
    if (is.null(Z)){
        Z = hatB / hatS
    }
      
    hatBS = transform_Z(Z, N)
    hatB = hatBS$hatB
    hatS = hatBS$hatS
    
    m1 <- CASE_train(hatB = hatB, hatS = hatS, R = R, N = N)
    
    a <- CASE_test(hatB = hatB, R = R, N = N, 
                      V = m1$V, U = m1$U, pi = m1$pi)
    t2 = Sys.time()
    a$time = difftime(t2, t1, units = "secs")
    if (cs){
        a$sets <- get_credible_sets(a$pvalue,  R = R)
    }
}
