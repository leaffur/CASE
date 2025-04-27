#' @title CASE Model Training
#' 
#' @description Fit the priors for the cis-eQTL effect sizes.
#' @author Chen Lin, Hongyu Zhao
#' @details TBD
#' @references TBD
#' @param Z M * C matrix of z scores.
#' @param R M * M matrix of LD.
#' @param hatB M * C matrix of the estimated effects. Alternative summary data (together with hatS) to be provided instead of Z.
#' @param hatS M * C matrix of standard errors of the estimated effects. Alternative summary data (together with hatB) to be provided instead of Z.
#' @param N either C vector of the sample size, or C * C matrix of the sample size (diagonal) and ovelaps  (off-diagonal). If provided with a vector, CASE assumes that each pair of traits overlaps with their minimal sample size.
#' @param V (optional) C * C covariance (correlation) matrix for the noise between traits. If not provided, the default is an identity matrix.
#' @param ... additional arguments.
#' @return A \code{"CASE_training"} object with the following elements:
#' \item{pi:}{L-vector, the prior probabilities of sharing patterns.}
#' \item{U:}{L-list of C * C matrix, the prior covariances of sharing patterns.}
#' \item{V:}{C * C matrix, the sample-adjusted phenotypical variance.}
#' @importFrom magrittr %>%
#' @importFrom stats pnorm qchisq cov2cor
#' @export
CASE_train <- function(Z = NULL, R, hatB = NULL, hatS = NULL, N, V = NULL, ...){
  args = list(...)
  
  cat("Start Prior fitting.", "\n")
  
  if (is.null(Z)){
    Z = hatB / hatS
  }
  Z = as.matrix(Z)

  n.iter = ifelse("n.iter" %in% names(args), args$h, 45)
  MC.max = ifelse("MC.max" %in% names(args), args$MC.max, 125)
  MC.sim = ((MC.max / 60)^((1:n.iter-1) / (n.iter - 1)) * 60) %>% round
  C <- ncol(Z)
  
  if (is.vector(N)){
    if (C == 1){
      N = matrix(N)
    }else if (length(N) == C){
      N = diag(N)
      for (i in 1:(C-1)){
        for (j in (i+1):C){
          N[j, i] = N[i, j] = min(N[i, i], N[j, j])
        }
      }
    }else if (length(N) == 1){
      N = matrix(N, C, C)
    }
  }

  hatBS = transform_Z(Z, N)
  hatB = hatBS$hatB
  hatS = hatBS$hatS
  
  if (is.null(V)){
    V = diag(rep(1, C))
  } else{
    V = cov2cor(V)
  }
  for (i in 1:C){
    for (j in 1:C){
      V[i, j] = V[i, j] * N[i, j] / (N[i, i] * N[j, j])
    }
  }

  # Initialization
  if ("pi.init" %in% names(args)){
    pi = args$pi.init
    U = args$U.orig
  } else{
    init = Initialize_pi_U(hatB, hatS, C, M = nrow(R))
    pi = init$pi
    U = init$U
  }
  L = length(pi)
  
  
  ## Train with only marginally significant SNPs
  M0 = nrow(R)
  ME_p <- 2 - 2 * pnorm(abs(hatB / hatS), 0, 1)
  significant_thres = ifelse("significant_thres" %in% names(args), args$significant_thres, 1e-1)
  idx <- which(ME_p <= significant_thres, arr.ind = TRUE)[, 1] %>% unique
  if (length(idx) == 1){
    idx = c(idx - 1, idx, idx + 1)
    idx = idx[idx >= 1 & idx <= nrow(hatB)]
  }
  hatB = hatB[idx, ]
  R = R[idx, idx]
  
  M1 = length(idx)
  pi.in = M1 / M0
  M <- nrow(R)
  if (M1 == 0){
    cat("No marginally significant variants in the inputs.", "\n")
    return(list(pi = 1, U = list(matrix(0, C, C)), V = V, n.iter = 0))
  }
  
  J <- 0
  g <- list()
  patterns.old = names(U)
  repeated_pattern = 0
  alpha = ifelse("alpha" %in% names(args), args$alpha, 0.05)

  # MCEM steps
  for (kk in 1:n.iter){
    # cat("\n iter:", kk, " ")

    # E-step
    ## MC step
    # init
    nsim = 1
    gg = matrix(paste(rep(0, C), collapse = ""), M, MC.sim[kk])
    BB = array(0, dim = c(M, C, MC.sim[kk]))
    
    # MC iter
    gBc = gB_coef(U, V)
    while (nsim < MC.sim[kk]){
      BB[, , nsim + 1] = gBupdate(B = BB[, , nsim], hatB = hatB,
                    R = R, pi = pi, h = args$h,
                    TT = gBc$TT, TT_det = gBc$TT_det, mu1 = gBc$mu1, Sigma1 = gBc$Sigma1, alpha = alpha)
      
      nsim  = nsim + 1
      gg[, nsim] = ifelse(as.matrix(BB[, , nsim]) != 0, 1, 0) %>% apply(1, paste, collapse = "")
    }
    
    samp.ind = ceiling(MC.sim[kk] / 3):MC.sim[kk]
    patterns = unique(c(gg[, samp.ind]))
    zero_pattern = paste(rep(0, C), collapse = "")
    if (zero_pattern %in% patterns) {
      patterns <- c(patterns[patterns != zero_pattern], zero_pattern)
    } else {
      patterns <- c(patterns, zero_pattern)
    }
    g = apply(gg[, samp.ind], 1, function(x) table(factor(x, levels = patterns))) %>% t 
    
    Sigma = vector("list", M)
    L = length(patterns)
    
    if (L <= 1){
      cat("Estimates no eQTL effects in the CASE prior fitting step.", "\n")
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
    
    for (k in 1:M){
      Sigma[[k]] = vector("list", L - 1)
      for (l in seq(L-1)){
        ind = which(gg[k, samp.ind] == patterns[l]) + samp.ind[1] - 1
        if (length(ind) > 0){
          if (C == 1){
            Sigma[[k]][[l]] = matrix(mean((BB[k, , ind])^2))
          }else{
            Sigma[[k]][[l]] = tcrossprod(BB[k, , ind]) / length(ind) 
          }
         
        } else{
          Sigma[[k]][[l]] = matrix(0, C, C)
        }
      }
    }
    
    # End of E-step
    
    # g_jkt
    if (kk <= 2){
      g = t(apply(g, 1, replace_nonzero_with_mean))
    }
    g <- g / sum(g)
    
    # terminations:
    pi.new  = colSums(g) * pi.in
    pi.new[L] = 1 - sum(pi.new[-L])
    
    if (length(pi) == L){
      tol = ifelse("tol" %in% names(args), args$tol, 1e-2)
      if (max(abs(pi.new - pi) / pi) <= tol){
        break
      }
    }
    
    # M-step
    # Update pi's
    if ("pi.fix" %in% names(args)){
      if (!args$pi.fix){
        pi = pi.new
      }
    }else{
      pi = pi.new
    }
    
    # Update U's
    U = vector("list", L)
    for (l in 1:(L-1)){
      U[[l]] <- matrix(0, C, C)
      for (k in 1:M){
        U[[l]] = U[[l]] + g[k, l] * Sigma[[k]][[l]]
      }
      U[[l]] <- U[[l]] / pi[l] * pi.in
      U[[l]][abs(U[[l]]) < 1e-5 / M] = 0
    }
    U[[L]] <- matrix(0, C, C)
    names(U) = patterns
    
    #########################################
    # Cut pi's
    if (kk <= 5){
      ll = which(pi < 1e-8)
    } else{
      ll = which(pi < 1e-2 / M * pi.in)
    }
    if (length(ll) > 0){
      U = U[-ll]
      pi = pi[-ll]
    }
    
    if (length(pi) <= 1){
      cat("Estimated no eQTL effects in the CASE prior fitting step.", "\n")
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
    
    pi = pi / sum(pi)
    L = length(U)
    
    # Update and cut U's
    ll = integer(0)
    dd = 0
    for (l in 1:(L-1)){
      ind = which(diag(U[[l]]) <= qchisq(1 - alpha, 1) * diag(V))
      U[[l]][ind, ] = 0
      U[[l]][, ind] = 0
      lp = 0
      flag = TRUE
      while (flag & (lp < L)){
        lp = lp + 1
        if (mean( (U[[l]] != 0) == (U[[lp]] != 0) ) == 1 & (l != lp) & !(lp %in% ll)){
          dd = dd + 1
          ll[dd] = l
          U[[lp]] = (pi[l] * U[[l]] + pi[lp] * U[[lp]]) / (pi[l] + pi[lp])
          pi[lp] = pi[l] + pi[lp]
          flag = FALSE
        }
      }
    }
    
    if (length(ll) > 0){
      U = U[-ll]
      pi = pi[-ll]
    }
    
    L = length(U)
    if (L <= 1){
      cat("Estimated no eQTL effects in the CASE prior fitting step.", "\n")
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
    
    if (setequal(patterns.old, names(pi))){
      repeated_pattern = repeated_pattern + 1
    } else{
      repeated_pattern = 0
      patterns.old = names(pi)
    }
    
    if (repeated_pattern >= 10){
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
  }
  
  return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
}


#' @title CASE Model Testing
#' 
#' @description Obtain the posterior probabilities and estimates of the cis-eQTL effect sizes.
#' @author Chen Lin, Hongyu Zhao
#' @details TBD
#' @references TBD
#' @param Z M * C matrix of z scores.
#' @param R M * M matrix of LD.
#' @param hatB M * C matrix of the estimated effects. Alternative summary data (together with hatS) to be provided instead of Z.
#' @param hatS M * C matrix of standard errors of the estimated effects. Alternative summary data (together with hatB) to be provided instead of Z.
#' @param N either 1 or C vector of the sample size, or C * C matrix of the sample size (diagonal) and overlaps  (off-diagonal). If provided with a vector, CASE assumes that each pair of traits overlaps with their minimal sample size.
#' @param CASE_training A \code{"CASE_training"} object.
#' @param ... additional arguments.
#' @return A \code{"CASE"} object with the following elements:
#' \item{pi:}{L-vector, the prior probabilities of sharing patterns.}
#' \item{U:}{L-list of C * C matrix, the prior covariances of sharing patterns.}
#' \item{V:}{C * C matrix, the sample-adjusted phenotypical variance.}
#' \item{pip:}{M * C matrix, posterior probability of having eQTL effects per SNP per cell type.}
#' \item{post_mean:}{M * C matrix, average posterior estimates of eQTL effects per SNP per cell type.}
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @export
CASE_test <- function(Z = NULL, R, hatB = NULL, hatS = NULL, N, CASE_training, ...){
  # Here V is V adjusted for sample sizes
  #### Testing ####
  args = list(...)
  cat("Start Posterior Analysis.", "\n")

  U = CASE_training$U
  V = CASE_training$V
  pi = CASE_training$pi
  
  if (is.null(Z)){
    Z = hatB / hatS
  }
  Z = as.matrix(Z)
  
  hatBS = transform_Z(Z, N)
  hatB = hatBS$hatB
  hatS = hatBS$hatS

  C <- ncol(hatB)
  M <- nrow(R)
  
  if (length(pi) <= 1){
    post_mean = matrix(0, M, C)
    pvalue = matrix(1, M, C)
    return(list(pi = pi, U = U, V = V, pvalue = pvalue, post_mean = post_mean))
  }
  
  L = length(U)
  ## MC step
  gBc = gB_coef(U, V)
  
  MC.sim = ifelse("MC.sim" %in% names(args), args$MC.sim, 41)
  MC.sample = ifelse("MC.sample" %in% names(args), args$MC.sample, 58)
  pp = pm = array(0, dim = c(M, C, MC.sample))

  for (ll in 1:MC.sample){
    nsim = 1
    BB = array(0, dim = c(M, C, MC.sim))
    while (nsim < MC.sim){
      gB = gBupdate(B = BB[, , nsim],
                    hatB = hatB,
                    R = R, pi = pi,
                    TT = gBc$TT, TT_det = gBc$TT_det,
                    mu1 = gBc$mu1, Sigma1 = gBc$Sigma1, alpha = args$alpha)

      nsim  = nsim + 1
      BB[, , nsim] = gB
    }

    pp[, , ll] = apply(BB[, , -(1:ceiling(MC.sim * 0.3))], 1:((C > 1) + 1), function(x) mean(x != 0))
    pm[, , ll] = apply(BB[, , -(1:ceiling(MC.sim * 0.3))], 1:((C > 1) + 1), mean)

    if (ll == 20){
      if (max(apply(pp[, , 1:ll], 1:((C > 1) + 1), sd)) <= 0.05){
        break
      }
    }
  }
  pip = apply(pp[, , 1:ll], 1:((C > 1) + 1), mean)
  post_mean = apply(pm[, , 1:ll], 1:((C > 1) + 1), mean)

  return(list(pi = pi, U = U, V = V, pip = pip, post_mean = post_mean))
}


##########

#' CASE Obtain Credible Sets
#'
#' Obtain credible sets for any multi-trait fine-mapping results.
#' @param pips (M * C),The pips of SNPs.
#' @param R M * M matrix of LD.
#' @param cor.min minimum correlation in the credible sets
#' @param coverage_thres threshold for the sum of PIPs.
#' @param ruled_out excluding SNPs with PIPs less than the threshold.
#' @return a length C list of credible sets.
#' @importFrom magrittr %>%
#' @export
get_credible_sets <- function(pips, R, cor.min = 0.5, coverage_thres = 0.95, ruled_out = 1e-4){
  cat("Start getting credible sets.", "\n")
  
  pips = as.matrix(pips)
  C = ncol(pips)
  css = vector("list", C)
  R1 = R
  
  for (ct in 1:C){
    p = pips[, ct]
    
    or = order(p, decreasing = TRUE)
    cs = list()
    coverage = numeric(0)
    purity = numeric(0)
    L = 0
    flag = rep(TRUE, nrow(R))
    
    for (kk in or){
      if (p[kk] < 0.05){
        break
      }
      
      if (p[kk] >= coverage_thres){
        L = L + 1
        cs[[L]] = kk
        flag[cs[[L]]] = FALSE
        purity[L] = 1
        coverage[L] = p[kk]
      } else{
        inds = which(abs(R[kk, ]) >= cor.min & flag)
        if (sum(p[inds]) >= coverage_thres){
          or_inds = order(p[inds])
          or_inds = or_inds[p[inds[or_inds]] >= ruled_out]
          if (sum(p[inds[or_inds]]) > coverage_thres){
            best_local_sets = select_first_valid_set(p[inds[or_inds]], R[inds[or_inds], inds[or_inds]],
                                                     coverage_thres, cor.min)
            if (!is.null(best_local_sets)){
              L = L + 1
              cs[[L]] = inds[or_inds[best_local_sets]]
              coverage[L] = min(1, sum(p[cs[[L]]]))
              mat = abs(R[cs[[L]], cs[[L]]])
              purity[L] = min(mat[upper.tri(mat)])
              flag[cs[[L]]] = FALSE
            }
          }
        }
      }
      flag[kk] = FALSE
    }
    
    css[[ct]] = list(cs = cs, purity_min_cor = purity, coverage = coverage)
  }
  return(css)
}
