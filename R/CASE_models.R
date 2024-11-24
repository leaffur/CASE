#' @title CASE Model Training
#' 
#' @description Fit the priors for the cis-eQTL effect sizes.
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
#' @importFrom magrittr %>%
#' @importFrom stats pnorm qchisq
#' @export
CASE_train <- function(Z = NULL, R, N, hatB = NULL, hatS = NULL, 
                     V = NULL, h = NULL,
                     pi.init = NULL, U.orig = NULL,  
                     V.fix = TRUE, pi.fix = FALSE, 
                     n.iter = 45, MC.max = 125, tol = 1e-2, MC.seq = "exp", 
                     significant_thres = 0.1, noise_threshold = 0.1){
  ## hatB: G list with matrix M_j * C as the eQTL summary stats for each gene.
  ## R: Glist with M_j * M_j matrix: LD matrix for each region.
  ## N: C * C matrix, sample size and overlapped structure for each cell type.
  ##    Can be C vector which assumes full overlaps
  ## h: G * C matrix of cis-heritability.
  ## V: C * C covariance for the noise between cell types.
  
  
  ## U.orig: t1 + t2 lists / array with matrices C * C, initial cell type correlations for EM algorithm.
  ## pi.init: t1 + t2 + 1 vector, initial guess for prob of correlations.
  ##          The last is prob of the delta mass of zero.
  if (is.null(Z)){
    Z = hatB / hatS
  }
  
  hatBS = transform_Z(Z, N)
  hatB = hatBS$hatB
  hatS = hatBS$hatS
  
 
  if (MC.seq == "fix"){
    MC.sim = rep(MC.max, n.iter)
  }else if (MC.seq == "linear"){
    MC.sim = seq(25, MC.max, length.out = n.iter) %>% round
  }else if (MC.seq == "exp"){
    MC.sim = ((MC.max / 60)^((1:n.iter-1) / (n.iter - 1)) * 60) %>% round
  }

  C <- ncol(hatB)
  
  if (is.null(dim(N))){
    N = diag(N)
    for (i in 1:(C-1)){
      for (j in (i+1):C){
        N[j, i] = N[i, j] = min(N[i, i], N[j, j])
      }
    }
  }
  
  if (is.null(V)){
    V = diag(rep(1, C))
    for (i in 1:C){
      for (j in 1:C){
        V[i, j] = V[i, j] * N[i, j] / (N[i, i] * N[j, j])
      }
    }
  }
  
  # Initialization
  M <- nrow(R)
  init = Initialize_pi_U(hatB, hatS, C, M)
  
  if (is.null(pi.init)){
    pi.init = init$pi
    U.orig = init$U
  }
  pi = pi.init
  U = U.orig
  L = length(pi)
  
  
  ## Train with only marginally significant SNPs
  M0 = nrow(R)
  ME_p <- 2 - 2 * pnorm(abs(hatB / hatS), 0, 1)
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
    return(list(pi = 1, U = list(matrix(0, C, C)), V = V, n.iter = 0))
  }
  
  
  # if (!V.fix){
  #   iRB <- vector("list", G)
  #   for (j in 1:G){
  #     iRB[[j]] = t(hatB[[j]]) %*% solve(R[[j]]) %*% hatB[[j]]
  #   }
  # }
  
  J <- 0
  g <- list()

  kk = 1
  patterns.old = names(U)
  repeated_pattern = 0
  
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
                    R = R, pi = pi, h = h,
                    TT = gBc$TT, TT_det = gBc$TT_det, mu1 = gBc$mu1, Sigma1 = gBc$Sigma1)
      
      nsim  = nsim + 1
      gg[, nsim] = ifelse(BB[, , nsim] != 0, 1, 0) %>% apply(1, paste, collapse = "")
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
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
    
    for (k in 1:M){
      Sigma[[k]] = vector("list", L - 1)
      for (l in seq(L-1)){
        ind = which(gg[k, samp.ind] == patterns[l]) + samp.ind[1] - 1
        if (length(ind) > 0){
          Sigma[[k]][[l]] = tcrossprod(BB[k, , ind]) / length(ind)
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
      if (max(abs(pi.new - pi) / pi) <= tol){
        break
      }
    }
    
    # M-step
    # Update pi's
    if (!pi.fix){
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
      # cat("pi delete", ll, "\n")
      U = U[-ll]
      pi = pi[-ll]
    }
    
    if (length(pi) <= 1){
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
    
    pi = pi / sum(pi)
    L = length(U)
    
    # Update and cut U's
    ll = integer(0)
    dd = 0
    for (l in 1:(L-1)){
      ind = which(diag(U[[l]]) <= qchisq(0.95, 1) * diag(V))
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
      # cat("U delete", ll, "\n")
      U = U[-ll]
      pi = pi[-ll]
    }
    
    L = length(U)
    if (L <= 1){
      return(list(pi = pi, U = U, V = V, n.iter = kk, pi.in = pi.in, M1 = M1))
    }
    
    # cat("\n", names(pi), "\n")
    # print(pi)
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
#' @param hatB M * C matrix of standard errors of the estimated effects. Alternative summary data (together with hatB) to be provided instead of Z.
#' @return A \code{"CASE"} object with the following elements:
#' \item{pi:}{L-vector, the prior probabilities of sharing patterns.}
#' \item{U:}{L-list of C * C matrix, the prior covariances of sharing patterns.}
#' \item{V:}{C * C matrix, the sample-adjusted phenotypical variance.}
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @export
CASE_test <- function(hatB = NULL, Z = NULL, R, N, V, U, pi, MC.sim = 41, MC.sample = 58){
  # Here V is V adjusted for sample sizes
  #### Testing ####
  if (is.null(Z)){
    Z = hatB / hatS
  }
  
  hatBS = transform_Z(Z, N)
  hatB = hatBS$hatB
  hatS = hatBS$hatS
  
  pi = pi / sum(pi)
  C <- ncol(hatB)
  M <- nrow(R)
  
  if (length(pi) <= 1){
    post_mean = matrix(0, M, C)
    pvalue = matrix(1, M, C)
    return(list(pi = pi, U = U, V = V, pvalue = pvalue, post_mean = post_mean))
  }
  
  # m_split = U_split(U, pi, V)
  # U = m_split$U
  # pi = m_split$pi
  
  L = length(U)
  ## MC step
  gBc = gB_coef(U, V)
  pp = pm = array(0, dim = c(M, C, MC.sample))

  for (ll in 1:MC.sample){
    nsim = 1
    BB = array(0, dim = c(M, C, MC.sim))
    while (nsim < MC.sim){
      gB = gBupdate(B = BB[, , nsim],
                    hatB = hatB,
                    R = R, pi = pi,
                    TT = gBc$TT, TT_det = gBc$TT_det,
                    mu1 = gBc$mu1, Sigma1 = gBc$Sigma1)

      nsim  = nsim + 1
      BB[, , nsim] = gB
    }

    pp[, , ll] = apply(BB[, , -(1:ceiling(MC.sim * 0.3))], 1:2, function(x) mean(x == 0))
    pm[, , ll] = apply(BB[, , -(1:ceiling(MC.sim * 0.3))], 1:2, mean)

    if (ll == 20 & max(apply(pp[, , 1:ll], 1:2, sd)) <= 0.05){
      break
    }
  }
  pvalue = apply(pp[, , 1:ll], 1:2, mean)
  post_mean = apply(pm[, , 1:ll], 1:2, mean)

  return(list(pi = pi, U = U, V = V, pvalue = pvalue, post_mean = post_mean))
}


##########

#' CASE Obtain Credible Sets
#'
#' Obtain credible sets for any multi-trait fine-mapping results.
#' @param pvalues (M * C),The pvalues of SNPs.
#' @return Credible Sets
#' @importFrom magrittr %>%
#' @export
get_credible_sets <- function(pvalues, R, cor.min = 0.5, pip = 0.95, ruled_out = 1 - 1e-4){
  C = ncol(pvalues)
  css = vector("list", C)
  R1 = R
  
  for (ct in 1:C){
    p = pvalues[, ct]
    
    or = order(p)
    cs = list()
    coverage = numeric(0)
    purity = numeric(0)
    L = 0
    flag = rep(TRUE, nrow(R))
    
    for (kk in or){
      if (p[kk] > 0.95){
        break
      }
      
      if (p[kk] <= 1 - pip){
        L = L + 1
        cs[[L]] = kk
        flag[cs[[L]]] = FALSE
        purity[L] = 1
        coverage[L] = 1 - p[kk]
      } else{
        inds = which(abs(R[kk, ]) >= cor.min & flag)
        if (sum(1 - p[inds]) >= pip){
          or_inds = order(p[inds])
          or_inds = or_inds[p[inds[or_inds]] <= ruled_out]
          if (sum(1-p[inds[or_inds]]) > pip){
            best_local_sets = select_first_valid_set(1 - p[inds[or_inds]], R[inds[or_inds], inds[or_inds]],
                                                     pip, cor.min)
            if (!is.null(best_local_sets)){
              L = L + 1
              cs[[L]] = inds[or_inds[best_local_sets]]
              coverage[L] = min(1, sum(1 - p[cs[[L]]]))
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