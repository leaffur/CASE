transform_Z <- function(Z, N){
  C = ncol(Z)
  hatB = hatS = Z
  if (length(dim(N)) == 2){
    N = diag(N)
  }
  for (c in 1:C){
    hatS[, c] <- 1 / sqrt(N[c] - 1 + Z[, c]^2)
    hatB[, c] <- Z[, c] * hatS[, c]
  }
  return(list(hatB = hatB, hatS = hatS))
}

replace_nonzero_with_mean <- function(row) {
  # Find non-zero and non-last elements
  non_zero_non_last_indices <- which(row != 0 & seq_along(row) != length(row))
  
  # Compute the mean of the non-zero and non-last elements
  if (length(non_zero_non_last_indices) > 0) {
    mean_value <- mean(row[non_zero_non_last_indices])
    row[non_zero_non_last_indices] <- mean_value
  }
  
  return(row)
}

select_first_valid_set <- function(p, R, threshold_sum = 0.95, min_corr = 0.5) {
  N <- length(p)
  
  # Helper function to check if adding a new element keeps the correlation valid
  is_valid_set <- function(selected, new_idx) {
    if (length(selected) == 0) {
      return(TRUE)
    }
    # Check if all correlations with the current selected set are >= min_corr
    return(all(abs(R[new_idx, selected]) >= min_corr))
  }
  
  # Recursive backtracking function to find the first valid set
  backtrack <- function(selected, current_sum, idx) {
    # If the current sum exceeds the threshold, return the selected set
    if (current_sum >= threshold_sum) {
      return(selected)
    }
    if (current_sum + sum(p[idx:N]) < threshold_sum){
      return(NULL)
    }
    if (idx < N) {
      if (length(selected) == 1){
        loops = 1:(N-idx+1)
      }else{
        loops = which(apply(abs(R[idx:N, selected]), 1, min) >= min_corr)
      }
      loops = loops + idx - 1
      
      if (length(loops) > 0){
        for (i in 1:length(loops)) {
          if (current_sum + sum(p[loops[i:length(loops)]]) < threshold_sum){
            return(NULL)
          }
          if (i < length(loops)){
            result <- backtrack(c(selected, loops[i]), current_sum + p[loops[i]], loops[i + 1])
          } else{
            if (current_sum + p[loops[i]] >= threshold_sum){
              return(c(selected, loops[i]))
            }
          }
          
          if (!is.null(result)) {
            return(result)
          }
        }
      }
    } else{
      if (all(abs(R[N, selected]) >= min_corr) & current_sum + p[N] >= threshold_sum){
        return(c(selected, N))
      }
    }
    return(NULL)
  }
  
  # Start backtracking with the first element already included
  return(backtrack(selected = 1, current_sum = p[1], idx = 2))
}

#' @importFrom MASS ginv
gB_coef <- function(U, V){
  C = nrow(U[[1]])
  L = length(U)
  
  TT = array(NA, dim = c(C, C, L))
  TT_det = rep(NA, L)
  mu1 = array(NA, dim = c(C, C, L))
  Sigma1 = array(NA, dim = c(C, C, L))
  
  for (l in 1:L){
    TT[, , l] = ginv(V + U[[l]])
    TT_det[l] = determinant(TT[, , l] %>% as.matrix, logarithm = TRUE)$modulus
    
    mu1[, , l] = U[[l]] %*% TT[, , l]
    Sigma1[, , l] = U[[l]] - U[[l]] %*% TT[, , l] %*% U[[l]]
    Sigma1[, , l] = (Sigma1[, , l] + t(Sigma1[, , l])) / 2
  }
  return(list(TT = TT, TT_det = TT_det, mu1 = mu1, Sigma1 = Sigma1))
}

#' @importFrom mvtnorm rmvnorm
#' @importFrom stats quantile rmultinom qnorm qchisq
gBupdate <- function(B, hatB, R, pi, h = NULL,
                     TT, TT_det, mu1, Sigma1, alpha = 0.05){
  
  if (is.null(alpha)){
    alpha = 0.05
  }
  
  B = B %>% as.matrix
  hatB = hatB %>% as.matrix
  
  C = dim(hatB)[2]
  M = nrow(B)
  L = length(pi)
  
  # update g by coord
  res = hatB - R %*% B + B # B + V
  
  prob = ((res %*% TT[, , L]) * res)
  inds1 = which(prob >= qchisq(1 - 5e-3 / 2, 1), arr.ind = TRUE)[, 1] %>% unique
  prob = rowSums(prob)
  inds2 = which(prob > max(quantile(prob, 0.9), qchisq(1 - 5e-3 / 2, C)))
  
  
  SNPS = union(inds1, inds2)
  prob = prob[SNPS] / 3
  
  or = integer(0)
  while (length(prob) > 0){
    prob.in = exp(prob - max(prob))
    index.in = which(prob.in >= 0.01)
    if (length(index.in) == 1){
      sel = index.in
    }else{
      Mc = sum(prob.in >= 0.1)
      sel = sample(index.in, Mc, replace = FALSE, prob = prob.in[index.in])
    }
    
    or = c(or, SNPS[sel])
    prob = prob[-sel]
    SNPS = SNPS[-sel]
  }
  
  for (i in or){
    res = hatB[i, ] - t(B) %*% R[, i] + B[i, ]
    
    pi1 = TT_det / 2 + log(pi)
    
    for (l in 1:L){
      pi1[l] = pi1[l] - (t(res) %*% TT[, , l] %*% res) / 2
    }
    
    pi1 = exp(pi1 - max(pi1))
    g = rmultinom(1, 1, pi1) %>% apply(2, function(x) which(x == 1))
    
    if (g == L){
      B[i, ] = 0
    } else{
      mu = c(mu1[, , g] %*% res)
      sigma0 = Sigma1[, , g] %>% as.matrix
      
      if (!is.null(h)){
        fake_thres = rbind(sqrt(h / 40), qnorm(1 - alpha / 2) / sqrt(diag(TT[, , L] %>% as.matrix))) %>% apply(2, max)
      } else{
        fake_thres = qnorm(1 - alpha / 2) / sqrt(diag(TT[, , L] %>% as.matrix))
      }
      fake_idx = which(abs(mu) < fake_thres)
      mu[fake_idx] = 0
      sigma0[, fake_idx] = 0
      sigma0[fake_idx, ] = 0
      
      B[i, ] = mvtnorm::rmvnorm(1, mean = mu, sigma = sigma0)
      if (any(abs(B[i, ]) < 1e-5)){
        B[i, ][abs(B[i, ]) < 1e-5] = 0
      }
    }
  }
  
  B[-or, ] = 0
  
  return(B)
}

canonical_patterns <- function(C) {
  single_bit_strings <- sapply(1:C, function(i) {
    binary_string <- rep(0, C)  # Create a vector of C zeroes
    binary_string[i] <- 1       # Set the ith position to 1
    paste(binary_string, collapse = "")  # Convert to string
  })
  
  all_one_string <- paste(rep(1, C), collapse = "")
  
  result <- c(single_bit_strings, all_one_string)
  return(result)
}

#' @importFrom stats p.adjust pnorm
adjust_sumstats <- function(hatB, hatS, method = "BH"){
  ME_p <- 2 - 2 * pnorm(abs(hatB / hatS), 0, 1)
  ME_p <- apply(ME_p, 2, function(x) p.adjust(x, method = method))
  return(ME_p)
}

Initialize_pi_U <- function(hatB, hatS, C, M, sig_threshold = 0.1){
  ME_p = adjust_sumstats(hatB, hatS)
  patterns <- ifelse(ME_p <= sig_threshold, 1, 0) %>% apply(1, paste, collapse = "")
  
  # Add canonical patterns and move zero pattern to the last.
  patterns = unique(c(patterns, canonical_patterns(C)))
  zero_pattern = paste(rep(0, C), collapse = "")
  if (zero_pattern %in% patterns) {
    patterns <- c(patterns[patterns != zero_pattern], zero_pattern)
  } else {
    patterns <- c(patterns, zero_pattern)
  }
  
  cat("Initialize with patterns: ", patterns, "\n")
  
  #################################################
  L = length(patterns)
  if (L > 6){
    pi = c(rep(5 / M / (L - 1), L - 1), (M - 5) / M)
  }else{
    pi = c(rep(1 / M, L - 1), (M - L + 1) / M)
  }
  
  U = list()
  for (l in patterns){
    U[[l]] = matrix(.1, C, C)
    diag(U[[l]]) = 1
    idx = which(as.numeric(strsplit(l, "")[[1]]) == 0)
    U[[l]][idx, ] = 0
    U[[l]][, idx] = 0
  }
  
  DD = apply(abs(hatB), 2, max)
  if (length(DD) > 1){
    U = lapply(U, function(x) (x * DD) %*% diag(DD))
  }else{
    U = lapply(U, function(x) x * DD^2)
  }
  
  return(list(pi = pi, U = U))
}
