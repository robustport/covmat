#' Calculate the average correlation
#'
#' @details
#' Calculates the average correlation based on the sample covariance matrix.
#' Assumes that the diagonal of the sample covariance matrix contains variances.
#'
#' @param sample_cov_matrix Sample covariance matrix
#' @param robust use robust estimator of average correlation
#' @return Average correlation across variables
#' @noRd
calculate_avg_correlation <- function(sample_cov_matrix, robust) {
  cor_matrix <- cov2cor(sample_cov_matrix)
  p <- ncol(cor_matrix)
  rho <- ifelse(robust, median(cor_matrix), mean(cor_matrix))
  avg_correlation <- (p*rho - 1) / (p - 1)
  return(avg_correlation)
}

#' Calculate the constant correlation shrinkage target matrix
#'
#' @details
#' Calculates the target matrix for shrinkage based on the sample covariance
#' matrix, number of variables, and average correlation.
#'
#' @param sample_cov_matrix Sample covariance matrix
#' @param robust use robust estimator of average correlation
#' @return Target matrix for shrinkage
#' @noRd
calculate_target <- function(sample_cov_matrix, robust) {
  avg_correlation <- calculate_avg_correlation(sample_cov_matrix, robust)
  vol <- sqrt(diag(sample_cov_matrix))
  p <- ncol(sample_cov_matrix)
  const_cor_matrix <- matrix(avg_correlation, nrow = p, ncol = p)
  diag(const_cor_matrix) <- 1
  diag_vol <- diag(vol)
  target <- diag_vol %*% const_cor_matrix %*% diag_vol
  return(target)
}

#' Calculate the pi_hat parameter
#'
#' @details
#' Calculates the pi_hat parameter based on the data matrix X, effective sample
#' size, number of variables, and the sample covariance matrix.
#'
#' @param X Data matrix
#' @param n Effective sample size
#' @param sample_cov_matrix Sample covariance matrix
#' @return pi_hat parameter
#' @noRd
calculate_pi_hat <- function(X, sample_cov_matrix, n, diag=FALSE) {
  X_sq <- X^2
  if (diag) {
   pi_hat <- sum(diag((t(X_sq) %*% X_sq) / n - sample_cov_matrix^2))
  } else {
   pi_hat <- sum((t(X_sq) %*% X_sq) / n - sample_cov_matrix^2)
  }

  return(pi_hat)
}

#' Calculate the gamma_hat parameter
#'
#' @details
#' Calculates the gamma_hat parameter based on the sample and target covariance
#' matrices. Uses Frobenius norm.
#'
#' @param sample Sample covariance matrix
#' @param target Target covariance matrix for shrinkage
#' @return gamma_hat parameter
#' @noRd
calculate_gamma_hat <- function(sample, target) {
  gamma_hat <- norm(sample - target, type = "F")^2
  return(gamma_hat)
}

#' Calculate the rho_hat parameter
#'
#' @details
#' Calculates the rho_hat parameter based on the data matrix X, effective sample
#' size, number of variables, sample covariance matrix, and average correlation.
#'
#' @param X Data matrix
#' @param n Effective sample size
#' @param p Number of variables
#' @param sample_cov_matrix Sample covariance matrix
#' @param n degrees of freedom
#' @param robust use robust estimator of average correlation
#' @return rho_hat parameter
#' @noRd
calculate_rho_hat <- function(X, sample_cov_matrix, n, robust) {
  p <- ncol(X)
  sample_var <- diag(sample_cov_matrix)
  sample_vol <- sqrt(sample_var)
  avg_correlation <- calculate_avg_correlation(sample_cov_matrix, robust)
  Z <- (t(X^3) %*% X) / n - outer(rep(1, p), sample_var) * sample_cov_matrix
  coeffs <- outer(sample_vol, 1 / sample_vol)
  temp_res <- coeffs * Z
  rho_off <- avg_correlation * (sum(temp_res) - sum(diag(temp_res)))
  rho_diag <- calculate_pi_hat(X, sample_cov_matrix, n, diag = TRUE)
  rho_hat <- rho_diag + rho_off
  return(rho_hat)
}

#' Calculate the shrinkage intensity
#'
#' @details
#' Calculates the shrinkage intensity based on pi_hat, rho_hat, gamma_hat, and
#' effective sample size.
#'
#' @param pi_hat pi_hat parameter
#' @param rho_hat rho_hat parameter
#' @param gamma_hat gamma_hat parameter
#' @param n Effective sample size
#' @return Shrinkage intensity
#' @param robust use robust estimator of average correlation
#' @noRd
calculate_shrinkage_intensity <- function(X, sample_cov_matrix, target, n, robust) {
  pi_hat <- calculate_pi_hat(X, sample_cov_matrix, n)
  gamma_hat <- calculate_gamma_hat(sample_cov_matrix, target)
  rho_hat <- calculate_rho_hat(X, sample_cov_matrix, n, robust)
  kappa_hat <- (pi_hat - rho_hat) / gamma_hat
  shrinkage <- max(0, min(1, kappa_hat / n))
  return(shrinkage)
}

#' Calculate the shrinkage estimator of the covariance matrix
#'
#' @details
#' Calculates the shrinkage estimator of the covariance matrix based on the
#' shrinkage intensity, target matrix, and sample covariance matrix.
#'
#' @param shrinkage Shrinkage intensity
#' @param target Target covariance matrix for shrinkage
#' @param sample Sample covariance matrix
#' @return Shrinkage estimator of the covariance matrix
#' @noRd
calculate_shrinkage_estimator <- function(shrinkage, target, sample) {
  sigma_hat <- shrinkage * target + (1 - shrinkage) * sample
  return(sigma_hat)
}

#' Calculate the shrinkage estimator of covariance using the procedure
#' described in (Ledoit, Olivier, and Michael Wolf. "Honey, I shrunk
#' the sample covariance matrix.".
#'
#' @details
#' This function calculates the Ledoit-Wolf shrinkage covariance estimator
#' to constant correlation unequal variance target matrix.
#'
#' @param X Data matrix.
#' @param demean remove mean from the features
#' @param robust use robust estimator of average correlation
#' @return A list containing the shrinkage estimator of the covariance matrix
#'         and the shrinkage intensity.
#' @export
#' @author Rohit Arora
lw_linear_shrinkage <- function(X, demean = TRUE, robust = FALSE) {  
  sample_mean <- colMeans(X)
  mean_outer_product <- tcrossprod(sample_mean)
  
  if (demean) {
    X <- scale(X, center = TRUE, scale = FALSE)
  }
  
  df <- ifelse(demean, nrow(X) - 1, nrow(X))
  sample_cov_matrix <- (t(X) %*% X) / df
  
  target <- calculate_target(sample_cov_matrix, robust)
  shrinkage_intensity <- calculate_shrinkage_intensity(X, 
                                        sample_cov_matrix, target, df, robust)
  sigma_hat <- calculate_shrinkage_estimator(shrinkage_intensity, 
                                             target, sample_cov_matrix)
  
  if (demean) {
    sigma_hat <- sigma_hat + mean_outer_product
  }
  
  return(list(cov = sigma_hat, shrinkage_intensity = shrinkage_intensity))
}