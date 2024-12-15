#' Calculation of constant required for consistency
#' 
#' @details
#' The biweight rho function needs a constant whose value is determined based on
#' the consistency equation. Solving the integral and assuming multivariate
#' normality we can get an analytical expression for the constant.
#' 
#' @param p number of variables
#' @param c  square root of the chi-square cutoff
#' @noRd
.gamma_cp <- function(p, c) {
  
  IncGamma <- function(a,z) {
    pgamma(z, a, lower.tail = FALSE)*gamma(a)
  }
  
  p*(c^6*gamma(p/2))/(gamma(p/2 + 1)*(6*c^4 -6*c^2*(p+2) + 2*(p+2)*(p+4)) - 
                        6*c^4*IncGamma(1 + p/2, c^2/2) + 
                        12*c^2*IncGamma(2 + p/2, c^2/2) - 
                        8*IncGamma(3 + p/2, c^2/2) + c^6*IncGamma(p/2, c^2/2))
}

#' Givens rotation matrix
#' 
#' @details
#' Givens rotation matrix is an orthogonal matrix with certain structure.
#' 
#' @param i row number
#' @param j column number
#' @param theta angle in the range [-pi/2, pi/2]
#' @param d number of columns
#' @noRd
.givens.rotation <- function(i, j, theta, d) {
  
  G <- diag(rep(1, d))
  G[i,i] <-  G[j,j] <- cos(theta)
  G[j,i] <- sin(theta); G[i,j] <- -G[j,i]
  
  G
}

#' Product of Givens rotoation matrices
#' 
#' @details
#' An orthogonal matrix of size d can be parameterized by d(d-1)/2 angles
#' 
#' @param d number of columns
#' @param theta angle in the range [-pi/2, pi/2]
#' @noRd
.orthogonal.matrix <- function(d, angles) {
  
  O <- diag(rep(1, d)); k <- 1
  
  for(i in 2:d) 
    for(j in 1:(i-1)) {
      O <- O %*% .givens.rotation(i, j, theta = angles[k], d)
      k <- k + 1
    }
  O
}


#' Smoothing matrix
#' 
#' @details
#' The smooting matrix can be parametrized using the d eigenvalues and d(d-1)/2
#' angles to represent the orthogonal matrix when carrying out a spectral
#' decomposition
#' 
#' @param d number of columns
#' @param params a vector of d(d+1)/2 parameters. First d are the eigenvalues and
#'          remaining d(d-1)/2 are the angles.
#' @noRd
.smoothing.matrix <- function(params, d) {
  end <- d*(d+1)/2
  theta.lambda <- params[1:d]; theta.angles <- params[(d+1):end]
 
  lambdas <- theta.lambda; angles <- theta.angles 
  W <- .orthogonal.matrix(d, angles)
  N <- W %*% diag(lambdas) %*% t(W)
  
  N
}

#' Huber function
#' 
#' @details
#' Calculate the output of a huber function given the value and the number of 
#' variables
#' 
#' @param x the value
#' @param p number of columns
#' @noRd
.huber <- function(x, p){
  k <- sqrt(qchisq(0.95, df=p))
  min(k, max(x, -k))
}

#' Biweight Rho function
#' 
#' @details
#' Calculate the output of a biweight rho function given the value and the number of 
#' variables
#' 
#' @param x the value
#' @param p number of columns
#' @noRd
.biweight <- function(x, p){
  
  c <- sqrt(qchisq(0.95, df=p))
  g <- .gamma_cp(p, c)
  
  if(abs(x) <= c) g*(1 - (1 - (x/c)^2)^3)
  else g
}

#' Internal helper function for computing the Mahalanobis distances
#' @noRd
.calculate_mahalanobis <- function(cleanValues, covts, dates) {
    xts(
        sapply(dates, function(d) {
            tryCatch({
                x <- cleanValues[d]
                sigma <- covts[d]
                sqrt(x %*% solve(sigma) %*% t(x))
            }, error = function(e) NA)
        }),
        order.by=dates
    )
}

#' Internal helper function for computing the eigenvalues
#' @noRd
.calculate_eigenvalues <- function(matrices, corr = FALSE, dates) {
    xts(
        t(sapply(dates, function(d) {
            tryCatch({
                sigma <- if(corr) matrices$corrts[d] else matrices$covts[d]
                sort(eigen(sigma)$values, decreasing = TRUE)
            }, error = function(e) rep(NA, ncol(sigma)))
        })),
        order.by=dates
    )
}

#' Internal helper function for computing the standard deviation
#' @noRd
.calculate_standard_deviations <- function(matrices, dates) {
    xts(
        t(sapply(dates, function(d) {
            tryCatch({
                sigma <- matrices$covts[d]
                sqrt(diag(sigma))
            }, error = function(e) rep(NA, ncol(sigma)))
        })),
        order.by=dates
    )
}

#' Internal helper function for the Objective
#' @noRd
.obj.helper <- function(smoothing.matrix, R, y.hat, Sigma.hat, startup_period, 
                        training_period, lambda) {
  
  d <- ncol(R); I <- diag(rep(1, d))
  
  forecast.error <- matrix(NA, nrow = training_period - startup_period, ncol = d)
  cleaned.val <- matrix(NA, nrow = training_period - startup_period, ncol = d)
  covMatList <- vector("list", length = training_period - startup_period)

  
  for(t in (startup_period + 1):training_period) {
    
    forecast.error[t - startup_period,] <- R[t,] - y.hat[t-1,]
    ferr <- t(forecast.error[t - startup_period,,drop=FALSE])
    
    I.Sigma.hat <- solve(Sigma.hat)
    temp <- as.numeric(t(ferr) %*% I.Sigma.hat %*% ferr)
    Sigma.hat <- lambda*.biweight(sqrt(temp), d)/temp * ferr %*% t(ferr) + 
      (1 - lambda)*Sigma.hat
    
    temp <- as.numeric(sqrt(t(ferr) %*% solve(Sigma.hat) %*% ferr))
    cleanVal <- cleaned.val[t - startup_period,] <- .huber(temp, d)/temp*ferr + 
      t(y.hat[t-1,,drop=FALSE])
    
    y.hat[t,] <- smoothing.matrix %*% cleanVal + 
      (I - smoothing.matrix) %*% t(y.hat[t-1,,drop=FALSE])
    
    if ((t - startup_period) >= 2*d) {
      h <- floor(0.75 * (t - startup_period))
      forecast.error.t <- forecast.error[1:(t - startup_period), , drop = FALSE]
      covMatList[[t - startup_period]] <- covMcd(forecast.error.t, nsamp = h)$cov
      rownames(covMatList[[t - startup_period]]) <- colnames(covMatList[[t - startup_period]]) <- colnames(R)
    }
  }
  
  time_periods <- index(R)[(startup_period + 1):training_period]
  names(covMatList) <- time_periods
  covMatList <- Filter(Negate(is.null), covMatList)
  smoothVal <- y.hat[(startup_period + 1):training_period,]
  list(smoothValues = smoothVal, cleanValues = cleaned.val, covMatList = covMatList)
}

#' Objective to calculate the optimal smoothing matrix
#' 
#' @details
#' The objective calculates the determinant of the one step ahead forecast errors.
#' The objective is a very noisy function and hence the optimization is time 
#' consuming and replicating paramaters across runs is difficult
#' 
#' @param params d(d+1)/2 paramaters to be optimzed
#' @param R data
#' @param y.hat fitted yalues for the data
#' @param Sigma.hat value of the covarianve matrix
#' @param startup_period length of samples required to calculate initial values
#' @param training_period length of samples required to calculate forecast errors
#'                for evalualating the objective
#' @param lambda known constant as described in the paper
#' @noRd
.obj <- function(params, R, y.hat, Sigma.hat, startup_period, training_period, 
                 lambda) {
  
  d <- ncol(R); smoothing.matrix <- .smoothing.matrix(params, d);
  
  fit <- .obj.helper(smoothing.matrix, R, y.hat, Sigma.hat, startup_period, 
              training_period, lambda)

  last_element <- tail(fit$covMatList, n = 1)[[1]]
  forecast.accuracy <- det(last_element)
  forecast.accuracy
}


#' Optimal Smoothing Matrix
#' 
#' @details
#' Calcuation of smoothing matrix is done by assuming that the smoothing matrix
#' is symmetrix and has a spectral decomposition. The orthogonal matrix in the 
#' decomposition is calculated using the product of givens rotation matrices and
#' requires d(d-1)/2 angles for a d dimensional matrix. The eigenvalues are 
#' restricted to lie in [0,1].
#' 
#' @importFrom optimx optimx
#' @importFrom parallel parRapply parSapply
#' @importFrom RobStatTM covRobRocke
#' 
#' @param R data
#' @param startup_period length of samples required to calculate initial values
#' @param training_period length of samples required to calculate forecast errors
#'                for evalualating the objective
#' @param seed random seed to replicate the starting values for optimization
#' @param trials number of strarting values to try for any optimization. 
#'            Large number of trials for high dimensions can be time consuming
#' @param method optimization method to use to evaluate an estimate of 
#'            smoothing matrix. Default is L-BFGS-B
#' @param lambda known constant as described in the paper. Defaulted to 0.2           
#' 
#' @export
#' @author Rohit Arora
#' 
smoothing.matrix <- function(R, startup_period = 10, training_period = 60 , 
                             seed = 9999, trials = 50, method = "L-BFGS-B",
                             lambda = 0.2) { 
  
  M <- nrow(R); d <- ncol(R)
  if(M < 4*d) stop("Not enough data for estimation")
  
  if(is.na(startup_period) || startup_period < 2*d) startup_period <- 2*d
  
  if(is.na(training_period) || training_period < (startup_period + 2*d)) 
    training_period <- startup_period + 2*d
  
  if ( M < (startup_period + training_period)) 
    stop("Insufficient data. Reset correct startup & training periods")
  
  startup.fit <- lapply(1:d, function(i) {
    lmrob(coredata(R[1:startup_period,i]) ~ as.matrix(1:startup_period))
  })
  
  y.hat <- matrix(NA, nrow = training_period, ncol = ncol(R))
  y.hat[1:startup_period,] <- do.call(cbind, lapply(startup.fit, fitted))
  
  res <- do.call(cbind, lapply(startup.fit, residuals))  
  Sigma.hat <- RobStatTM::covRobRocke(res)$cov
  
  set.seed(seed)
  
  lower <- c(rep(0,d), rep(-pi/2, d*(d-1)/2))
  upper <- c(rep(1,d), rep(pi/2, d*(d-1)/2))
  nlower <- length(lower); width <- upper - lower
  
  Umin <- matrix(rep.int(lower, trials), nrow = trials, ncol=nlower, byrow=T)
  start <- (Umin + matrix(rep.int(width, trials), nrow = trials, 
                          ncol=nlower, byrow=T)*lhs::maximinLHS(n = trials, k = nlower))
  
  # Convert rows of 'start' matrix to a list of parameter vectors
  start_list <- split(start, row(start))
  
  objmin <- lapply(start_list, function (x)
    try(optim(x, .obj, lower = lower, upper = upper, method = method, R = R , 
              y.hat = y.hat, Sigma.hat = Sigma.hat, startup_period = startup_period, 
              training_period = training_period, lambda = lambda), silent=FALSE))
  
  fit <- objmin[[which.min(sapply(objmin, '[[', "value"))]]
  
  params <- unlist(fit$par[1:(d*(d+1)/2)])
  N <- .smoothing.matrix(params, d)
  list(smooth.mat = N, minObj = fit, startup_period = startup_period, 
       training_period = training_period)
}



#' Robust Multivariate Exponential Smoothing
#' 
#' @description
#' Calculate a robust estimate of the covariance matrix while also smoothing and 
#' cleaning the data using the procedure described in 
#' (Croux, Gelper, and Mahieu, 2010).
#' 
#' @importFrom RobStatTM covRobRocke
#' 
#' @param R An xts object containing the data.
#' @param corr A boolean indicating whether additional correlation matrix is needed. If true,
#'              additional correlation matrix is returned.
#' @param smoothMat Optional. The optimal smoothing matrix. If missing, it is estimated. 
#'                  The procedure may be very slow for high-dimensional data. Also,
#'                  the objective function being very noisy, optimization across
#'                  multiple runs may lead to different smoothing matrices.
#' @param startup_period Integer. Length of samples required to calculate initial values.
#' @param training_period Integer. Length of samples required to calculate forecast errors
#'                        for evaluating the objective if the smoothing matrix is estimated.
#' @param seed Integer. Random seed to replicate the starting values for optimization.
#' @param trials Integer. Number of starting values to try for any optimization. 
#'               A large number of trials for high dimensions can be time-consuming.
#' @param method Character. Optimization method to use to evaluate an estimate of 
#'               the smoothing matrix. Default is "L-BFGS-B".
#' @param lambda Numeric. Known constant as described in the paper. Defaults to 0.2.
#' 
#' @return A list containing:
#' \item{smoothMat}{The optimal smoothing matrix.}
#' \item{smoothValues}{The smoothed values as an xts object.}
#' \item{cleanValues}{The cleaned values as an xts object.}
#' \item{covMatList}{A list of covariance matrices for each time period.}
#' \item{mahalanobis_distances}{A list of distances using covariance for each time period.}
#' \item{eigenvalues}{A list of eigenvalues for covariance/correlation for each time period.}
#' \item{standard_deviations}{A list of standard deviations when corr is TRUE for each time period.}
#' 
#' @export
#' @author Rohit Arora
robustMultExpSmoothing <- function(R, corr = FALSE, smoothMat = NA, startup_period = 10, 
                         training_period = 60 , seed = 9999, trials = 50, 
                         method = "L-BFGS-B", lambda = 0.2) {
  
  if(!is.xts(R)) stop("Only xts objects are allowed")
  
  M <- nrow(R); d <- ncol(R)
  
  if (is.null(dim(smoothMat))) {
    
    if (is.na(trials) || is.na(method) || is.na(lambda))
      stop("Incorrect parameters for smoothing matrix estimation")
    
    smoothfit <-  smoothing.matrix(R, startup_period = startup_period, 
                                   training_period = training_period, 
                                   seed = seed, trials = trials, 
                                   method = method, lambda = lambda) 
    
    smoothMat <- smoothfit$smooth.mat
    startup_period <- smoothfit$startup_period
    training_period <- smoothfit$training_period
  }
  else {
    if(!isSymmetric(smoothMat)) stop("Smoothing matrix must be symmetric")
    
    eigVal <- eigen(smoothMat, symmetric = TRUE, only.values = TRUE)$values
    if (any((eigVal < 0 | eigVal > 1))) 
      stop("Eigenvalues of smoothing matrix must be in [0,1]")
    
    if(is.na(startup_period) || startup_period < 2*d) startup_period <- 2*d
  }
  
  y.hat <- matrix(NA, nrow = M, ncol = ncol(R))
  startup.fit <- lapply(1:d, function(i) {
    lmrob(coredata(R[1:startup_period,i]) ~ as.matrix(1:startup_period))
  })
  
  y.hat[1:startup_period,] <- do.call(cbind, lapply(startup.fit, fitted))
  
  res <- do.call(cbind, lapply(startup.fit, residuals))  
  Sigma.hat <- RobStatTM::covRobRocke(res)$cov
  
  fit <- .obj.helper(smoothMat, R = R, y.hat = y.hat, Sigma.hat = Sigma.hat, 
                     startup_period = startup_period, 
                     training_period = M, lambda = lambda)
  
  smoothValues <- fit$smoothValues; cleanValues <- fit$cleanValues; covMatList <- fit$covMatList  
  smoothValues <- xts(smoothValues, order.by = index(R)[(startup_period+1):M])
  cleanValues <- xts(cleanValues, order.by = index(R)[(startup_period+1):M])

  covts <- create_matrix_xts(covMatList)

  result <- list(
    smoothMat = smoothMat,
    smoothValues = smoothValues,
    cleanValues = cleanValues,
    covts = covts
  )
  
  if (corr) {
    result$corrts <- create_matrix_xts(lapply(covMatList, cov2cor))
  }

  dates <- as.Date(names(covMatList))
  result$mahalanobis_distances <- .calculate_mahalanobis(cleanValues, result$covts, dates)
  result$eigenvalues <- .calculate_eigenvalues(result, corr, dates)
  result$standard_deviations <- .calculate_standard_deviations(result, dates)
  
  return(result)
}

#' Plot Changes in Correlation/Covariance Matrices
#' 
#' @description
#' Creates a heatmap visualization of percentage changes between adjacent time periods
#' for correlation or covariance matrices. For correlation matrices, only the upper
#' triangular elements (excluding diagonal) are shown. For covariance matrices, the
#' diagonal elements are included.
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal labs theme
#' 
#' @param smoothfit A list containing 'corrts' and 'covts' elements from robustMultExpSmoothing
#'                  output, each being a MatrixXTS object of correlation/covariance matrices
#' @param corr Logical. If TRUE, plots correlation changes (excluding diagonal).
#'             If FALSE, plots covariance changes including diagonal elements.
#' @param start_date Character. Start date in "YYYY-MM-DD" format. If NULL, starts
#'                   from first available date.
#' @param end_date Character. End date in "YYYY-MM-DD" format. If NULL, ends at
#'                 last available date.
#' @param symbols Character vector. Asset symbols to include. If NULL, uses all
#'                available symbols.
#' 
#' @return A ggplot object showing a heatmap visualization with:
#' \item{x-axis}{Time periods}
#' \item{y-axis}{Matrix entries (labeled by symbol pairs if symbols provided)}
#' \item{colors}{Green (positive changes), Red (negative changes), White (minimal changes)}
#' 
#' @examples
#' # Plot correlation changes for specific symbols and date range
#' plot_matrix_changes(smoothfit, corr = TRUE, 
#'                    start_date = "2014-07-07", 
#'                    end_date = "2014-08-30",
#'                    symbols = c("AAPL", "MSFT"))
#'                    
#' # Plot full covariance matrix changes
#' plot_matrix_changes(smoothfit, corr = FALSE)
#' 
#' @export
plot_matrix_changes <- function(smoothfit, corr = FALSE, 
                                start_date = NULL, end_date = NULL,
                                symbols = NULL) {
  
  data_to_use <- if(corr) {
    if (is.null(smoothfit$corrts)) {
      stop("Correlation matrices not available in smoothfit object")
    }
    smoothfit$corrts
  } else {
    if (is.null(smoothfit$covts)) {
      stop("Covariance matrices not available in smoothfit object")
    }
    smoothfit$covts
  }

  
  # Get matrices for the date range
  if (!is.null(start_date) && !is.null(end_date)) {
    date_range <- paste0(start_date, "/", end_date)
    matrices <- data_to_use[date_range]
  } else {
    matrices <- data_to_use[index(data_to_use)]
  }
  
  # Filter by symbols if provided
  if (!is.null(symbols)) {
    all_symbols <- colnames(matrices[[1]])
    symbol_indices <- which(all_symbols %in% symbols)
    if (length(symbol_indices) > 0) {
      matrices <- lapply(matrices, function(mat) {
        mat[symbol_indices, symbol_indices, drop = FALSE]
      })
    }
  }
  
  # Function to extract upper triangle (excluding diagonal for correlations)
  get_upper_triangle <- function(mat) {
    n <- nrow(mat)
    indices <- which(upper.tri(mat, diag = !corr), arr.ind = TRUE)
    mat[indices]
  }
  
  # Convert list of matrices to array for easier manipulation
  dates <- names(matrices)
  n_entries <- ncol(matrices[[1]])
  n_upper <- if(corr) {
    (n_entries * (n_entries - 1)) / 2  # excluding diagonal
  } else {
    (n_entries * (n_entries + 1)) / 2  # including diagonal
  }
  
  flat_data <- do.call(rbind, lapply(matrices, get_upper_triangle))
  
  # Calculate percentage differences
  pct_diff <- 100 * (flat_data[-1,] - flat_data[-nrow(flat_data),]) / 
    flat_data[-nrow(flat_data),]
  
  # Cap outliers
  cap_outliers <- function(x) {
    lower_cap <- quantile(x, 0.01, na.rm = TRUE)
    upper_cap <- quantile(x, 0.99, na.rm = TRUE)
    pmin(pmax(x, lower_cap), upper_cap)
  }
  
  capped_pct_diff <- as.vector(t(apply(pct_diff, 2, cap_outliers)))
  
  # Create data frame for plotting
  df <- data.frame(
    Entry = rep(1:n_upper, nrow(pct_diff)),
    Time = rep(as.Date(dates[-1]), each = n_upper),
    Value = capped_pct_diff
  )
  
  # Create entry labels if symbols provided
  if (!is.null(symbols)) {
    n <- length(symbols)
    indices <- which(upper.tri(matrix(0, n, n), diag = !corr), arr.ind = TRUE)
    entry_labels <- paste(symbols[indices[,1]], symbols[indices[,2]])
  }
  
  # Create plot
  p <- ggplot(df, aes(x = Time, y = Entry, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "darkgreen",
      midpoint = 0,
      name = "% Change",
      limits = c(-100, 100)
    ) +
    theme_minimal() +
    labs(
      title = paste("Percentage Changes Between Adjacent Time Periods -",
                    if(corr) "Correlations" else "Covariances"),
      x = "Time",
      y = "Matrix Entry"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "gray", size = 0.5),
      panel.grid.minor = element_line(color = "lightgray", size = 0.25)
    )
  
  # Add appropriate y-axis labels
  if (!is.null(symbols)) {
    p <- p + scale_y_continuous(
      breaks = 1:length(entry_labels),
      labels = entry_labels,
      minor_breaks = seq(0.5, length(entry_labels) + 0.5, by = 1)
    )
  } else {
    p <- p + scale_y_continuous(
      breaks = 1:n_upper,
      minor_breaks = seq(0.5, n_upper + 0.5, by = 1)
    )
  }
  
  return(p)
}

