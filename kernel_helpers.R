# Helper function to compute test statistics for a range of kernel weights
SKAT_2Kernel_Optimal_Get_Q <- function(Z1, Z2, res, r.all, n.Resampling = 0, res.out = NULL, res.moments = NULL) {
  n.r <- length(r.all)
  p.m <- dim(Z1)[2]

  Q.r <- rep(0, n.r)
  Q.r.res <- NULL
  Q.sim <- NULL

  temp1 <- t(res) %*% Z1
  temp2 <- t(res) %*% Z2

  # Calculate Q-statistics for all weights
  for (i in 1:n.r) {
    r.corr <- r.all[i]
    Q1 <- (1 - r.corr) * rowSums(temp1^2)
    Q2 <- r.corr * rowSums(temp2^2)
    Q.r[i] <- Q1 + Q2
  }
  Q.r <- Q.r / 2

  # For resampling, calculate Q-statistics for resampled data
  if (n.Resampling > 0) {
    temp1 <- t(res.out) %*% Z1
    temp2 <- t(res.out) %*% Z2
    Q.r.res <- matrix(rep(0, n.Resampling * n.r), ncol = n.r)
    for (i in 1:n.r) {
      r.corr <- r.all[i]
      Q1 <- (1 - r.corr) * rowSums(temp1^2)
      Q2 <- r.corr * rowSums(temp2^2)
      Q.r.res[, i] <- Q1 + Q2
    }
    Q.r.res <- Q.r.res / 2
  }

  # If moments are available, calculate Q-statistics for moments
  if (!is.null(res.moments)) {
    temp1 <- t(res.moments) %*% Z1
    temp2 <- t(res.moments) %*% Z2
    n.moments <- dim(res.moments)[2]
    Q.sim <- matrix(rep(0, n.moments * n.r), ncol = n.r)
    for (i in 1:n.r) {
      r.corr <- r.all[i]
      Q1 <- (1 - r.corr) * rowSums(temp1^2)
      Q2 <- r.corr * rowSums(temp2^2)
      Q.sim[, i] <- Q1 + Q2
    }
    Q.sim <- Q.sim / 2
  }

  re <- list(Q.r = Q.r, Q.r.res = Q.r.res, Q.sim = Q.sim)
  return(re)
}

# Helper function to calculate p-values from kernel test statistics
SKAT_2Kernel_Ortho_Optimal_Get_Pvalue <- function(Q.all, Z1.1, Z2.1, r.all) {
  n.r <- length(r.all)
  n.q <- dim(Q.all)[1]
  p.m <- dim(Z1.1)[2]

  # Get cumulants needed for method of moments for all grid values
  c1.all <- SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r(Z1.1, Z2.1, r.all)

  # Get information for each Q-statistic and calculate p-values
  Each_Info <- SKAT_2Kernel_Ortho_Optimal_Each_Q(Q.all, r.all, c1.all)
  pmin.q <- Each_Info$pmin.q
  pmin <- Each_Info$pmin
  pval <- rep(0, n.q)

  multi <- 3
  if (length(r.all) < 3) {
    multi <- 2
  }

  # If p-values are negative, adjust and use minimum p-value
  for (i in 1:n.q) {
    pval.each <- Each_Info$pval[i, ]
    IDX <- which(pval.each > 0)

    pval1 <- min(pval.each) * multi
    if (pval[i] <= 0 || length(IDX) < length(r.all)) {
      pval[i] <- pval1
    }

    if (pval[i] == 0 && length(IDX) > 0) {
      pval[i] <- min(pval.each[IDX])
    }
  }

  return(list(p.value = pval, p.val.each = Each_Info$pval))
}

# Helper function to get parameters for method of moments for each weight (r)
SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r <- function(Z1.1, Z2.1, r.all) {
  c1 <- matrix(rep(0, 4 * length(r.all)), ncol = length(r.all))

  A1 <- crossprod(Z1.1)  # Zg' P0 %*% P0' Zg
  B1 <- crossprod(Z2.1)

  A2 <- A1 %*% A1
  B2 <- B1 %*% B1

  A11 <- crossprod(Z1.1, Z2.1)
  A22 <- tcrossprod(A11)
  B22 <- crossprod(A11)
  B333 <- crossprod(A11, A1) %*% A11

  # Compute cumulants
  c1[1, ] <- sum(diag(A1)) * (1 - r.all) + sum(diag(B1)) * r.all
  c1[2, ] <- sum(diag(A2)) * (1 - r.all)^2 + sum(diag(B2)) * r.all^2 + sum(diag(A22)) * 2 * (1 - r.all) * r.all
  c1[3, ] <- sum(A2 * A1) * (1 - r.all)^3 + sum(B2 * B1) * r.all^3 + sum(A22 * A1) * 3 * (1 - r.all)^2 * r.all + sum(B1 * B22) * 3 * (1 - r.all) * r.all^2
  c1[4, ] <- sum(A2 * A2) * (1 - r.all)^4 + sum(B2 * B2) * r.all^4 + sum(A22 * A2) * 4 * (1 - r.all)^3 * r.all + sum(B2 * B22) * 4 * (1 - r.all) * r.all^3 + sum(B1 * B333) * 4 * (1 - r.all)^2 * r.all^2 + sum(B22 * B22) * 2 * (1 - r.all)^2 * r.all^2

  return(c1)
}

# Helper function to evaluate each Q-statistic and return p-values
SKAT_2Kernel_Ortho_Optimal_Each_Q <- function(Q.all, r.all, c1.all) {
  n.r <- length(r.all)
  n.q <- nrow(Q.all)

  pval <- matrix(rep(0, n.r * n.q), ncol = n.r)
  pmin.q <- matrix(rep(0, n.r * n.q), ncol = n.r)
  param.mat <- NULL

  # Compute p-values for each Q-statistic
  for (i in 1:n.r) {
    Q <- Q.all[, i]
    r.corr <- r.all[i]

    c1 <- c1.all[, i]
    param.temp <- SKAT:::Get_Liu_Params_Mod(c1)

    muQ <- param.temp$muQ
    varQ <- param.temp$sigmaQ^2
    df <- param.temp$l

    # Get p-value using method of moments
    Q.Norm <- (Q - muQ) / sqrt(varQ) * sqrt(2 * df) + df
    pval[, i] <- pchisq(Q.Norm, df = df, lower.tail = FALSE)

    param.mat <- rbind(param.mat, c(muQ, varQ, df))
  }

  # Get quantiles for each mixture
  pmin <- apply(pval, 1, min)
  for (i in 1:n.r) {
    muQ <- param.mat[i, 1]
    varQ <- param.mat[i, 2]
    df <- param.mat[i, 3]

    q.org <- qchisq(1 - pmin, df = df)
    q.q <- (q.org - df) / sqrt(2 * df) * sqrt(varQ) + muQ
    pmin.q[, i] <- q.q
  }

  param.m <- list(par.moments = list(param.mat[1, ], param.mat[n.r, ]))

  out <- list(pmin = pmin, pval = pval, pmin.q = pmin.q, param.m = param.m)
  return(out)
}
