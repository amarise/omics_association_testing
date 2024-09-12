# Load required libraries
library(GENESIS)  # For genetic analysis and null model fitting
library(SKAT)     # For kernel-based association testing
library(CompQuadForm)  # For computing Davies' method and saddlepoint approximation

# Main function for performing kernel PCA association testing for two omics data types
# Inputs:
#   w: Vector of weights used to balance kernels for omics1 and omics2
#   y: Outcome variable (phenotype or trait)
#   X: Covariate matrix (Samples in rows, covariates in columns)
#   omics1: First omics data matrix (Samples in rows, features in columns)
#   omics2: Second omics data matrix (Samples in rows, features in columns)
#   grm: Genetic relatedness matrix (covariance structure for related individuals)
# Outputs:
#   Returns a list with:
#     - grid_results: Data frame containing test statistics and p-values for each weight in w
#     - p_val: Overall p-value for the test
#     - test.stat: Test statistic corresponding to the optimal weight
#     - w: Optimal weight for kernel combination
kernelPCA <- function(w, y, X, omics1, omics2, grm) {
  
  # Create a data frame to store results for different weights (w)
  grid_results <- data.frame(w,
                             test.stat = numeric(length(w)),  # Test statistics
                             test.pval = numeric(length(w)))  # P-values
  
  # Fit the null model (without including the omics data)
  nullmod <- fitNullModel(data.frame(y, X, row.names = rownames(grm)),
                          outcome = 'y', covars = c('X'), cov.mat = grm,  # Covariate matrix for relatedness
                          family = "gaussian", start = NULL, drop.zeros = F,
                          return.small = F, verbose = F)
  
  # Construct projection matrix (P0)
  P0 <- -tcrossprod(nullmod$CXCXI,nullmod$CX); diag(P0) <- 1 + diag(P0)
  P0 <- nullmod$cholSigmaInv %*% P0

  # Calculate kernel matrices for omics1 and omics2
  K1 <- tcrossprod(omics1); K2 <- tcrossprod(omics2)
  
  # Compute standard deviation for quadratic forms of omics1 and omics2
  sdQOmics1 <- sd_quad(tcrossprod(crossprod(omics1,P0)))
  sdQOmics2 <- sd_quad(tcrossprod(crossprod(omics2,P0)))
  
  # Compute eta (balance between omics1 and omics2)
  eta <- sdQOmics1 / (sdQOmics1 + sdQOmics2)
  
  # Adjust kernel matrices with eta
  K1 <- (1  - eta) * K1
  K2 <- eta * K2
  
  # Eigen decomposition of kernel matrices
  eigenK1 <- eigen(K1)
  eigenK2 <- eigen(K2)
  
  # Filter out small eigenvalues for omics1
  wK <- which(eigenK1$values > 1e-9)
  Z1 <- eigenK1$vectors[,wK] %*% diag(sqrt(eigenK1$values[wK]))
  
  # Filter out small eigenvalues for omics2
  wK <- which(eigenK2$values > 1e-9)
  Z2 <- eigenK2$vectors[,wK] %*% diag(sqrt(eigenK2$values[wK]))
  
  # Projection of Z1 and Z2 to adjust for covariates
  temp_proj <- -tcrossprod(tcrossprod(nullmod$model.matrix, nullmod$CXCXI), nullmod$cholSigmaInv)
  diag(temp_proj) <- 1 + diag(temp_proj)
  Z1_star <- temp_proj %*% Z1
  
  # Projection of Z2 based on Z1 (double projection)
  proj_part <- crossprod(Z1_star, tcrossprod(nullmod$cholSigmaInv))
  Pg0 <- try(-Z1_star %*% solve(proj_part %*% Z1_star) %*% proj_part)

  if (class(Pg0) == "try-error") {
    # Handle error in projection calculation for Z2
    print(paste0("Calculating Z2 projection without double projection at ", Sys.time()))
    proj_part <- crossprod(Z1, tcrossprod(nullmod$cholSigmaInv))
    Pg0 <- -Z1 %*% solve(proj_part %*% Z1) %*% proj_part
  }
  
  diag(Pg0) <- 1 + diag(Pg0)
  Z2_star <- Pg0 %*% Z2

  # Update eta after projection and recompute Z1 and Z2_star
  old_eta <- eta
  sdQOmics1 <- sdQOmics1 * (1 - old_eta)
  sdQOmics2 <- sd_quad(tcrossprod(crossprod(Z2_star,P0)))
  eta <- sdQOmics1 / (sdQOmics1 + sdQOmics2)
  
  Z1 <- sqrt(1 - eta) * Z1
  Z2_star <- sqrt(eta) * Z2_star

  # Compute projected Z matrices for omics1 and omics2
  Z1.1 <- crossprod(P0, Z1)
  Z2.1 <- crossprod(P0, Z2_star)

  # Get test statistics (Q) and p-values for each weight in grid (w)
  res.out = NULL
  out.Q <- SKAT_2Kernel_Optimal_Get_Q(Z1, Z2_star, nullmod$fit$resid.PY, w, n.Resampling = 0 , res.out)
  
  # Collect p-values and test statistics for each kernel
  Q.all <- rbind(out.Q$Q.r, out.Q$Q.r.res)
  out <- SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1 / sqrt(2), Z2.1 / sqrt(2), w)

  # Store results and return the best result
  grid_results$test.pval[length(w):1] <- out$p.val.each
  grid_results$test.stat[length(w):1] <- Q.all[1,]
  index <- which.min(grid_results$test.pval)
  
  # Return the result with the smallest p-value
  list(grid_results = grid_results,
       p_val = out$p.value,
       test.stat = grid_results$test.stat[index],
       w = grid_results$w[index])
}