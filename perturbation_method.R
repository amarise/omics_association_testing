# Load required libraries
library(GENESIS)  # For genetic analysis and null model fitting
library(SKAT)     # For kernel-based association testing
library(CompQuadForm)  # For computing Davies' method and saddlepoint approximation

# Function for performing perturbation-based association testing for two omics data types
perturbation <- function(w, y, X, omics1, omics2, grm, num_reps = 10000) {
  
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

  # Perturbation method begins
  pert_Lambda <- list()
  bigV <- NULL
  
  for (i in 1:length(w)) {
    K <- w[i]*K1 + (1-w[i])*K2
    Q <-  as.numeric(crossprod(nullmod$fit$resid.PY, K) %*% nullmod$fit$resid.PY)
    V <- crossprod(P0, K) %*% P0
    Veig <- eigen(V, symmetric = T)
    wK <- which(Veig$values > 1e-10)
    pert_Lambda[[i]] <- Veig$values[wK]
    bigV <- rbind(bigV, t(Veig$vectors[,wK]))
    
    # Use Davies' method to compute the p-value
    f <- CompQuadForm::davies(q=Q, lambda=Veig$values[wK], acc = 1e-9)
    
    # Handle potential issues with Davies' method and switch to saddlepoint approximation if necessary
    if((f$ifault > 0) | (f$Qq < 1e-9) | (f$Qq > 1)) {
      f$Qq <- survey:::saddle(Q, Veig$values[wK])
      if(is.na(f$Qq)) f$Qq <- 1
    }
    
    grid_results[i, 'test.stat'] <- Q
    grid_results[i, 'test.pval'] <- f$Qq
  }

  # Perturbation steps for generating random data and recalculating p-values
  bigSigma <- tcrossprod(bigV)
  rm(bigV)
  bigSigmaEig <- eigen(bigSigma, symmetric = T)
  wK <- which(bigSigmaEig$values < 1e-10)
  bigSigmaEig$values[wK] <- 0
  R <- bigSigmaEig$vectors %*% diag(sqrt(bigSigmaEig$values))
  
  min_p_stars <- rep(NA, num_reps + 1)
  min_p_stars[1] <- as.numeric(min(grid_results$test.pval))
  
  start <- cumsum(c(1, unlist(lapply(pert_Lambda, length))[-length(pert_Lambda)]))
  end <- cumsum(unlist(lapply(pert_Lambda, length)))
  
  for (j in 1:num_reps) {
    r <- rnorm(ncol(R))
    r_star <- R %*% r

    p_star_min <- Inf
    for (k in 1:length(w)) {
      r_star_d <- r_star[start[k]:end[k]]
      if (length(r_star_d) == 1) {
        Q_star_d <- r_star_d^2 * pert_Lambda[[k]]
      } else {
        Q_star_d <- crossprod(r_star_d, diag(pert_Lambda[[k]])) %*% r_star_d
      }
      
      f <- CompQuadForm::davies(q=Q_star_d, lambda=pert_Lambda[[k]], acc = 1e-9)
      if((f$ifault > 0) | (f$Qq < 1e-9) | (f$Qq > 1)) {
        f$Qq <- survey:::saddle(Q_star_d, pert_Lambda[[k]])
        if(is.na(f$Qq)) f$Qq <- 1
      }
      if (f$Qq < p_star_min) p_star_min <- f$Qq
    }
    min_p_stars[j+1] <- p_star_min
  }

  list(grid_results = grid_results,
       p_val = mean(min_p_stars <= min(grid_results$test.pval)),
       test.stat = grid_results$test.stat[which.min(grid_results$test.pval)],
       w = grid_results$w[which.min(grid_results$test.pval)])
}

# Helper function to compute standard deviation for quadratic forms
sd_quad <- function(M) {
  sqrt(2 * sum(M * M))
}
