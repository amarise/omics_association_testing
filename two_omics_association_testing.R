library(GENESIS)
library(SKAT)

kernelPCA <- function(w, y, X, G, M, grm)
{
  grid_results <- data.frame(w,
                              test.stat = numeric(length(w)),
                              test.pval = numeric(length(w)))
  nullmod <- fitNullModel(data.frame(y, X, row.names = rownames(grm)),
                          outcome = 'y', covars = c('X'), cov.mat = grm,
                          family = "gaussian", start = NULL, drop.zeros = F,
                          return.small = F, verbose = F)
  P0 <- -tcrossprod(nullmod$CXCXI,nullmod$CX); diag(P0) <- 1 + diag(P0);
  P0 <- nullmod$cholSigmaInv %*% P0

  K1 <- tcrossprod(G); K2 <- tcrossprod(M);
  sdQG <- sd_quad(tcrossprod(crossprod(G,P0)));
  sdQM <- sd_quad(tcrossprod(crossprod(M,P0)));
  eta <- sdQG / (sdQG + sdQM);
  K1 <- (1  - eta) * K1; K2 <- eta * K2;
  eigenK1 <- eigen(K1); eigenK2 <- eigen(K2);
  wK <- which(eigenK1$values > 1e-9)
  Zg <- eigenK1$vectors[,wK] %*% diag(sqrt(eigenK1$values[wK]))
  wK <- which(eigenK2$values > 1e-9)
  Zm <- eigenK2$vectors[,wK] %*% diag(sqrt(eigenK2$values[wK]))

  temp_proj <- -tcrossprod(tcrossprod(nullmod$model.matrix, nullmod$CXCXI), nullmod$cholSigmaInv); diag(temp_proj) <- 1+diag(temp_proj)
  Zg_star <- temp_proj %*% Zg
  proj_part <- crossprod(Zg_star, tcrossprod(nullmod$cholSigmaInv))
  Pg0 <- try(-Zg_star %*% solve(proj_part %*% Zg_star) %*% proj_part);
  if (class(Pg0) == "try-error")
  {
    print(paste0("Calculating Zm projection without double projection at ",Sys.time()))
    proj_part <- crossprod(Zg, tcrossprod(nullmod$cholSigmaInv)) #gxn %*% nxn = gxn
    Pg0 <- -Zg %*% solve(proj_part %*% Zg) %*% proj_part #nxg %*% gxg = nxg
  }
  diag(Pg0) <- 1+diag(Pg0)
  Zm_star <- Pg0 %*% Zm

  old_eta <- eta
  sdQG <- sdQG * (1 - old_eta)
  sdQM <- sd_quad(tcrossprod(crossprod(Zm_star,P0)));
  eta <- sdQG / (sdQG + sdQM);
  Zg <- sqrt(1 - eta) * Zg
  Zm_star <- sqrt(eta) * Zm_star

  Z1.1 <- crossprod(P0, Zg)
  Z2.1 <- crossprod(P0, Zm_star)

  res.out = NULL
  out.Q <- SKAT_2Kernel_Optimal_Get_Q(Zg, Zm_star, nullmod$fit$resid.PY, w, n.Resampling = 0 , res.out)
  Q.all <- rbind(out.Q$Q.r, out.Q$Q.r.res)
  out <- SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1/sqrt(2),
                                               Z2.1/sqrt(2), w)

  grid_results$test.pval[length(w):1] <- out$p.val.each
  grid_results$test.stat[length(w):1] <- Q.all[1,]
  index <- which.min(grid_results$test.pval)
  list(grid_results =grid_results,
       p_val = out$p.value,
       test.stat = grid_results$test.stat[index],
       w = grid_results$w[index])
}

SKAT_2Kernel_Optimal_Get_Q <- function(Z1, Z2, res, r.all, n.Resampling, res.out, res.moments=NULL)
{
  n.r<-length(r.all)
  p.m<-dim(Z1)[2]

  Q.r<-rep(0,n.r)
  Q.r.res<-NULL
  Q.sim<-NULL

  temp1<-t(res) %*% Z1
  temp2<-t(res) %*% Z2
  for(i in 1:n.r){
    r.corr<-r.all[i]
    Q1<-(1-r.corr) * rowSums(temp1^2)
    Q2<-r.corr * rowSums(temp2^2)
    Q.r[i]<-Q1 + Q2
  }
  Q.r = Q.r /2
  if(n.Resampling > 0){

    temp1<-t(res.out) %*% Z1
    temp2<-t(res.out) %*% Z2
    Q.r.res<-matrix(rep(0,n.Resampling *n.r),ncol=n.r)
    for(i in 1:n.r){
      r.corr<-r.all[i]
      Q1<-(1-r.corr) * rowSums(temp1^2)
      Q2<-r.corr * rowSums(temp2^2)
      Q.r.res[,i]<-Q1 + Q2
    }
    Q.r.res = Q.r.res/2
  }

  if(!is.null(res.moments)){

    temp1<-t(res.moments) %*% Z1
    temp2<-t(res.moments) %*% Z2
    n.moments<-dim(res.moments)[2]
    Q.sim<-matrix(rep(0,n.moments *n.r),ncol=n.r)
    for(i in 1:n.r){
      r.corr<-r.all[i]
      Q1<-(1-r.corr) * rowSums(temp1^2)
      Q2<-r.corr * rowSums(temp2^2)
      Q.sim[,i]<-Q1 + Q2
    }
    Q.sim = Q.sim/2

  }

  re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.sim=Q.sim)
  return(re)
}

SKAT_2Kernel_Ortho_Optimal_Get_Pvalue <- function(Q.all, Z1.1, Z2.1, r.all)
{

  n.r <- length(r.all)
  n.q <- dim(Q.all)[1]
  p.m <- dim(Z1.1)[2]

  # Get cumulants needed for method of moments for all grid values
  c1.all <- SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r(Z1.1, Z2.1, r.all)

  # Get Mixture param
  # These are the moment-based parameters needed for method of
  # moments as in Liu2009. Only for K1 and K2_star (no mixtures, yet)
  # This part is unnecessary...it is redundant to what is calculated for
  # Each_Info
  # param.m <- SKAT_2Kernel_Ortho_Optimal_Param(Z1.1, Z2.1,r.all)

  # The param.m input doesn't get used
  # This function gets p-vals for each test stat using MoM
  Each_Info <- SKAT_2Kernel_Ortho_Optimal_Each_Q(Q.all, r.all, c1.all)
  pmin.q <- Each_Info$pmin.q
  pmin <- Each_Info$pmin
  pval <- rep(0,n.q)
  param.m <- Each_Info$param.m
  #pmin1<<-pmin

  multi<-3
  if(length(r.all) < 3){
    multi<-2
  }

  # if any p-vals are negative, then reassign the final p-val to be multiple of min p-value from grid search
  for(i in 1:n.q){
    pval.each <- Each_Info$pval[i,]
    IDX <- which(pval.each > 0)

    pval1 <- min(pval.each) * multi
    if(pval[i] <= 0 || length(IDX) < length(r.all)){
      pval[i] <- pval1
    }

    # if pval==0, use nonzero min each.pval as p-value
    if(pval[i] == 0){
      if(length(IDX) > 0){
        pval[i] = min(pval.each[IDX])
      }
    }
  }
  return(list(p.value=pval,p.val.each=Each_Info$pval))
}

SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r <- function(Z1.1, Z2.1, r.all)
{
  c1 <- matrix(rep(0,4* length(r.all)), ncol=length(r.all))

  A1 <- crossprod(Z1.1) # Zg' P0 %*% P0' Zg
  B1 <- crossprod(Z2.1)

  A2 <- A1 %*% A1
  B2 <- B1 %*% B1

  A11 <- crossprod(Z1.1, Z2.1)
  A22 <- tcrossprod(A11)
  B22 <- crossprod(A11)
  B333 <- crossprod(A11, A1) %*% A11


  #c1[1,]<-sum(Z1.1^2) * (1-r.all) + sum(Z2.1^2) * r.all # (1-w)*tr(Zg' P0 P0 Zg) + w*tr(Zm_star' P0 P0 Zm_star) // tr(Zg' P0 P0 Zg) = sum of eigenvals of Zg' P0 P0 Zg
  c1[1,] <- sum(diag(A1)) * (1-r.all) + sum(diag(B1)) * r.all
  #c1[2,]<-sum(A1^2) * (1-r.all)^2 + sum(B1^2) * (r.all)^2 + sum(A11^2) * 2 * (1-r.all) * r.all # (1-w)*tr(Zg' P0 P0 Zg) + w*tr(Zm_star' P0 P0 Zm_star)
  c1[2,] <- sum(diag(A2)) * (1-r.all)^2 + sum(diag(B2)) * (r.all)^2 + sum(diag(A22)) * 2 * (1-r.all) * r.all
  c1[3,] <- sum(A2 * A1) * (1-r.all)^3 + sum(B2 * B1) * (r.all)^3 + sum(A22 * A1) * 3 * (1-r.all)^2 * r.all + sum(B1 * B22) * 3 * (1-r.all) * r.all^2
  c1[4,] <- sum(A2 * A2) * (1-r.all)^4 + sum(B2 * B2) * (r.all)^4 + sum(A22 * A2) * 4 * (1-r.all)^3 * r.all + sum(B2 * B22) * 4 * (1-r.all) * r.all^3 + sum(B1 * B333) * 4 * (1-r.all)^2 * r.all^2 + sum(B22 * B22) * 2 * (1-r.all)^2 * r.all^2

  return(c1)
}

SKAT_2Kernel_Ortho_Optimal_Each_Q <- function(Q.all, r.all, c1.all)
{
  n.r <- length(r.all)
  n.q <- nrow(Q.all)

  pval <- matrix(rep(0,n.r*n.q),ncol=n.r)
  pmin.q <- matrix(rep(0,n.r*n.q),ncol=n.r)
  param.mat <- NULL

  for(i in 1:n.r){
    Q <- Q.all[,i]
    r.corr <- r.all[i]

    c1 <- c1.all[,i]
    param.temp <- SKAT:::Get_Liu_Params_Mod(c1)

    muQ <- param.temp$muQ
    varQ <- param.temp$sigmaQ^2
    df <- param.temp$l

    # get pvalue
    Q.Norm <- (Q - muQ)/sqrt(varQ) * sqrt(2*df) + df # eq 4 from Liu
    pval[,i] <- pchisq(Q.Norm,  df = df, lower.tail=FALSE)

    #pval[,i]<-SKAT:::Get_PValue.Lambda(lambda.temp,Q)$p.value

    param.mat <- rbind(param.mat,c(muQ,varQ,df))
  }

  # get quantiles for each mixture
  pmin <- apply(pval,1,min)
  for(i in 1:n.r){

    muQ <- param.mat[i,1]
    varQ <- param.mat[i,2]
    df <- param.mat[i,3]

    q.org <- qchisq(1-pmin,df=df)
    q.q <- (q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ # set q.org = t*sigmaX + muX and solve for q.q
    pmin.q[,i] <- q.q
  }

  param.m <- list()
  param.m$par.moments <- list()
  param.m$par.moments[[1]] <- param.mat[1,]
  param.m$par.moments[[2]] <- param.mat[n.r,]
  # pmin : min p-value
  # pmin.q : q-values of min p.
  out <- list(pmin=pmin,pval=pval,pmin.q=pmin.q,param.m=param.m)
  return(out)
}

sd_quad <- function(M)
{
  sqrt(2 * sum(M * M))
}

tr <- function(X) sum(diag(X))

perturbation <- function(w, y, X, G, M, grm, num_reps = 10000)
{
  grid_results <- data.frame(w,
                             test.stat = numeric(length(w)),
                             test.pval = numeric(length(w)))
  nullmod <- fitNullModel(data.frame(y, X, row.names = rownames(grm)),
                          outcome = 'y', covars = c('X'), cov.mat = grm,
                          family = "gaussian", start = NULL, drop.zeros = F,
                          return.small = F, verbose = F)
  P0 <- -tcrossprod(nullmod$CXCXI,nullmod$CX); diag(P0) <- 1 + diag(P0);
  P0 <- nullmod$cholSigmaInv %*% P0

  K1 <- tcrossprod(G); K2 <- tcrossprod(M);
  sdQG <- sd_quad(tcrossprod(crossprod(G,P0)));
  sdQM <- sd_quad(tcrossprod(crossprod(M,P0)));
  eta <- sdQG / (sdQG + sdQM);
  K1 <- (1  - eta) * K1; K2 <- eta * K2;

  pert_Lambda <- list()
  bigV <- NULL
  for (i in 1:length(w))
  {
    K <- w[i]*K1 + (1-w[i])*K2
    Q <-  as.numeric(crossprod(nullmod$fit$resid.PY, K) %*% nullmod$fit$resid.PY)
    V <- crossprod(P0, K) %*% P0
    Veig <- eigen(V, symmetric = T)
    wK <- which(Veig$values > 1e-10)
    pert_Lambda[[i]] <- Veig$values[wK]
    bigV <- rbind(bigV,t(Veig$vectors[,wK]))
    f <- CompQuadForm::davies(q=Q, lambda=Veig$values[wK], acc = 1e-9)
    if((f$ifault > 0) | (f$Qq < 1e3*.Machine$double.eps) | (f$Qq > 1))
    {
      ## try saddlepoint
      f$Qq <- survey:::saddle(Q, Veig$values[wK])
      if(is.na(f$Qq)) f$Qq <- 1
    }
    grid_results[i, 'test.stat'] <- Q
    grid_results[i, 'test.pval'] <- f$Qq
  }
  bigSigma <- tcrossprod(bigV)
  rm(bigV)
  bigSigmaEig <- eigen(bigSigma, symmetric = T)
  wK <- which(bigSigmaEig$values < 1e-10)
  bigSigmaEig$values[wK] <- 0
  R <- bigSigmaEig$vectors %*% diag(sqrt(bigSigmaEig$values))
  min_p_stars <- rep(NA,num_reps+1)
  min_p_stars[1] <- as.numeric(min(grid_results$test.pval))
  start <- cumsum(c(1,unlist(lapply(pert_Lambda,length))[-length(pert_Lambda)]))
  end <- cumsum(unlist(lapply(pert_Lambda,length)))
  for (j in 1:num_reps)
  {
    r <- rnorm(ncol(R))
    r_star <- R %*% r

    p_star_min <- Inf
    for (k in 1:length(w))
    {
      r_star_d <- r_star[start[k]:end[k]]
      if (length(r_star_d) == 1)
      {
        Q_star_d <- r_star_d^2 * pert_Lambda[[k]]
      } else {
        Q_star_d <-crossprod(r_star_d, diag(pert_Lambda[[k]])) %*% r_star_d
      }
      f <- CompQuadForm::davies(q=Q_star_d, lambda=pert_Lambda[[k]], acc = 1e-9)
      if((f$ifault > 0) | (f$Qq < 1e3*.Machine$double.eps) | (f$Qq > 1))
      {
        ## try saddlepoint
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