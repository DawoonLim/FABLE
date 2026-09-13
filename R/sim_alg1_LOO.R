# =============================================================
# ================ algorithm1 Leave-One-Out ===================
# =============================================================
# caffeinate -d
Rcpp::sourceCpp("src/updated-FABLE-functions.cpp") # revise; tausq_est = R_PosInf;
library(MASS)

algorithm1_LOO <- function(
    Y,
    k             = NULL,
    MCMC          = 1000,
    tausq         = NULL,
    gamma0        = 1,
    delta0_sq     = 1,
    rho           = 1,
    align         = TRUE,   # V_{-i} -> overall V Procrustes (SVD ratation prevent)
    seed          = 1,
    verbose       = TRUE
) {
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  
  ## overall SVD: order in V
  svd_Y <- svd(Y)
  U <- svd_Y$u[, 1:k, drop = FALSE]
  D <- svd_Y$d[1:k]
  V <- svd_Y$v[, 1:k, drop = FALSE]          # p × k
  
  D_mat        <- if (k == 1) as.matrix(D) else diag(D, k)
  #sigma_hat_sq <- colSums((Y - U %*% D_mat %*% t(V))^2) / n
  #if (is.null(Sigma0)) Sigma0 <- sigma_hat_sq
  Sigma0 <- Sigma0  # True
  
  ## tausq = infinity
  if (is.null(tausq)) {
    YtU   <- sweep(V, 2, D, "*")
    tausq <- max(mean(colSums(t(YtU)^2) / n / (k * sigma_hat_sq)), 1e-6)
    if (verbose) cat("  Empirical Bayes tausq =", round(tausq, 4), "\n")
  }
  
  ## ------------------------------------------------------------------
  ## per-i LOO factor posterior:  f_i | a_i ~ N(Fmean_i, Fvar_i)
  ##   Y_{-i}=U_{-i}D_{-i}V_{-i}^T,  Vtil = V_{-i} R (정렬),  a_i = Vtil^T y_i / sqrt(p)
  ##   Sigma_a = Vtil^T (Y_{-i}^T Y_{-i})_k Vtil /((n-1)p) = R^T D_{-i}^2 R /((n-1)p)
  ##   S_i     = Vtil^T diag(Sigma0/p) Vtil
  ##   CC_i    = (Sigma_a - S_i)_+ ,  C_i C_i^T = CC_i
  ##   M_i     = C_i^T (CC_i+S_i)^{-1} ,  Fvar_i = I - C_i^T (CC_i+S_i)^{-1} C_i
  ## ------------------------------------------------------------------
  if (verbose) cat("Precompute LOO posteriors (", n, " SVDs)...\n", sep = "")
  Fmean <- matrix(0, n, k)
  Fvar  <- vector("list", n) 
  
  for (i in seq_len(n)) {
    if (verbose && i %% 50 == 0) cat("  i =", i, "/", n, "\r")
    sv <- svd(Y[-i, , drop = FALSE], nu = k, nv = k)
    Vi <- sv$v[, 1:k, drop = FALSE]
    Di <- sv$d[1:k]
    
    if (align) {                                  # orthogonal Procrustes: Vi %*% R ~ V
      ss <- svd(crossprod(Vi, V))
      R  <- ss$u %*% t(ss$v)
    } else R <- diag(k)
    Vtil <- Vi %*% R
    
    a_i     <- crossprod(Vtil, Y[i, ]) / sqrt(p)                       # k
    Sigma_a <- t(R) %*% diag(Di^2, k) %*% R / ((n - 1) * p)           # k × k
    S_i     <- crossprod(Vtil, sweep(Vtil, 1, Sigma0 / p, "*"))       # Vtil^T diag(Sigma0/p) Vtil
    
    eig  <- eigen(Sigma_a - S_i, symmetric = TRUE)
    CC_i <- eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors))
    C_i  <- t(chol(CC_i))
    
    W      <- solve(CC_i + S_i, C_i)               # (CC+S)^{-1} C
    Fvar_i <- diag(k) - crossprod(C_i, W)  # I - C^T(CC+S)^{-1}C
    
    Fmean[i, ] <- crossprod(W, a_i)                # M_i a_i = C^T(CC+S)^{-1} a_i
    Fvar[[i]]  <- Fvar_i   
  }
  if (verbose) cat("\n")
  
  ## ------------------------------------------------------------------
  ## Sampling — Step 2,3 (sigma^2, Lambda, Psi)
  ## ------------------------------------------------------------------
  gamma_n     <- gamma0 + n
  yy          <- colSums(Y^2)
  Psi_mean    <- matrix(0, p, p)
  Psi_samples <- array(0, dim = c(MCMC, p, p))
  
  for (t in seq_len(MCMC)) {
    if (verbose && t %% 200 == 0) cat("sample", t, "/", MCMC, "\r")
    
    F_tilde <- matrix(0, n, k)
    for (i in seq_len(n)) F_tilde[i, ] <- MASS::mvrnorm(1, mu = Fmean[i, ], Sigma = Fvar[[i]])
    
    FtF_reg <- crossprod(F_tilde) + diag(1 / tausq, k)
    K       <- solve(FtF_reg)
    
    sigma2_tilde <- numeric(p)
    Lambda_tilde <- matrix(0, p, k)
    for (j in seq_len(p)) {
      mu_j             <- K %*% crossprod(F_tilde, Y[, j])
      delta2_j         <- gamma0 * delta0_sq + yy[j] - sum(mu_j * (FtF_reg %*% mu_j))
      sigma2_tilde[j]  <- 1 / rgamma(1, shape = gamma_n / 2, rate = delta2_j / 2)
      Lambda_tilde[j,] <- MASS::mvrnorm(1, mu = mu_j, Sigma = rho^2 * sigma2_tilde[j] * K)
    }
    
    Psi_tilde        <- tcrossprod(Lambda_tilde) + diag(sigma2_tilde, p)
    Psi_mean         <- Psi_mean + Psi_tilde / MCMC
    Psi_samples[t,,] <- Psi_tilde
  }
  if (verbose) cat("\ncomplete.\n")
  
  list(
    Psi_mean = Psi_mean, Psi_samples = Psi_samples,
    F_post_mean = Fmean, F_post_var = Fvar,
    tausq = tausq, k = k, rho = rho, align = align
  )
}

n = 100 # 500, 1000
p = 150 # 500, 1000
pi0 = 0.5
alpha = 0.05


lambdasd = 0.5
#relevantIndices = c(1:100)

# dir.name = NA #set directory to save
# if(!is.na(dir.name)) {dir.create(dir.name)}

set.seed(1) #set the seed here

relevantIndices = sample(1:p, size = min(p, 100), replace = FALSE) # size = 100
pSub = length(relevantIndices)


k = 10

R = 10

# Sample storage 
#covStor_oracle   <- matrix(0, nrow = R, ncol = pSub * (pSub + 1) / 2) 
#widthStor_oracle <- rep(0, R) 

covStor_alg1     <- matrix(NA, nrow = R, ncol = pSub*(pSub+1)/2)
widthStor_alg1   <- numeric(R)

#covStor = matrix(0, nrow = R, ncol = pSub * (pSub+1)/2)
#widthStor = rep(0, R)




Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) #pi0 = P(zero)
Lambda = Lambda * BinMat

Sigma0 = runif(p, 0.5, 5)

gamma0 = 1
delta0sq = 1  
MC = 1000
Psi0 = Matrix::tcrossprod(Lambda) + diag(Sigma0)  # True cov.

r = 1
for (r in 1:R) {
  
  print(paste0("Replicate: ", r))
  
  set.seed(2001 + r)
  
  M <- matrix(rnorm(n * k), nrow = n, ncol = k)
  E <- matrix(rnorm(n * p), nrow = n, ncol = p)
  E <- sweep(E, 2, sqrt(Sigma0), "*")
  Y <- (M %*% t(Lambda)) + E
  
  svdmod  <- svd(Y)
  U_Y     <- svdmod$u
  V_Y     <- svdmod$v
  svalsY  <- svdmod$d
  
  kEst        <- k
  varInflation <- 1
  
  # ------------------------------------------------------------------
  # algorithm1_LOO 
  # ------------------------------------------------------------------
  res_alg1 <- algorithm1_LOO(
    Y=Y,
    k             = k,
    MCMC          = MC,
    tausq         = Inf,
    gamma0        = 1,
    delta0_sq     = 1,
    rho           = 1,
    align         = TRUE, 
    seed          = 2001 + r,
    verbose       = TRUE  
  )
  
  Psi_sub_arr_alg1 <- aperm(
    res_alg1$Psi_samples[, relevantIndices, relevantIndices],
    c(2, 3, 1)   # MCMC × pSub × pSub  →  pSub × pSub × MCMC
  )
  
  lowPsi_alg1  <- apply(Psi_sub_arr_alg1, c(1, 2), quantile, probs = alpha / 2)
  highPsi_alg1 <- apply(Psi_sub_arr_alg1, c(1, 2), quantile, probs = 1 - alpha / 2)
  
  truePsi0Sub <- Psi0[relevantIndices, relevantIndices]
  trueVec     <- truePsi0Sub[upper.tri(truePsi0Sub, diag = TRUE)]
  lowVec_alg1  <- lowPsi_alg1[ upper.tri(lowPsi_alg1,  diag = TRUE)]
  highVec_alg1 <- highPsi_alg1[upper.tri(highPsi_alg1, diag = TRUE)]
  
  covStor_alg1[r, ]  <- as.numeric((lowVec_alg1 <= trueVec) & (trueVec <= highVec_alg1))
  widthStor_alg1[r]  <- mean(highVec_alg1 - lowVec_alg1)
  
  
  
  ## Save information
  
  #if(!is.na(dir.name)) {
  #  
  #  write.csv(covStor[r,], file = paste0(dir.name, "/", "coverage_rep=", r,  "_n=", n, "_p=", p, "_lambdasd=", lambdasd, "_k=", k, "_pi0=", pi0, ".csv"))
  #  write.csv(widthStor[r], file = paste0(dir.name, "/", "width_rep=", r,  "_n=", n, "_p=", p, "_lambdasd=", lambdasd, "_k=", k, "_pi0=", pi0, ".csv"))
  #  write.csv(mean(covStor[r,]), file = paste0(dir.name, "/", "avg_coverage_rep=", r,  "_n=", n, "_p=", p, "_lambdasd=", lambdasd, "_k=", k, "_pi0=", pi0, ".csv"))
  #  
  #}
  
  # =============================================================
  # Result aggregation
  # =============================================================
  
  cat("=== Oracle ===\n")
  cat("avg coverage:", round(mean(covStor_oracle[1:r, ]), 3),
      " avg width:", round(mean(widthStor_oracle[1:r]), 3), "\n")
  
  
}# end of replication



summarize_simulation <- function(
    covStor_alg1, widthStor_alg1,
    R, n, p) { 
  fmt <- function(x) {
    q <- quantile(x, probs = c(0.025, 0.975))
    sprintf("%.3f [%.3f, %.3f]", mean(x), q[1], q[2])
  }
  
  methods <- list(
    list(label = "Algorithm1_LOO", cov = rowMeans(covStor_alg1[1:R, ]), width = widthStor_alg1[1:R])
  )
  
  cat("=================================================================\n")
  cat(sprintf("  R = %d replications (n=%d, p=%d)\n", R, n, p)) 
  cat("=================================================================\n")
  cat(sprintf("  %-20s %-30s %s\n", "Method", "Coverage", "Width"))
  cat("-----------------------------------------------------------------\n")
  for (m in methods) {
    cat(sprintf("  %-20s %-30s %s\n", m$label, fmt(m$cov), fmt(m$width)))
  }
  cat("=================================================================\n")
}
summarize_simulation(
  covStor_alg1, widthStor_alg1,
  R = R, n = n, p = p)


