# =============================================================
# ===== algorithm1 vs algorithm1_oracle vs FABLE ==============
# =============================================================
Rcpp::sourceCpp("src/updated-FABLE-functions.cpp") # revise; tausq_est = R_PosInf;
library(MASS)


# =============================================================================
# Algorithm 1
# True Sigma, True k, tausq infinity
# =============================================================================
#' @param Y       n by p data matrix (centering)
#' @param k       the number of factor 
#' @param MCMC      the number of samples
#' @param tausq    loading prior variance
#' @param gamma0  IG prior shape
#' @param delta0_sq IG prior scale
#' @param rho     coverage correction (default=1)
#' @param store_samples store the whole samples of Ψ
#' @param seed reproducibility
algorithm1 <- function(
    Y,
    k             = NULL,
    MCMC          = 1000,
    tausq          = NULL,
    gamma0        = 1,
    delta0_sq     = 1,
    rho           = 1,
    seed          = 1,
    verbose       = TRUE
) {
  set.seed(seed)
  n <- nrow(Y)
  p <- ncol(Y)
  
  # ------------------------------------------------------------------
  # Step 0: SVD, k is known
  # ------------------------------------------------------------------
  if (verbose) cat("Step 0: SVD and k...\n")
  svd_Y <- svd(Y)
  U <- svd_Y$u[, 1:k, drop = FALSE]   # n × k
  D <- svd_Y$d[1:k]                   # length k
  V <- svd_Y$v[, 1:k, drop = FALSE]   # p × k
  
  # ------------------------------------------------------------------
  # Step 1a: σ̂²_j = ||(I - UU^T) y^(j)||² / n
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1a: hat Sigma...\n")
  D_mat        <- if (k == 1) as.matrix(D) else diag(D, k)
  UDVt         <- U %*% D_mat %*% t(V)
  sigma_hat_sq <- colSums((Y - UDVt)^2) / n       # length p
  
  # ------------------------------------------------------------------
  # Step 1b: Ŝ = V^T diag(σ̂²/p) V
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1b: S_hat...\n")
  #S_hat <- crossprod(V, sweep(V, 1, sigma_hat_sq / p, "*"))   # estimate
  S_hat <- crossprod(V, sweep(V, 1, Sigma0 / p, "*"))   # True
  
  # ------------------------------------------------------------------
  # Step 1c: ĈĈ^T = (D²/(np) - Ŝ)₊,   C̃ s.t. C̃C̃^T = ĈĈ^T
  #   (1) eigenvalue 0  → CC_hat (PSD projection)
  #   (2) chol() C_hat: C_hat = t(chol(CC_hat))  →  C_hat %*% t(C_hat) = CC_hat
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1c: CC_hat and C_hat...\n")
  eig_cc   <- eigen(diag(D^2 / (n * p), k) - S_hat, symmetric = TRUE)
  evals_cc <- pmax(eig_cc$values, 0)
  CC_hat   <- eig_cc$vectors %*% diag(evals_cc, k) %*% t(eig_cc$vectors)
  # C_hat    <- t(chol(CC_hat + diag(1e-10, k)))
  C_hat    <- t(chol(CC_hat))
  # ------------------------------------------------------------------
  # Step 1d: f_i | a_i ~ N_k(M_post a_i,  F_post_var), M_post= C^t (CC^t+Shat)^{-1}
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1d: conditional posterior parameters...\n")
  A           <- Y %*% V / sqrt(p)                                           # n × k
  M_post      <- t(C_hat) %*% solve(CC_hat + S_hat)                           # k × k
  F_post_mean <- A %*% t(M_post)                                               # n × k
  F_post_var      <- diag(k) - t(C_hat) %*% solve(CC_hat + S_hat) %*% C_hat      # k × k
  
  # ------------------------------------------------------------------
  # Step 1e: τ² Empirical Bayes
  # ------------------------------------------------------------------
  if (is.null(tausq)) {
    YtU  <- sweep(V, 2, D, "*")
    tausq <- mean(colSums(t(YtU)^2) / n / (k * sigma_hat_sq))
    tausq <- max(tausq, 1e-6)
    if (verbose) cat("  Empirical Bayes tausq =", round(tausq, 4), "\n")
  }
  
  # ------------------------------------------------------------------
  # Sampling
  # ------------------------------------------------------------------
  if (verbose) cat("Sampling =", MCMC, "...\n")
  gamma_n     <- gamma0 + n
  yy          <- colSums(Y^2)
  Psi_mean    <- matrix(0, p, p)
  Psi_samples <- array(0, dim = c(MCMC, p, p))   # MCMC × p × p
  
  if (verbose) cat("Sampling =", MCMC, "...\n")
  gamma_n     <- gamma0 + n
  yy          <- colSums(Y^2)
  Psi_mean    <- matrix(0, p, p)
  Psi_samples <- array(0, dim = c(MCMC, p, p))   # MCMC × p × p
  
  for (t in seq_len(MCMC)) {
    if (verbose && t %% 200 == 0) cat("sample", t, "/", MCMC, "\r")
    
    # -----------------------------------------------------------------
    # Step 1: f_i | a_i ~ N_k(C^T(CC^T+S)^{-1} a_i, F_post_var),  i = 1,...,n
    # -----------------------------------------------------------------
    F_tilde <- matrix(0, n, k)
    for (i in seq_len(n)) {
      F_tilde[i, ] <- MASS::mvrnorm(1, mu = M_post %*% A[i, ], Sigma = F_post_var)
    }
    
    # -----------------------------------------------------------------
    # Step 2: σ̃²_j | ... ~ IG
    # -----------------------------------------------------------------
    FtF_reg <- crossprod(F_tilde) + diag(1 / tausq, k)   # k × k  (= K^{-1})
    K       <- solve(FtF_reg)                             # k × k
    
    sigma2_tilde <- numeric(p)
    Lambda_tilde <- matrix(0, p, k)
    for (j in seq_len(p)) {
      mu_j            <- K %*% crossprod(F_tilde, Y[, j])            # (F~^T F~ + I/τ²)^{-1} F~^T y^(j)
      delta2_j        <- gamma0 * delta0_sq + yy[j] - sum(mu_j * (FtF_reg %*% mu_j)) # γ0δ0² + ||y^(j)||² - μ_j^T K^{-1} μ_j
      sigma2_tilde[j] <- 1 / rgamma(1, shape = gamma_n / 2, rate = delta2_j / 2)
      Lambda_tilde[j, ] <- MASS::mvrnorm(1, mu = mu_j,
                                         Sigma = rho^2 * sigma2_tilde[j] * K)
    }
    
    # -----------------------------------------------------------------
    # Step 3: Ψ̃ = Λ̃Λ̃^T + diag(σ̃²)
    # -----------------------------------------------------------------
    Psi_tilde        <- tcrossprod(Lambda_tilde) + diag(sigma2_tilde, p)
    Psi_mean         <- Psi_mean + Psi_tilde / MCMC
    Psi_samples[t,,] <- Psi_tilde
  }
  if (verbose) cat("\ncomplete.\n")
  
  list(
    Psi_mean     = Psi_mean,
    Psi_samples  = Psi_samples,
    F_post_mean  = F_post_mean, # posterior mean of tilde f_i
    F_post_var       = F_post_var,
    S_hat        = S_hat,
    CC_hat       = CC_hat,
    sigma_hat_sq = sigma_hat_sq,
    tausq         = tausq,
    k            = k,
    rho          = rho
  )
}



# =============================================================================
# Algorithm 1_Oracle
# True Sigma, True k
# =============================================================================
algorithm1_oracle <- function(Y, Lambda0, Sigma0, k,
                              gamma0 = 1, delta0_sq = 1, rho2 = 1,
                              mcmc = 1, seed = 1) {
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  
  samples <- vector("list", mcmc)
  
  # tilde F
  LtSinv  <- t(Lambda0 / Sigma0)  # Lambda^t Sigma^{-1} = t(Lambda0) %*% solve(diag(Sigma0))
  A0      <- diag(k) + LtSinv %*% Lambda0  #I_k + Lamdba^t Sigma^{-1} Lambda  
  A0_inv  <- solve(A0)
  M_mat   <- A0_inv %*% LtSinv # A^{-1} Lambda^t Sigma^{-1}
  F_post_mean <- Y %*% t(M_mat)
  
  for (s in seq_len(mcmc)) {
    F_tilde <- matrix(0, n, k)
    for (i in seq_len(n)) {
      mu_i         <- M_mat %*% Y[i, ]
      F_tilde[i, ] <- MASS::mvrnorm(1, mu = mu_i, Sigma = A0_inv)
    }
    
    # Lambda and sigmasq
    K    <- solve(crossprod(F_tilde))
    Fty  <- crossprod(F_tilde, Y)
    Mu   <- K %*% Fty
    Kinv <- solve(K)
    
    gn    <- gamma0 + n
    gn_d2 <- numeric(p)
    for (j in seq_len(p)) {
      gn_d2[j] <- gamma0 * delta0_sq + sum(Y[, j]^2) -
        as.numeric(t(Mu[, j]) %*% Kinv %*% Mu[, j])
    }
    
    sigma2 <- numeric(p)
    for (j in seq_len(p)) {
      sigma2[j] <- 1 / rgamma(1, shape = gn / 2, rate = gn_d2[j] / 2)
    }
    
    Lambda_samp <- matrix(0, p, k)
    for (j in seq_len(p)) {
      Lambda_samp[j, ] <- MASS::mvrnorm(1, mu = Mu[, j], Sigma = rho2 * sigma2[j] * K)
    }
    
    # aggregate
    Psi <- tcrossprod(Lambda_samp) + diag(sigma2)
    
    samples[[s]] <- list(F_tilde = F_tilde, F_post_mean = F_post_mean, Lambda_est = Lambda_samp,
                         Sigma_est = diag(sigma2), Psi_est = Psi)
  }
  
  samples
}





# =============================================================================
# Simulation (coverage_FABLE.R)
# =============================================================================
n = 500 # 500, 1000
p = 500 # 500, 1000
pi0 = 0.5
alpha = 0.05


lambdasd = 0.5
#relevantIndices = c(1:100)

# dir.name = NA #set directory to save
# if(!is.na(dir.name)) {dir.create(dir.name)}

set.seed(1) #set the seed here

relevantIndices = sample(1:p, size = 100, replace = FALSE)
pSub = length(relevantIndices)


k = 10

R = 10

# Sample storage 
covStor_oracle   <- matrix(0, nrow = R, ncol = pSub * (pSub + 1) / 2) 
widthStor_oracle <- rep(0, R) 

covStor_alg1     <- matrix(NA, nrow = R, ncol = pSub*(pSub+1)/2)
widthStor_alg1   <- numeric(R)

covStor = matrix(0, nrow = R, ncol = pSub * (pSub+1)/2)
widthStor = rep(0, R)


# posterior samples of F (n x k) 
# mean (n x k)
mdiff_L1_stor  <- numeric(R)
mdiff_L2_stor  <- numeric(R) 
mdiff_max_stor <- numeric(R)
mdiff_frob_stor  <- numeric(R)

gram_L1_stor   <- numeric(R)
gram_frob_stor <- numeric(R)
ggram_L2_stor   <- numeric(R)
gram_max_stor  <- numeric(R)
gram_rel_stor  <- numeric(R)

# variance (k × k) 
vdiff_L1_stor  <- numeric(R)
vdiff_L2_stor  <- numeric(R)   
vdiff_max_stor <- numeric(R)
vdiff_frob_stor <- numeric(R) 
vdiff_rel_stor <- numeric(R)

# trace_alg1_stor <- numeric(R)
#trace_ora_stor  <- numeric(R)
#mse_alg1_stor   <- numeric(R) # mean_i ||mu_alg1_i - M_i||^2  (vs true)
#mse_ora_stor    <- numeric(R) # mean_i ||mu_ora_i  - M_i||^2



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
  # FABLE
  # ------------------------------------------------------------------
  CPPSamplingOutput <- CPPFABLESampler(Y, gamma0, delta0sq, MC,
                                       U_Y, V_Y, svalsY, kEst, varInflation)
  
  t11 <- proc.time()
  CPPPostProcess <- CCFABLEPostProcessingSubmatrix(CPPSamplingOutput, alpha, relevantIndices)
  t12 <- proc.time()
  print(t12 - t11)
  
  truePsi0     <- Psi0[relevantIndices, relevantIndices]
  lowPsi0      <- CPPPostProcess$LowerQuantileMatrix
  highPsi0     <- CPPPostProcess$UpperQuantileMatrix
  truePsi0Vec  <- truePsi0[upper.tri(truePsi0,  diag = TRUE)]
  lowPsi0Vec   <- lowPsi0[ upper.tri(lowPsi0,   diag = TRUE)]
  highPsi0Vec  <- highPsi0[upper.tri(highPsi0,  diag = TRUE)]
  
  covStor[r, ] <- as.numeric((lowPsi0Vec <= truePsi0Vec) & (truePsi0Vec <= highPsi0Vec))
  widthStor[r] <- mean(highPsi0Vec - lowPsi0Vec)
  print(c(mean(covStor[r, ]), widthStor[r]))
  
  # ------------------------------------------------------------------
  # algorithm1_oracle 
  # ------------------------------------------------------------------
  res_oracle <- algorithm1_oracle(Y         = Y,
                                  Lambda0   = Lambda,
                                  Sigma0    = Sigma0,
                                  k         = kEst,
                                  gamma0    = gamma0,
                                  delta0_sq = delta0sq,
                                  rho2      = varInflation,
                                  mcmc      = MC,
                                  seed      = 2001 + r)
  # res_oracle[[s]]$Psi_est : p × p
  
  # submatrix extraction
  Psi_sub_arr <- array(
    unlist(lapply(res_oracle, function(s) s$Psi_est[relevantIndices, relevantIndices])),
    dim = c(pSub, pSub, MC)
  )
  
  # quantile truncation
  lowPsi_oracle  <- apply(Psi_sub_arr, c(1, 2), quantile, probs = alpha / 2)
  highPsi_oracle <- apply(Psi_sub_arr, c(1, 2), quantile, probs = 1 - alpha / 2)
  
  truePsi0Sub  <- Psi0[relevantIndices, relevantIndices]
  trueVec      <- truePsi0Sub[upper.tri(truePsi0Sub, diag = TRUE)]
  lowVec_oracle   <- lowPsi_oracle[upper.tri(lowPsi_oracle,  diag = TRUE)]
  highVec_oracle  <- highPsi_oracle[upper.tri(highPsi_oracle, diag = TRUE)]
  
  covStor_oracle[r, ]  <- as.numeric((lowVec_oracle <= trueVec) & (trueVec <= highVec_oracle))
  widthStor_oracle[r]  <- mean(highVec_oracle - lowVec_oracle)
  
  # ------------------------------------------------------------------
  # algorithm1; rho  = sqrt(varInflation)  
  # ------------------------------------------------------------------
# set.seed(2001 + r)
  res_alg1 <- algorithm1(
    Y         = Y,
    k         = kEst,
    MCMC      = MC,
    tausq     = Inf,
    gamma0    = gamma0,
    delta0_sq = delta0sq,
    rho       = 1,
    seed = 2001 + r,
    verbose   = FALSE
  )
  # res_alg1$Psi_samples : MCMC × p × p
  
  # ------------------------------------------------------------------
  # F posterior comparison
  # ------------------------------------------------------------------
  # oracle result
  LtSinv  <- t(Lambda0 / Sigma0)  # Lambda^t Sigma^{-1} = t(Lambda0) %*% solve(diag(Sigma0))
  A0      <- diag(k) + LtSinv %*% Lambda0  #I_k + Lamdba^t Sigma^{-1} Lambda  
  A0_inv  <- solve(A0)
  M_mat   <- A0_inv %*% LtSinv # A^{-1} Lambda^t Sigma^{-1}
  
  F_mean_ora <- Y %*% t(M_mat)
  F_var_ora <- A0_inv
  
  
  # algorithm1 result
  F_mean_alg1  <- res_alg1$F_post_mean   # n × k
  F_var_alg1 <- res_alg1$F_post_var        # k × k
  
  # posterior mean diff
  Mdiff <- F_mean_alg1 - F_mean_ora 
  #trace_mean_alg1_stor[r] <- sum(diag(F_mean_alg1))
  #trace_mean_ora_stor[r]  <- sum(diag(F_mean_ora))
  mdiff_L1_stor[r]  <- norm(Mdiff, "1")
  mdiff_L2_stor[r]  <- norm(Mdiff, "2")
  mdiff_max_stor[r] <- norm(Mdiff, "M")
  mdiff_frob_stor[r] <- norm(Mdiff, "F") 
  #mse_alg1_stor[r]  <- mean(rowSums((F_mean_alg1 - M)^2))
  #mse_ora_stor[r]   <- mean(rowSums((F_mean_ora  - M)^2))
  
  # Gram matrix diff (F %*% Ft); (n × n)
  PP_gram <- tcrossprod(F_mean_alg1)  # P: posterior mean of algorithm1 
  QQ_gram <- tcrossprod(F_mean_ora)   # Q: posterior mean of oracle
  GG_diff <- PP_gram - QQ_gram
  
  gram_L1_stor[r]   <- norm(GG_diff, "1")
  ggram_L2_stor[r]   <- norm(GG_diff, "2")
  gram_max_stor[r]  <- norm(GG_diff, "M")
  gram_frob_stor[r] <- norm(GG_diff, "F")
  gram_rel_stor[r] <- norm(GG_diff, "F") / norm(QQ_gram, "F") # GG^t / QQ^t
  
  # posterior variance diff
  Vdiff <- F_var_alg1 - F_var_ora  
  #trace_var_alg1_stor[r] <- sum(diag(F_var_alg1))
  #trace_var_ora_stor[r]  <- sum(diag(F_var_ora))
  vdiff_L1_stor[r]   <- norm(Vdiff, "1")        
  vdiff_L2_stor[r]   <- norm(Vdiff, "2")            
  vdiff_max_stor[r]  <- norm(Vdiff, "M")          
  vdiff_frob_stor[r] <- norm(Vdiff, "F")  
  vdiff_rel_stor[r]  <- norm(Vdiff, "F") / norm(F_var_ora, "F") 
  
  cat(sprintf(paste0(
    "[Rep %d]\n",
    "  F_mean diff | L1: %8.2f  L2: %8.2f  Max: %8.2f  Frob: %8.2f\n",
    "  F_var  diff | L1: %8.2f  L2: %8.2f  Max: %8.2f  Frob: %8.2f\n"
  ),
  r,
  mdiff_L1_stor[r], mdiff_L2_stor[r], mdiff_max_stor[r], mdiff_frob_stor[r],
  vdiff_L1_stor[r], vdiff_L2_stor[r], vdiff_max_stor[r], vdiff_frob_stor[r]
  ))
  
  
  # ------------------------------------------------------------------
  # Coverage / width — submatrix quantile
  #   Psi_samples: MCMC × p × p -> relevantIndices -> aperm pSub × pSub × MCMC apply
  # ------------------------------------------------------------------
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
  
  print(paste0("Rep ", r, " | FABLE cov: ",  round(mean(covStor[r, ]), 3),
               " width: ", round(widthStor[r], 3),
               " | Oracle cov: ", round(mean(covStor_oracle[r, ]), 3),
               " width: ", round(widthStor_oracle[r], 3),
               " | Alg1 cov: ",   round(mean(covStor_alg1[r, ]), 3),
               " width: ",        round(widthStor_alg1[r], 3)))
  
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
  cat("=== FABLE ===\n")
  cat("avg coverage:", round(mean(covStor[1:r, ]), 3),
      " avg width:", round(mean(widthStor[1:r]), 3), "\n")
  
  cat("=== Oracle ===\n")
  cat("avg coverage:", round(mean(covStor_oracle[1:r, ]), 3),
      " avg width:", round(mean(widthStor_oracle[1:r]), 3), "\n")
  
  cat("=== Algorithm1 ===\n")
  cat("avg coverage:", round(mean(covStor_alg1[1:r, ]), 3),
      " avg width:",   round(mean(widthStor_alg1[1:r]), 3), "\n")

}# end of replication

# =============================================================
# ================== Result Summary ===========================
# =============================================================
summarize_simulation <- function(covStor, widthStor,
                                 covStor_oracle, widthStor_oracle,
                                 covStor_alg1,   widthStor_alg1,
                                 R, n, p) { 
  fmt <- function(x) {
    q <- quantile(x, probs = c(0.025, 0.975))
    sprintf("%.4f [%.4f, %.4f]", mean(x), q[1], q[2])
  }
  
  methods <- list(
    list(label = "FABLE",             cov = rowMeans(covStor[1:R, ]),        width = widthStor[1:R]),
    list(label = "Algorithm1_Oracle", cov = rowMeans(covStor_oracle[1:R, ]), width = widthStor_oracle[1:R]),
    list(label = "Algorithm1",        cov = rowMeans(covStor_alg1[1:R, ]),   width = widthStor_alg1[1:R])
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
summarize_simulation(covStor, widthStor,
                     covStor_oracle, widthStor_oracle,
                     covStor_alg1,   widthStor_alg1,
                     R = R, n = n, p = p)


cat("\n========================================================================\n")
cat("  Algorithm1 vs Oracle — posterior comparison (mean over R reps)\n")
cat("========================================================================\n")
cat(sprintf("  %-22s  %8s  %8s  %8s  %8s  %8s\n",
            "Metric", "L1", "L2", "Max", "Frob", "Rel(Frob)"))
cat("------------------------------------------------------------------------\n")
cat(sprintf("  %-22s  %8.4f  %8.4f  %8.4f  %8.4f  %8s\n",
            "F_mean diff (n×k)",
            mean(mdiff_L1_stor[1:R]), mean(mdiff_L2_stor[1:R]),
            mean(mdiff_max_stor[1:R]), mean(mdiff_frob_stor[1:R]), "—"))
cat(sprintf("  %-22s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
            "Gram diff (n×n)",
            mean(gram_L1_stor[1:R]),  mean(ggram_L2_stor[1:R]),
            mean(gram_max_stor[1:R]), mean(gram_frob_stor[1:R]),
            mean(gram_rel_stor[1:R])))
cat(sprintf("  %-22s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
            "F_var diff (k×k)",
            mean(vdiff_L1_stor[1:R]), mean(vdiff_L2_stor[1:R]),
            mean(vdiff_max_stor[1:R]), mean(vdiff_frob_stor[1:R]),
            mean(vdiff_rel_stor[1:R])))
cat("========================================================================\n")
