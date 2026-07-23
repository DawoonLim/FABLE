# usethis::create_from_github("DawoonLim/FABLE")
devtools::load_all()



library(Matrix)
# =============================================================================
# Ancillary functions
# =============================================================================

#' Positive semi-definite projection: (A)_+
#' negative eigenvalues truncate zero
psd_proj <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  eig$vectors %*% diag(pmax(eig$values, 0), nrow(A)) %*% t(eig$vectors)
}

#' matrix square root: B s.t. B %*% t(B) = A
mat_sqrt <- function(A, tol = 1e-10) {
  eig <- eigen(A, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  eig$vectors %*% diag(sqrt(vals), nrow(A))
}

#' Inverse-Gamma sampling: X ~ IG(shape, scale)
rig <- function(n, shape, scale) 1 / rgamma(n, shape = shape, rate = scale)


# =============================================================================
# Algorithm 1
# =============================================================================
#' @param Y       n by p data matrix (centering)
#' @param k       the number of factor (If NULL, use FABLE_code.R: RankEstimator)
#' @param N0      the number of samples
#' @param tau2    loading prior variance (If NULL, Empirical Bayes)
#' @param gamma0  IG prior shape
#' @param delta0_sq IG prior scale
#' @param rho     coverage correction (default=1)
#' @param store_samples store the whole samples of Ψ

algorithm1 <- function(
    Y,
    k           = NULL,
    N0          = 1000,
    tau2        = NULL,
    gamma0      = 1,
    delta0_sq   = 1,
    rho         = 1,
    store_samples = FALSE,
    verbose     = TRUE
) {
  n <- nrow(Y)
  p <- ncol(Y)
  
  # ------------------------------------------------------------------
  # Step 0: SVD, k is known
  # ------------------------------------------------------------------
  if (verbose) cat("Step 0: SVD and k...\n")
  
  svd_Y <- svd(Y) 
  
  if (is.null(k)) {
    kMax <- min(which(cumsum(svd_Y$d) / sum(svd_Y$d) >= 0.95))
    k    <- RankEstimator(Y, svd_Y, kMax)   # FABLE_code.R
    if (verbose) cat("  RankEstimator(JIC) k =", k, "\n")
  }
  
  U <- svd_Y$u[, 1:k, drop = FALSE]   # n × k
  D <- svd_Y$d[1:k]                   # k x k
  V <- svd_Y$v[, 1:k, drop = FALSE]   # p × k
  
  # ------------------------------------------------------------------
  # Step 1a: Σ̂ (FABLE_code.R)
  # σ̂²_j = ||(I - UU^T)y^(j)||² / n
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1a: hat Σ ...\n")
  
  if (k == 1) {
    D_mat <- as.matrix(D)
  } else {
    D_mat <- diag(D, k)
  }
  
  UDVt         <- U %*% D_mat %*% t(V)
  sigma_hat_sq <- colSums((Y - UDVt)^2) / n   # p vector
  
  # ------------------------------------------------------------------
  # Step 1b: Ŝ = V^T Σ̂ V / p
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1b: Ŝ ...\n")
  
  S_hat <- crossprod(V, sweep(V, 1, sigma_hat_sq / p, "*"))   # k x k
  
  # ------------------------------------------------------------------
  # Step 1c: ĈĈ^T = (D²/(np) - Ŝ)₊  Equation (8)
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1c: ĈĈ^T...\n")
  
  CC_hat <- psd_proj(diag(D^2 / (n * p), k) - S_hat)      # k × k, PSD
  C_hat  <- mat_sqrt(CC_hat)                              # k × k, Cholesky factor
  
  # ------------------------------------------------------------------
  # Step 1d: Equation (6)
  # f_i | a_i ~ N_k(M_post %*% a_i,  V_post)
  # M_post = C^T (CC^T + Ŝ)^{-1}  →  C^{-1} as S→0 
  # V_post = I_k - C^T (CC^T + Ŝ)^{-1} C  →  0 as S→0
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1d: conditional posterior parameter...\n")
  
  A        <- Y %*% V / sqrt(p)                         # n × k : A = YV/√p
  CC_plus_S_inv <- solve(CC_hat + S_hat)                # k × k
  M_post   <- t(C_hat) %*% CC_plus_S_inv               # k × k
  V_post   <- psd_proj(diag(k) - M_post %*% C_hat)     # k × k
  V_post_sqrt <- mat_sqrt(V_post)                        # k × k
  
  F_post_mean <- A %*% t(M_post)                        # n × k
  
  # ------------------------------------------------------------------
  # Step 1e: τ² Empirical Bayes (FABLE_code.R)
  # τ̂² = mean(||UDV^T y^(j)||² / n) / (k × σ̂²_j)
  # ------------------------------------------------------------------
  if (is.null(tau2)) {
    YtU  <- sweep(V, 2, D, "*")        # p×k : V diag(D) = (U D)^T Y 
    tau2 <- mean(colSums(t(YtU)^2) / n / (k * sigma_hat_sq))
    tau2 <- max(tau2, 1e-6)
    if (verbose) cat("  Empirical Bayes τ² =", round(tau2, 4), "\n")
  }
  
  # ------------------------------------------------------------------
  # sampling
  # ------------------------------------------------------------------
  if (verbose) cat("sampling (N0 =", N0, ")...\n")
  
  gamma_n <- gamma0 + n
  yy      <- colSums(Y^2)   # p vector
  Psi_mean <- matrix(0, p, p)
  if (store_samples) Psi_samples <- array(0, dim = c(N0, p, p))
  
  for (t in seq_len(N0)) {
    if (verbose && t %% 200 == 0) cat("sample", t, "/", N0, "\r")
    
    # ---------------------------------------------------------------
    # Step 1: f̃_i | a_i 
    # ---------------------------------------------------------------
    F_tilde <- F_post_mean + matrix(rnorm(n * k), n, k) %*% t(V_post_sqrt)
    
    # ---------------------------------------------------------------
    # Step 2: NIG (λ̃_j, σ̃²_j)
    # K = (F̃^T F̃ + I/τ²)^{-1} solve
    # ---------------------------------------------------------------
    FtF_reg <- crossprod(F_tilde) + diag(1 / tau2, k)   # k×k
    K       <- solve(FtF_reg)                              # k×k
    
    Mu         <- K %*% crossprod(F_tilde, Y)             # k×p
    mu_KinvMu  <- colSums(Mu * (FtF_reg %*% Mu))          # p vector
    
    gn_delta2  <- pmax(gamma0 * delta0_sq + yy - mu_KinvMu, 1e-10)
    
    sigma2_tilde <- rig(p, shape = gamma_n / 2, scale = gn_delta2 / 2)
    
    K_sqrt       <- mat_sqrt(K)
    Lambda_tilde <- t(Mu + K_sqrt %*%
                        sweep(matrix(rnorm(k * p), k, p),
                              2, rho * sqrt(sigma2_tilde), "*"))   # p×k
    
    # ---------------------------------------------------------------
    # Step 3: Ψ̃ = Λ̃ Λ̃^T + diag(σ̃²)
    # ---------------------------------------------------------------
    Psi_tilde <- tcrossprod(Lambda_tilde) + diag(sigma2_tilde, p)
    Psi_mean  <- Psi_mean + Psi_tilde / N0
    if (store_samples) Psi_samples[t, , ] <- Psi_tilde
  }
  if (verbose) cat("\n complete.\n")
  
  res <- list(Psi_mean = Psi_mean, F_post_mean = F_post_mean,
              V_post = V_post, S_hat = S_hat, CC_hat = CC_hat,
              sigma_hat_sq = sigma_hat_sq, tau2 = tau2, k = k, rho = rho)
  if (store_samples) res$Psi_samples <- Psi_samples
  return(res)
}


# =============================================================================
# simulation: Algorithm 1 vs FABLE_code.R (CCFABLE_DirectSampler / FABLEPostmean)
# =============================================================================

#' @param n             # of observation
#' @param p             # of variable
#' @param k_true        # of latent variable
#' @param N0            the number of sample
#' @param n_pairs       Coverage  (u,v) 
#' @param seed          reproducibility

run_simulation <- function(
    n = 150, p = 100, k_true = 3,
    N0 = 1000, n_pairs = 50, seed = 42
) {
  set.seed(seed)
  cat("=== simulation: n =", n, ", p =", p, ", k =", k_true, "===\n\n")
  
  Lambda0   <- matrix(rnorm(p * k_true), p, k_true)
  F0        <- matrix(rnorm(n * k_true), n, k_true)
  sigma0_sq <- runif(p, 0.5, 2)
  E         <- matrix(rnorm(n * p) * rep(sqrt(sigma0_sq), each = n), n, p)
  Y         <- scale(F0 %*% t(Lambda0) + E, center = TRUE, scale = FALSE)
  
  Psi_true <- Lambda0 %*% t(Lambda0) + diag(sigma0_sq, p)
  
  # ------------------------------------------------------------------
  # FABLE_code.R : CCFABLE_DirectSampler + FABLEPostmean
  # ------------------------------------------------------------------
  cat("── FABLE...\n")
  if (!exists("CCFABLE_DirectSampler"))
    stop("run source('FABLE_code.R')")
  
  t_fable <- proc.time()
  # coverage computation
  fable_samples <- CCFABLE_DirectSampler(Y, gamma0 = 1, delta0sq = 1, MC = N0)
  # posterior mean
  fable_mean    <- FABLEPostmean(Y, gamma0 = 1, delta0sq = 1)
  t_fable <- (proc.time() - t_fable)["elapsed"]
  
  # ------------------------------------------------------------------
  # Algorithm 1
  # ------------------------------------------------------------------
  cat("── Algorithm 1...\n")
  t_alg1 <- proc.time()
  res1   <- algorithm1(Y, k = k_true, N0 = N0,
                       store_samples = TRUE, verbose = TRUE)
  t_alg1 <- (proc.time() - t_alg1)["elapsed"]
  
  # ------------------------------------------------------------------
  # 2 norm difference
  # ------------------------------------------------------------------
  err_fable <- norm(fable_mean - Psi_true, "2") / norm(Psi_true, "2")
  err_alg1  <- norm(res1$Psi_mean - Psi_true, "2") / norm(Psi_true, "2")
  
  cat("\n── relative spectral error ──\n")
  cat("  FABLE (FABLEPostmean) :", round(err_fable, 4), "\n")
  cat("  Algorithm 1            :", round(err_alg1,  4), "\n")
  
  cat("\n── computation time ──\n")
  cat("  FABLE  :", round(t_fable, 2), "sec.\n")
  cat("  Alg. 1 :", round(t_alg1,  2), "sec.\n")
  
  # ------------------------------------------------------------------
  # Coverage (n_pairs)
  # ------------------------------------------------------------------
  set.seed(seed + 1)
  idx <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  sel <- idx[sample(nrow(idx), min(n_pairs, nrow(idx))), , drop = FALSE]
  
  compute_coverage <- function(samples_array, Psi_true, sel, alpha = 0.05) {
    cov_vec <- width_vec <- numeric(nrow(sel))
    for (i in seq_len(nrow(sel))) {
      u <- sel[i, 1]; v <- sel[i, 2]
      smp <- samples_array[, u, v]
      ci  <- quantile(smp, c(alpha / 2, 1 - alpha / 2))
      cov_vec[i]   <- as.numeric(Psi_true[u, v] >= ci[1] & Psi_true[u, v] <= ci[2])
      width_vec[i] <- ci[2] - ci[1]
    }
    list(coverage = mean(cov_vec), width = mean(width_vec))
  }
  
  cv_fable <- compute_coverage(fable_samples, Psi_true, sel)
  cv_alg1  <- compute_coverage(res1$Psi_samples, Psi_true, sel)
  
  cat("\n── Coverage (95% CI, ", nrow(sel), "pairs) ──\n", sep = "")
  cat("  FABLE (CC-FABLE, ρ=b̄) : coverage =", round(cv_fable$coverage, 3),
      "| width =", round(cv_fable$width, 4), "\n")
  cat("  Algorithm 1 (ρ=1)      : coverage =", round(cv_alg1$coverage, 3),
      "| width =", round(cv_alg1$width, 4), "\n")


  cat("\n── Algorithm 1 ──\n")
  cat("  k:", res1$k, "/ τ²:", round(res1$tau2, 4), "\n")
  cat("  ||Ŝ||_F =", round(norm(res1$S_hat, "F"), 4))
  cat("  Var(f̃_i|a_i) :", round(diag(res1$V_post), 4), "\n")
  
  invisible(list(
    fable_mean    = fable_mean,
    fable_samples = fable_samples,
    res_alg1      = res1,
    Psi_true      = Psi_true,
    err_fable     = err_fable,
    err_alg1      = err_alg1,
    cv_fable      = cv_fable,
    cv_alg1       = cv_alg1
  ))
}



result <- run_simulation(n = 300, p = 270, k_true = 3, N0 = 1000, n_pairs = 50, seed = 1)

