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


.one_replicate <- function(Y, Psi_true, sel, N0, k_true, verbose = FALSE) {
  
  p <- ncol(Y)
  
  # ── FABLEPosteriorSampler / FABLEPosteriorMean (wrapper) ────────
  t_fable   <- proc.time()
  fable_smp  <- FABLEPosteriorSampler(Y, gamma0 = 1, delta0sq = 1, MC = N0)
  fable_mean <- FABLEPosteriorMean(Y, gamma0 = 1, delta0sq = 1)$ FABLEPostMean
  t_fable    <- (proc.time() - t_fable)["elapsed"]
  
  # ── Algorithm 1 ────────────────────────────────────────────────────
  t_alg1 <- proc.time()
  res1   <- algorithm1(Y, k = k_true, N0 = N0,
                       store_samples = TRUE, verbose = verbose)
  t_alg1 <- (proc.time() - t_alg1)["elapsed"]
  
  # ── Error ──────────────────────────────────────────────────────
  err_fable <- norm(fable_mean       - Psi_true, "2") / norm(Psi_true, "2")
  err_alg1  <- norm(res1$Psi_mean   - Psi_true, "2") / norm(Psi_true, "2")
  
  # ── Coverage: (u,v) CI ─────────────────────────────────────
  n_sel <- nrow(sel)
  cov_f <- cov_a <- wid_f <- wid_a <- numeric(n_sel)
  
  for (i in seq_len(n_sel)) {
    u <- sel[i, 1]; v <- sel[i, 2]
    truth <- Psi_true[u, v]
    
    # FABLE CI
    # Armadillo column-major: seq(j, p*k, by=p)
    idx_u <- seq(u, p * k_true, by = p)   # length k_est
    idx_v <- seq(v, p * k_true, by = p)
    
    lambda_u <- fable_smp$CCFABLESamples$LambdaSamples[, idx_u, drop = FALSE]  # MC × k
    lambda_v <- fable_smp$CCFABLESamples$LambdaSamples[, idx_v, drop = FALSE]  # MC × k
    
    # Ψ_{uv}^(m) = λ_u^{(m)T} λ_v^{(m)} + σ²_u^{(m)} · 1(u=v)
    psi_uv_f <- rowSums(lambda_u * lambda_v) +
      fable_smp$CCFABLESamples$SigmaSqSamples[, u] * (u == v)
    
    ci_f     <- quantile(psi_uv_f, c(0.025, 0.975))
    cov_f[i] <- as.numeric(truth >= ci_f[1] & truth <= ci_f[2])
    wid_f[i] <- ci_f[2] - ci_f[1]
    
    # Algorithm 1 CI
    ci_a     <- quantile(res1$Psi_samples[, u, v], c(0.025, 0.975))
    cov_a[i] <- as.numeric(truth >= ci_a[1] & truth <= ci_a[2])
    wid_a[i] <- ci_a[2] - ci_a[1]
  }
  
  list(
    err_fable  = err_fable,  err_alg1  = err_alg1,
    cov_fable  = cov_f,      cov_alg1  = cov_a,
    wid_fable  = wid_f,      wid_alg1  = wid_a,
    time_fable = t_fable,    time_alg1 = t_alg1,
    S_hat_norm = norm(res1$S_hat, "F"),
    V_post_diag = diag(res1$V_post)
  )
}


# =============================================================================
# simulation: Algorithm 1 vs FABLE_code.R (CCFABLE_DirectSampler / FABLEPostmean)
# =============================================================================

#' @param n             # of observation
#' @param p             # of variable
#' @param k_true        # of latent variable
#' @param R             replicate  (default: 100)
#' @param N0            the number of sample
#' @param n_pairs       Coverage  (u,v) 
#' @param seed          reproducibility

run_simulation <- function(
    n       = 500,
    p       = 1000,
    k_true  = 10,
    pi0     = 0.5,   # spike-and-slab 
    v0      = 0.5,   # slab sd  → slab variance = v0^2 = 0.25
    R       = 100,    
    N0      = 1000,
    n_pairs = 100,
    seed    = 1
) {
  
  if (!exists("CCFABLE_DirectSampler"))
    stop("run source('FABLE_code.R')")
  
  cat("=== Replicated simulation (FABLE paper setup) ===\n")
  cat(sprintf("  n=%d | p=%d | k=%d | pi0=%.2f | v0=%.2f | R=%d | N0=%d\n\n",
              n, p, k_true, pi0, v0, R, N0))
  
  # ------------------------------------------------------------------
  # (Λ₀, Σ₀): fixed — 논문 Section 4.1
  # Λ₀[j,l] iid~ π₀·δ₀ + (1-π₀)·N(0, v0²)
  # σ²₀j    iid~ Uniform(0.5, 5)
  # ------------------------------------------------------------------
  set.seed(seed)
  
  spike_mask <- matrix(runif(p * k_true) < pi0, p, k_true)
  slab_draw  <- matrix(rnorm(p * k_true, mean = 0, sd = v0), p, k_true)
  Lambda0    <- ifelse(spike_mask, 0, slab_draw)   # p × k
  
  sigma0_sq  <- runif(p, min = 0.5, max = 5)       
  
  Psi_true   <- tcrossprod(Lambda0) + diag(sigma0_sq, p)
  
  # 논문 식 Pav 확인
  Rav <- k_true * (1 - pi0) * v0^2 / mean(sigma0_sq)
  Pav <- Rav / (1 + Rav)
  cat(sprintf("  Pav (proportion of variance explained) ≈ %.1f%%\n\n",
              100 * Pav))
  
  # ------------------------------------------------------------------
  # (u,v)
  # ------------------------------------------------------------------
  idx <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  sel <- idx[sample(nrow(idx), min(n_pairs, nrow(idx))), , drop = FALSE]
  n_sel <- nrow(sel)

  mat_cov_f  <- matrix(NA, R, n_sel)
  mat_cov_a  <- matrix(NA, R, n_sel)
  mat_wid_f  <- matrix(NA, R, n_sel)
  mat_wid_a  <- matrix(NA, R, n_sel)
  vec_err_f  <- numeric(R)
  vec_err_a  <- numeric(R)
  vec_time_f <- numeric(R)
  vec_time_a <- numeric(R)
  
  # ------------------------------------------------------------------
  # Λ₀, Σ₀ fixed, Y re-generate
  # ------------------------------------------------------------------
  for (r in seq_len(R)) {
    cat(sprintf("  replicate %3d / %d\r", r, R))
    set.seed(seed + r)
    
    F0 <- matrix(rnorm(n * k_true), n, k_true)           # η_i ~ N(0, I_k)
    E  <- matrix(rnorm(n * p) * rep(sqrt(sigma0_sq), each = n), n, p)
    Y  <- scale(F0 %*% t(Lambda0) + E, center = TRUE, scale = FALSE)
    
    #res_r <- tryCatch(
    #  .one_replicate(Y, Psi_true, sel, N0, k_true, verbose = FALSE),
    #  error = function(e) {
    #    message("\n  replicate ", r, " error: ", conditionMessage(e)); NULL
    #  }
    #)
    
    res_r <- tryCatch(
      withCallingHandlers(
        .one_replicate(Y, Psi_true, sel, N0, k_true, verbose = FALSE),
        error = function(e) {
          cat("\n── call stack ──\n")
          calls <- sys.calls()
          for (i in seq_along(calls)) cat(i, ":", deparse(calls[[i]])[1], "\n")
          cat("────────────────\n")
        }
      ),
      error = function(e) {
        message("replicate ", r, " error: ", conditionMessage(e)); NULL
      }
    )
    
    
    if (is.null(res_r)) next
    
    mat_cov_f[r, ]  <- res_r$cov_fable
    mat_cov_a[r, ]  <- res_r$cov_alg1
    mat_wid_f[r, ]  <- res_r$wid_fable
    mat_wid_a[r, ]  <- res_r$wid_alg1
    vec_err_f[r]    <- res_r$err_fable    # ||Ψ̂ - Ψ₀|| / ||Ψ₀||  (spectral)
    vec_err_a[r]    <- res_r$err_alg1
    vec_time_f[r]   <- res_r$time_fable
    vec_time_a[r]   <- res_r$time_alg1
  }
  cat("\n\n")
  

  summarize <- function(x) {
    c(mean = mean(x, na.rm = TRUE),
      lo   = unname(quantile(x, 0.025, na.rm = TRUE)),
      hi   = unname(quantile(x, 0.975, na.rm = TRUE)))
  }
  prt <- function(label, x, fmt = "%6.4f") {
    s <- summarize(x)
    cat(sprintf(paste0("%-36s ", fmt, "  [", fmt, " – ", fmt, "]\n"),
                label, s["mean"], s["lo"], s["hi"]))
  }
  
  cov_r_fable <- rowMeans(mat_cov_f, na.rm = TRUE)
  cov_r_alg1  <- rowMeans(mat_cov_a, na.rm = TRUE)
  wid_r_fable <- rowMeans(mat_wid_f, na.rm = TRUE)
  wid_r_alg1  <- rowMeans(mat_wid_a, na.rm = TRUE)
  
  cat("─────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-36s %6s  [%6s – %6s]\n", "", "Mean", "2.5%", "97.5%"))
  cat("─────────────────────────────────────────────────────────────\n")
  prt("Relative spectral error  FABLE",    vec_err_f)
  prt("Relative spectral error  Algorithm1", vec_err_a)
  cat("─────────────────────────────────────────────────────────────\n")
  prt("Coverage FABLE (rho=b-bar)",  cov_r_fable, fmt = "%6.3f")
  prt("Coverage Algorithm1 (rho=1)", cov_r_alg1,  fmt = "%6.3f")
  cat("─────────────────────────────────────────────────────────────\n")
  prt("Interval width  FABLE",     wid_r_fable)
  prt("Interval width  Algorithm1", wid_r_alg1)
  cat("─────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-36s %6.2f sec\n", "Runtime FABLE (mean)",      mean(vec_time_f)))
  cat(sprintf("%-36s %6.2f sec\n", "Runtime Algorithm1 (mean)", mean(vec_time_a)))
  cat("─────────────────────────────────────────────────────────────\n")
  
  invisible(list(
    Psi_true     = Psi_true,    sel          = sel,
    cov_r_fable  = cov_r_fable, cov_r_alg1   = cov_r_alg1,
    wid_r_fable  = wid_r_fable, wid_r_alg1   = wid_r_alg1,
    err_fable    = vec_err_f,   err_alg1     = vec_err_a,
    mat_cov_fable = mat_cov_f,  mat_cov_alg1 = mat_cov_a
  ))
}


result <- run_simulation(n = 100, p = 90, k_true = 10,
                         pi0 = 0.5, v0 = 0.5,
                         R = 100, N0 = 1000, seed = 1)











# ------------------------------------------------------------------
# True covariance
# n = 500, 1000 & p = 500, 1000, R = 100, MCMC = 1000
# coverage of a randomly chosen 100 by 100 submatrix of Psi0
# ------------------------------------------------------------------
set.seed(1)
n = 1000
p = 500
lambdasd = 0.5
pi0 = 0.5
k = 10

Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) 
Lambda = Lambda * BinMat

Sigma0 = runif(p, 0.5, 5)

M = matrix(rnorm(n*k), nrow = n, ncol = k)
E = matrix(rnorm(n*p), nrow = n, ncol = p)
E = sweep(E, 2, sqrt(Sigma0), "*")

Y = (M %*% t(Lambda)) + E

# FABLEPostMean = FABLEPosteriorMean(Y, gamma0 = 1, delta0sq = 1, maxProp = 0.95)
# FABLESamples = FABLEPosteriorSampler(Y, gamma0 = 1, delta0sq = 1, maxProp = 0.95, MC = 1000)





# If tausq is infinity, then tausq inverse would be ignored.
algorithm1_oracle <- function(Y, Lambda0, Sigma0, k,
                         gamma0 = 1, delta0_sq = 1, rho2 = 1, seed = 1) {
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  
  # tilde F
  LtSinv  <- t(Lambda0 / Sigma0)  # Lamdba0^t %*% Sigma inv
  A0      <- diag(k) + LtSinv %*% Lambda0  #  I_k + Lamdba0 %*% Sigma inv %*% Lamdba0
  A0_inv  <- solve(A0)
  M_mat   <- A0_inv %*% LtSinv  # mean of tilde f_i
  
  F_tilde <- matrix(0, n, k)
  for (i in seq_len(n)) {
    mu_i         <- M_mat %*% Y[i, ]
    F_tilde[i, ] <- MASS::mvrnorm(1, mu = mu_i, Sigma = A0_inv)
  }
  
  # Lamdba and sigmasq
  K    <- solve(crossprod(F_tilde))
  Fty  <- crossprod(F_tilde, Y)
  Mu   <- K %*% Fty  # 
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
  
  list(F_tilde = F_tilde, Lambda_est = Lambda_samp, Sigma_est = diag(sigma2), Psi_est = Psi)
}

res <- algorithm1_oracle(Y = Y, Lambda0 = Lambda, Sigma0 = Sigma0, k = k)
str(res)





















# number of sample added
algorithm1_oracle <- function(Y, Lambda0, Sigma0, k,
                              gamma0 = 1, delta0_sq = 1, rho2 = 1,
                              mcmc = 1, seed = 1) {
  set.seed(seed)
  n <- nrow(Y); p <- ncol(Y)
  
  samples <- vector("list", mcmc)
  
  for (s in seq_len(mcmc)) {
    
    # tilde F
    LtSinv  <- t(Lambda0 / Sigma0)
    A0      <- diag(k) + LtSinv %*% Lambda0
    A0_inv  <- solve(A0)
    M_mat   <- A0_inv %*% LtSinv
    
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
    
    samples[[s]] <- list(F_tilde = F_tilde, Lambda_est = Lambda_samp,
                         Sigma_est = diag(sigma2), Psi_est = Psi)
  }
  
  samples
}

res <- algorithm1_oracle(Y = Y, Lambda0 = Lambda, Sigma0 = Sigma0, k = k, mcmc = 10)
str(res[[1]]); str(res[[2]])


aggregate_samples <- function(res) {
  S <- length(res)
  
  F_tilde_mean  <- Reduce("+", lapply(res, `[[`, "F_tilde"))  / S
  Lambda_mean   <- Reduce("+", lapply(res, `[[`, "Lambda_est")) / S
  Sigma_mean    <- Reduce("+", lapply(res, `[[`, "Sigma_est")) / S
  Psi_mean      <- Reduce("+", lapply(res, `[[`, "Psi_est"))  / S
  
  list(F_tilde  = F_tilde_mean,
       Lambda_est = Lambda_mean,
       Sigma_est  = Sigma_mean,
       Psi_est    = Psi_mean)
}

agg  <- aggregate_samples(res)
str(agg)



















# number of replication added
run_simulation <- function(n, p, k, Lambda_true, Sigma0_true,
                           rep = 100, mcmc = 1000,
                           lambdasd = 1, pi0 = 0.5,
                           gamma0 = 1, delta0_sq = 1, rho2 = 1,
                           seed = 1) {
  set.seed(seed)
  
  results <- vector("list", rep)
  
  for (b in seq_len(rep)) {
    
    # data generation
    M <- matrix(rnorm(n * k), nrow = n, ncol = k)
    E <- matrix(rnorm(n * p), nrow = n, ncol = p)
    E <- sweep(E, 2, sqrt(Sigma0_true), "*")
    Y <- (M %*% t(Lambda_true)) + E
    
    # MCMC samples
    res <- algorithm1_oracle(Y        = Y,
                             Lambda0  = Lambda_true,
                             Sigma0   = Sigma0_true,
                             k        = k,
                             gamma0   = gamma0,
                             delta0_sq = delta0_sq,
                             rho2     = rho2,
                             mcmc = mcmc,
                             seed      = b)        # 
    
    # posterior mean
    agg <- aggregate_samples(res)
    
    results[[b]] <- list(Y   = Y,
                         agg = agg)
  }
  
  results
}

sim <- run_simulation(n        = n,
                      p        = p,
                      k        = k,
                      Lambda_true = Lambda,
                      Sigma0_true = Sigma0,
                      rep    = 10,
                      mcmc = 10)

# Check each rep. est.
str(sim[[1]]$agg$Psi_est); str(sim[[2]]$agg$Psi_est); str(sim[[3]]$agg$Psi_est)
