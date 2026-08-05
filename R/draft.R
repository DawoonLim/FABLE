# usethis::create_from_github("DawoonLim/FABLE")
# devtools::load_all()
Rcpp::sourceCpp("src/updated-FABLE-functions.cpp")

library(MASS)
# =============================================================================
# Algorithm 1
# True Sigma, True k, tausq infinity
# =============================================================================
#' @param Y       n by p data matrix (centering)
#' @param k       the number of factor (If NULL, use FABLE_code.R: RankEstimator)
#' @param MCMC      the number of samples
#' @param tausq    loading prior variance (If NULL, Empirical Bayes)
#' @param gamma0  IG prior shape
#' @param delta0_sq IG prior scale
#' @param rho     coverage correction (default=1)
#' @param store_samples store the whole samples of Ψ
algorithm1 <- function(
    Y,
    k             = NULL,
    MCMC          = 1000,
    tausq          = NULL,
    gamma0        = 1,
    delta0_sq     = 1,
    rho           = 1,
    verbose       = TRUE
) {
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
  # Step 1d: f_i | a_i ~ N_k(M_post a_i,  V_post)
  # ------------------------------------------------------------------
  if (verbose) cat("Step 1d: conditional posterior parameters...\n")
  A           <- Y %*% V / sqrt(p)                              # n × k
  M_post      <- t(C_hat) %*% solve(CC_hat + S_hat)             # k × k
  F_post_mean <- A %*% t(M_post)                               # n × k
  
  eig_vp   <- eigen(diag(k) - M_post %*% C_hat, symmetric = TRUE)
  evals_vp <- pmax(eig_vp$values, 0)
  V_post   <- eig_vp$vectors %*% diag(evals_vp, k) %*% t(eig_vp$vectors)   # k × k
  
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
  
  for (t in seq_len(MCMC)) {
    if (verbose && t %% 200 == 0) cat("sample", t, "/", MCMC, "\r")
    
    # -----------------------------------------------------------------
    # Step 1: f̃_i | a_i ~ N_k(M_post a_i, V_post),  i = 1, ..., n
    # -----------------------------------------------------------------
    F_tilde <- matrix(0, n, k)
    for (i in seq_len(n)) {
      F_tilde[i, ] <- MASS::mvrnorm(1, mu = F_post_mean[i, ], Sigma = V_post)
    }
    
    # -----------------------------------------------------------------
    # Step 2: σ̃²_j | ... ~ IG
    # -----------------------------------------------------------------
    FtF_reg      <- crossprod(F_tilde) + diag(1 / tausq, k)   # k × k
    K            <- solve(FtF_reg)                             # k × k
    Mu           <- K %*% crossprod(F_tilde, Y)               # k × p
    mu_KinvMu    <- colSums(Mu * (FtF_reg %*% Mu))            # length p
    # gn_delta2    <- pmax(gamma0 * delta0_sq + yy - mu_KinvMu, 1e-10)
    gn_delta2    <- gamma0 * delta0_sq + yy - mu_KinvMu
    sigma2_tilde <- 1 / rgamma(p, shape = gamma_n / 2, rate = gn_delta2 / 2)
    
    # -----------------------------------------------------------------
    # Step 2 (cont): λ̃_j | σ̃²_j, F̃ ~ N_k(μ_j, ρ² σ̃²_j K),  j = 1, ..., p
    # -----------------------------------------------------------------
    Lambda_tilde <- matrix(0, p, k)
    for (j in seq_len(p)) {
      Lambda_tilde[j, ] <- MASS::mvrnorm(1,
                                         mu    = Mu[, j],
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
    F_post_mean  = F_post_mean,
    V_post       = V_post,
    S_hat        = S_hat,
    CC_hat       = CC_hat,
    sigma_hat_sq = sigma_hat_sq,
    tausq         = tausq,
    k            = k,
    rho          = rho
  )
}
# =============================================================================













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
