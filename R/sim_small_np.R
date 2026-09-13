# =============================================================
# ========= FABLE_mean vs FABLE_samples vs FABLE ==============
# =============================================================
# caffeinate -d
Rcpp::sourceCpp("src/updated-FABLE-functions.cpp") # revise; tausq_est = R_PosInf;
library(MASS)

n = 100   # 500, 1000
p = 150    # 500, 1000


pi0 = 0.5
alpha = 0.05


lambdasd = 0.5
#relevantIndices = c(1:100)

# dir.name = NA #set directory to save
# if(!is.na(dir.name)) {dir.create(dir.name)}

set.seed(1) #set the seed here

relevantIndices = sample(1:p, size = 100, replace = FALSE)  # size = 100, 50
pSub = length(relevantIndices)


k = 10

R = 10

Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) #pi0 = P(zero)
Lambda = Lambda * BinMat

Sigma0 = runif(p, 0.5, 5)

gamma0 = 1
delta0sq = 1  
MC = 1000
Psi0 = Matrix::tcrossprod(Lambda) + diag(Sigma0)  # True cov.



# ---- storage: one coverage matrix + one width vector per method -------------
ncolCov          <- pSub * (pSub + 1) / 2
covStor_mean     <- matrix(0, R, ncolCov)
covStor_draws    <- matrix(0, R, ncolCov)
covStor_fable    <- matrix(0, R, ncolCov)
widthStor_mean   <- rep(0, R)
widthStor_draws  <- rep(0, R)
widthStor_fable  <- rep(0, R)


r = 1
for (r in 1:R) {
  
  print(paste0("Replicate: ", r))
  set.seed(2001 + r)
  
  M <- matrix(rnorm(n * k), n, k)
  E <- sweep(matrix(rnorm(n * p), n, p), 2, sqrt(Sigma0), "*")
  Y <- (M %*% t(Lambda)) + E
  
  svdmod <- svd(Y)
  U_Y <- svdmod$u; V_Y <- svdmod$v; svalsY <- svdmod$d
  kEst <- k
  
  FABLEHypPars      = FABLEHyperParameters(Y, U_Y, V_Y, svalsY, kEst)
  covCorrectEntries = CPPcov_correct_matrix(FABLEHypPars$SigmaSqEstimate,
                                            FABLEHypPars$G)
  varInflation      = mean(covCorrectEntries)^2
  
  
  # ------------------------------------------------------------------
  # FABLE (True k, tausq = infinity) for est. Lambda & Sigma
  # ------------------------------------------------------------------
  CPPSamplingOutput <- CPPFABLESampler(Y, gamma0, delta0sq, MC,
                                       U_Y, V_Y, svalsY, kEst, varInflation)
  
  LambdaDraws <- CPPSamplingOutput$LambdaSamples    # MC x (k*p)
  SigmaDraws  <- CPPSamplingOutput$SigmaSqSamples   # MC x p
  
  # ---------------------------------------------------------------------------
  # (1) Algorithm1_Oracle on FABLE posterior mean
  # ---------------------------------------------------------------------------
  Lambda_FABLE <- matrix(colMeans(LambdaDraws), nrow = p, ncol = k, byrow = TRUE)
  Sigma_FABLE  <- as.numeric(CPPSamplingOutput$SigmaSqEstimatePostMean)
  
  res_mean <- algorithm1_oracle(Y = Y, Lambda0 = Lambda_FABLE, Sigma0 = Sigma_FABLE,
                                k = kEst, gamma0 = gamma0, delta0_sq = delta0sq,
                                rhosq = varInflation, mcmc = MC, seed = 2001 + r)
  
  Psi_sub_arr <- array(
    unlist(lapply(res_mean, function(s) s$Psi_est[relevantIndices, relevantIndices])),
    dim = c(pSub, pSub, MC))

  # quantile truncation
  lowPsi0  = apply(Psi_sub_arr, c(1, 2), quantile, probs = alpha / 2)
  highPsi0 = apply(Psi_sub_arr, c(1, 2), quantile, probs = 1 - alpha / 2)
  
  truePsi0 = Psi0[relevantIndices, relevantIndices]
  truePsi0Vec = truePsi0[upper.tri(truePsi0, diag = TRUE)]
  lowPsi0Vec  = lowPsi0[upper.tri(lowPsi0, diag = TRUE)]
  highPsi0Vec = highPsi0[upper.tri(highPsi0, diag = TRUE)]
  
  covStor_mean[r,]  = as.numeric((lowPsi0Vec <= truePsi0Vec) & (truePsi0Vec <= highPsi0Vec))
  widthStor_mean[r] = mean(highPsi0Vec - lowPsi0Vec)
  
  # ---------------------------------------------------------------------------
  # (2) Algorithm1_Oracle on FABLE draws
  # ---------------------------------------------------------------------------
  Psi_arr_draws <- array(0, dim = c(pSub, pSub, MC))
  Psi_arr_fable <- array(0, dim = c(pSub, pSub, MC))
  
  for (m in seq_len(MC)) {
    Lambda_m <- matrix(LambdaDraws[m, ], nrow = p, ncol = k, byrow = TRUE)
    Sigma_m  <- as.numeric(SigmaDraws[m, ])
    
    # (2) Algorithm1_Oracle, one draw at a time
    res_m <- algorithm1_oracle(Y = Y, Lambda0 = Lambda_m, Sigma0 = Sigma_m,
                               k = kEst, gamma0 = gamma0, delta0_sq = delta0sq,
                               rhosq = varInflation, mcmc = 1, seed = 2001 + r + m)
    Psi_arr_draws[, , m] <- res_m[[1]]$Psi_est[relevantIndices, relevantIndices]
  }
  
  lowPsi0  = apply(Psi_arr_draws, c(1, 2), quantile, probs = alpha / 2)
  highPsi0 = apply(Psi_arr_draws, c(1, 2), quantile, probs = 1 - alpha / 2)
  
  truePsi0 = Psi0[relevantIndices, relevantIndices]
  truePsi0Vec = truePsi0[upper.tri(truePsi0, diag = TRUE)]
  lowPsi0Vec  = lowPsi0[upper.tri(lowPsi0, diag = TRUE)]
  highPsi0Vec = highPsi0[upper.tri(highPsi0, diag = TRUE)]
  
  covStor_draws[r,]  = as.numeric((lowPsi0Vec <= truePsi0Vec) & (truePsi0Vec <= highPsi0Vec))
  widthStor_draws[r] = mean(highPsi0Vec - lowPsi0Vec)
  
  
  # ---------------------------------------------------------------------------
  # (3) plain FABLE
  # ---------------------------------------------------------------------------
  CPPSamplingOutput = CPPFABLESampler(Y, 
                                      gamma0, 
                                      delta0sq, 
                                      MC,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      varInflation)
  
  ##Check coverage of 100*100 submatrix
  CPPPostProcess = CCFABLEPostProcessingSubmatrix(CPPSamplingOutput,
                                                  alpha,
                                                  relevantIndices)
  truePsi0 = Psi0[relevantIndices, relevantIndices]
  lowPsi0 = CPPPostProcess$LowerQuantileMatrix
  highPsi0 = CPPPostProcess$UpperQuantileMatrix
  
  truePsi0Vec = truePsi0[upper.tri(truePsi0, diag = TRUE)]
  lowPsi0Vec = lowPsi0[upper.tri(lowPsi0, diag = TRUE)]
  highPsi0Vec = highPsi0[upper.tri(highPsi0, diag = TRUE)]
  
  covStor_fable[r,] = as.numeric((lowPsi0Vec <= truePsi0Vec) & (truePsi0Vec <= highPsi0Vec))
  widthStor_fable[r] = mean(highPsi0Vec - lowPsi0Vec)
  
}

# =============================================================================
# ================== Result Summary (all three methods) =======================
# =============================================================================
summarize_simulation <- function(
    covStor_mean,  widthStor_mean,
    covStor_draws, widthStor_draws,
    covStor_fable, widthStor_fable,
    R, n, p) {
  fmt <- function(x) {
    q <- quantile(x, probs = c(0.025, 0.975))
    sprintf("%.3f [%.3f, %.3f]", mean(x), q[1], q[2])
  }
  
  methods <- list(
    list(label = "FABLE mean (VarInf)",  cov = rowMeans(covStor_mean[1:R, ]),  width = widthStor_mean[1:R]),
    list(label = "FABLE samples (VarInf)", cov = rowMeans(covStor_draws[1:R, ]), width = widthStor_draws[1:R]),
    list(label = "FABLE(plain)",       cov = rowMeans(covStor_fable[1:R, ]), width = widthStor_fable[1:R])
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
  covStor_mean,  widthStor_mean,
  covStor_draws, widthStor_draws,
  covStor_fable, widthStor_fable,
  R = R, n = n, p = p)

