source("extras/automatic_RockovaGeorge.R")

autoSparse <- function(Y, KUpper = 20) { #only run for lambda0 = 1, 5, 10
  
  # INITIALIZATIONS
  
  n = nrow(Y)
  p = ncol(Y)
  
  K <- KUpper
  G <- p
  
  startB <- matrix(rnorm(G*K),G,K)
  alpha <- 1/G
  lambda1 <- 0.001
  epsilon <- 0.05
  
  # PXEM: Dynamic Posterior Exploration (Approximate M-step)
  
  start <- list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))
  
  result_1 <- FACTOR_ROTATE(Y,1,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
  
  result_5 <- FACTOR_ROTATE(Y,5,lambda1,result_1,K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
  
  result_10 <- FACTOR_ROTATE(Y,10,lambda1,result_5,K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
  
  #in older results, do till lambda0 = 10 and then require for general lambda0
  
  # result_20 <- FACTOR_ROTATE(Y,20,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
  
  # result_30 <- FACTOR_ROTATE(Y,30,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
  
  # final result
  
  #result_lambda <- FACTOR_ROTATE(Y,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
  #result_lambda <- result_30
  
  # Extract Estimates #
  
  LLprime_ROTATE_1 = result_1$B %*% t(result_1$B)
  Sigma_ROTATE_1 = diag(result_1$sigma)
  
  Psi_ROTATE_1 = LLprime_ROTATE_1 + Sigma_ROTATE_1
  
  LLprime_ROTATE_5 = result_5$B %*% t(result_5$B)
  Sigma_ROTATE_5 = diag(result_5$sigma)
  
  Psi_ROTATE_5 = LLprime_ROTATE_5 + Sigma_ROTATE_5
  
  LLprime_ROTATE_10 = result_10$B %*% t(result_10$B)
  Sigma_ROTATE_10 = diag(result_10$sigma)
  
  Psi_ROTATE_10 = LLprime_ROTATE_10 + Sigma_ROTATE_10
  
  allEstimates <- list("EstLambda1" = Psi_ROTATE_1,
                       "EstLambda5" = Psi_ROTATE_5,
                       "EstLambda10" = Psi_ROTATE_10)
  
  return(allEstimates)
  
}

covarianceOutput <- function(ROTATEModel) {
  
  # Extract Estimates #
  
  LLprime_ROTATE = ROTATEModel$B %*% t(ROTATEModel$B)
  Sigma_ROTATE = diag(ROTATEModel$sigma)
  
  Psi_ROTATE = LLprime_ROTATE + Sigma_ROTATE
  
  return(Psi_ROTATE)
  
}

autoSparse_paper <- function(Y, K) {
  
  # INITIALIZATIONS
  
  n = nrow(Y)
  p = ncol(Y)
  
  #K <- 20
  G <- p
  
  startB <- matrix(rnorm(G*K),G,K)
  alpha <- 1/G
  lambda1 <- 0.001
  epsilon <- 0.05
  
  start <- list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))
  allResults <- vector("list", 11)
  allResults[[1]] <- start
  
  # PXEM: Dynamic Posterior Exploration (Approximate M-step)
  
  for(k in 2:11) {
    
    lambda0 <- lambda1 + ((k-2) * 2)
    
    print(paste0("Tempering with lambda0 = ", lambda0))
    
    if(any(is.na(allResults[[k-1]]))) {
      
      allResults[[k]] <- NA
      
    } else {
      
      allResults[[k]] <- FACTOR_ROTATE(Y,lambda0,lambda1,allResults[[k-1]],K,epsilon,alpha,TRUE,TRUE,100,TRUE, plot = FALSE)
      
    }
    
  }
  
  # final result
  
  result_lambda <- allResults[[11]]
  
  # Extract Estimates #
  
  if(any(is.na(result_lambda))) {
    
    return(NA)
    
  } else {
    
    LLprime_ROTATE = result_lambda$B %*% t(result_lambda$B)
    Sigma_ROTATE = diag(result_lambda$sigma)
    
    Psi_ROTATE = LLprime_ROTATE + Sigma_ROTATE
    
    return(Psi_ROTATE)
  
  }
  
}