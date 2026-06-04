library(FABLE)
set.seed(1)
n = 100
p = 500
lambdasd = 0.5
pi0 = 0.85
k = 10

Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) 
Lambda = Lambda * BinMat

Sigma0 = runif(p, 0.5, 5)
Psi0 = (Lambda %*% t(Lambda)) + diag(Sigma0)
normPsi0 = norm(Psi0, type = "2")

R = 50
gamma0Grid = c(0.1, 0.5, 1, 10)
delta0sqGrid = c(0.1, 0.5, 1, 10)
errorMatrixMean = matrix(nrow = length(gamma0Grid), ncol = length(delta0sqGrid))
errorMatrixLower = matrix(nrow = length(gamma0Grid), ncol = length(delta0sqGrid))
errorMatrixUpper = matrix(nrow = length(gamma0Grid), ncol = length(delta0sqGrid))

for(rowInd in 1:length(gamma0Grid)) {
  
  for(colInd in 1:length(delta0sqGrid)) {
    
    gamma0 = gamma0Grid[rowInd]
    delta0sq = delta0sqGrid[colInd]
    
    errorCollect = rep(0, R)
    
    for(nrep in 1:R) {
      
      set.seed(2001 + nrep)
      
      M = matrix(rnorm(n*k), nrow = n, ncol = k)
      E = matrix(rnorm(n*p), nrow = n, ncol = p)
      E = sweep(E, 2, sqrt(Sigma0), "*")
      
      Y = (M %*% t(Lambda)) + E
      
      FABLEPostMean = FABLEPosteriorMean(Y, 
                                         gamma0 = gamma0, 
                                         delta0sq = delta0sq, 
                                         maxProp = 0.99)
      errorValue = norm(FABLEPostMean$FABLEPostMean - Psi0, type = "2") / normPsi0
      
      errorCollect[nrep] = errorValue
      
    }
    
    errorMatrixMean[rowInd, colInd] = mean(errorCollect)
    errorMatrixLower[rowInd, colInd] = quantile(errorCollect, 0.025)
    errorMatrixUpper[rowInd, colInd] = quantile(errorCollect, 0.975)
    
  }
  
}

print(errorMatrixMean)
print(errorMatrixLower)
print(errorMatrixUpper)

