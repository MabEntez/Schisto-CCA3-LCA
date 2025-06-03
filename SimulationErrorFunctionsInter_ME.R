## Generate a truncated distribution - from Keith Goldfeld in R-bloggers
rnormt <- function(n, range, mu, s = 1) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
}

## Function to calculate the individual level infection intensity 
modelData <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- c(rep(0,nInd-nInf),rgamma(nInf,shape = sh[ID],rate = rt[ID]))
  
  batch <- sample(1:3, nInd, replace = TRUE)  # Draw random column for each individual

  Intensityday <- matrix(NA, nrow = nInd, ncol = 3)  # Preallocate output matrix

  for (x in 1:3) {
    Intensityday[, x] <- c(
      rnormt(nInd - nInf, range = c(-3, Inf), mu = 0, s = sqrt(intraCCAsd[ID, batch[1:(nInd - nInf)]]^2 + interCCAsd[ID]^2)),
      rnormt(nInf, range = c(0, Inf), 
             mu = 9 / (1 + exp(-multim1_CCA[ID] * (Infectlvl[(nInd - nInf + 1):nInd] - multim2_CCA[ID]))), 
             s = sqrt(intraCCAsd[ID, batch[(nInd - nInf + 1):nInd]]^2 + interCCAsd[ID]^2))
    )
  }
  
  mIntensity1 <- round(Intensityday[,1])
  mIntensity2 <- rowMeans(Intensityday[,1:2])
  mIntensity3 <- rowMeans(Intensityday)
  
  return(c(length(which(mIntensity1>1))/nInd,
           length(which(mIntensity2>1))/nInd,
           length(which(mIntensity3>1))/nInd,
           length(which(mIntensity1>1.5))/nInd,
           length(which(mIntensity2>1.5))/nInd,
           length(which(mIntensity3>1.5))/nInd,
           length(which(mIntensity1>2))/nInd,
           length(which(mIntensity2>2))/nInd,
           length(which(mIntensity3>2))/nInd,
           length(which(mIntensity1>3))/nInd,
           length(which(mIntensity2>3))/nInd,
           length(which(mIntensity3>3))/nInd))
}

## Function to simulate sens and spec
modelSensSpecCCA <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- c(rep(0,nInd-nInf),rgamma(nInf,shape = sh[ID],rate = rt[ID]))
  
  #One batch number per individual per day
  batch_matrix <- matrix(sample(1:3, nInd * 3, replace = TRUE), nrow = nInd, ncol = 3)
  Intensityday <- matrix(NA, nrow = nInd, ncol = 3)
  
  for (x in 1:3) {
    # For non-infected individuals (rows 1 to nInd - nInf)
    noninf_values <- rnormt(nInd - nInf, 
                            range = c(-3, Inf), 
                            mu = 0, 
                            s = sqrt(intraCCAsd[ID, batch_matrix[1:(nInd - nInf), x]]^2 + interCCAsd[ID]^2))
    
    # For infected individuals (rows nInd - nInf + 1 to nInd)
    inf_values <- rnormt(nInf, 
                         range = c(0, Inf), 
                         mu = 9 / (1 + exp(-multim1_CCA[ID] * (Infectlvl[(nInd - nInf + 1):nInd] - multim2_CCA[ID]))), 
                         s = sqrt(intraCCAsd[ID, batch_matrix[(nInd - nInf + 1):nInd, x]]^2 + intraCCAsd[ID]^2))
    
    # Combine the two parts into the column for Intensityday
    Intensityday[, x] <- c(noninf_values, inf_values)
  }
  
  mIntensity1 <- round(Intensityday[,1])
  mIntensity2 <- rowMeans(Intensityday[,1:2])
  mIntensity3 <- rowMeans(Intensityday)
  
  return(c(length(which(mIntensity1[(nInd-nInf+1):nInd]>=1,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<1,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>=1.5,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<1.5,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>=2,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<2,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>=3,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<3,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=1,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<1,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=1.5,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<1.5,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=2,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<2,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=3,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<3,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=1,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<1,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=1.5,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<1.5,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=2,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<2,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=3,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<3,1,0))/(nInd-nInf)))
}


modelSensSpecREC <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- c(rep(0,nInd-nInf),rgamma(nInf,shape = sh[ID],rate = rt[ID]))
  
  # Create a batch matrix with dimensions: one batch number per individual per column
  batch_matrix <- matrix(sample(1:3, nInd * 3, replace = TRUE), nrow = nInd, ncol = 3)
  
  # Preallocate the output matrix
  Intensityday <- matrix(NA, nrow = nInd, ncol = 3)
  
  for (x in 1:3) {
    # For non-infected individuals (rows 1 to nInd - nInf)
    noninf_values <- rnormt(nInd - nInf, 
                            range = c(-3, Inf), 
                            mu = 0, 
                            s = sqrt(intraRECCCAsd[ID, batch_matrix[1:(nInd - nInf), x]]^2 + interRECCCAsd[ID]^2))
    
    # For infected individuals (rows nInd - nInf + 1 to nInd)
    inf_values <- rnormt(nInf, 
                         range = c(0, Inf), 
                         mu = 9 / (1 + exp(-multim1_REC[ID] * (Infectlvl[(nInd - nInf + 1):nInd] - multim2_REC[ID]))), 
                         s = sqrt(intraRECCCAsd[ID, batch_matrix[(nInd - nInf + 1):nInd, x]]^2 + intraRECCCAsd[ID]^2))
    
    # Combine the two parts into the column for Intensityday
    Intensityday[, x] <- c(noninf_values, inf_values)
  }
  mIntensity1 <- round(Intensityday[,1])
  mIntensity2 <- rowMeans(Intensityday[,1:2])
  mIntensity3 <- rowMeans(Intensityday)
  
  return(c(length(which(mIntensity1[(nInd-nInf+1):nInd]>=1,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<1,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>=1.5,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<1.5,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>=2,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<2,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>=3,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<3,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=1,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<1,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=1.5,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<1.5,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=2,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<2,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>=3,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<3,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=1,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<1,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=1.5,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<1.5,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=2,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<2,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>=3,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<3,1,0))/(nInd-nInf)))
}

modelSensSpecKK <- function(nInd, ID) {
  # Number of infected individuals
  nInf <- rbinom(1, nInd, prob = P[ID])
  
  # Simulate infection levels: non-infected get 0, infected get a gamma draw
  Infectlvl <- c(rep(0, nInd - nInf),
                 rgamma(nInf, shape = sh[ID], rate = rt[ID]))
  
  # Generate 6 intensity draws per individual.
  # For non-infected individuals, the probability is K[ID] / (0 + K[ID]) = 1, hence rnbinom returns 0.
  IntensityDraws <- sapply(1:6, function(x) {
    c(
      rnbinom(n = (nInd - nInf), size = K[ID], prob = K[ID] / (0 + K[ID])),
      rnbinom(n = nInf, size = K[ID],
              prob = K[ID] / (Infectlvl[(nInd - nInf + 1):nInd] + K[ID]))
    )
  })
  
  # Calculate daily intensities by averaging pairs of draws:
  day1 <- rowMeans(IntensityDraws[, 1:2])
  day2 <- rowMeans(IntensityDraws[, 1:4])
  day3 <- rowMeans(IntensityDraws[, 1:6])
  
  # For each day, compute sensitivity (infected individuals with value >= 1)
  # and specificity (non-infected individuals with value < 1)
  sens_day1 <- length(which(day1[(nInd - nInf + 1):nInd] >= 1)) / nInf
  spec_day1 <- length(which(day1[1:(nInd - nInf)] < 1)) / (nInd - nInf)
  
  sens_day2 <- length(which(day2[(nInd - nInf + 1):nInd] >= 1)) / nInf
  spec_day2 <- length(which(day2[1:(nInd - nInf)] < 1)) / (nInd - nInf)
  
  sens_day3 <- length(which(day3[(nInd - nInf + 1):nInd] >= 1)) / nInf
  spec_day3 <- length(which(day3[1:(nInd - nInf)] < 1)) / (nInd - nInf)
  
  # Return a vector of sensitivity and specificity for each day:
  # c(sens_day1, spec_day1, sens_day2, spec_day2, sens_day3, spec_day3)
  return(c(sens_day1, spec_day1, sens_day2, spec_day2, sens_day3, spec_day3))
}

## Function to simulate diagnostic prev
modelprevCCA <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- c(rep(0,nInd-nInf),rgamma(nInf,shape = sh[ID],rate = rt[ID]))
  
  #One batch number per individual per day
  batch <- sample(1:3, nInd * 3, replace = TRUE)
  Intensityday <- matrix(NA, nrow = nInd, ncol = 3)
  
  noninf_values <- rnormt(nInd - nInf, 
                          range = c(-3, Inf), 
                          mu = 0, 
                          s = sqrt(intraCCAsd[ID, batch[1:(nInd - nInf)]]^2 + interCCAsd[ID]^2))
  
  # For infected individuals (rows nInd - nInf + 1 to nInd)
  inf_values <- rnormt(nInf, 
                       range = c(0, Inf), 
                       mu = 9 / (1 + exp(-multim1_CCA[ID] * (Infectlvl[(nInd - nInf + 1):nInd] - multim2_CCA[ID]))), 
                       s = sqrt(intraCCAsd[ID, batch[(nInd - nInf + 1):nInd]]^2 + intraCCAsd[ID]^2))
  
  # Combine the two parts into the column for Intensityday
  Intensityday <- c(noninf_values, inf_values)
  
  mIntensity1 <- Intensityday
  
  return(c(length(which(mIntensity1>=1,1,0))/nInd,
           length(which(mIntensity1>=1.5,1,0))/nInd,
           length(which(mIntensity1>=2,1,0))/nInd,
           length(which(mIntensity1>=3,1,0))/nInd))
}


modelprevREC <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- c(rep(0,nInd-nInf),rgamma(nInf,shape = sh[ID],rate = rt[ID]))
  
  #One batch number per individual per day
  batch <- sample(1:3, nInd * 3, replace = TRUE)
  Intensityday <- matrix(NA, nrow = nInd, ncol = 3)
  
  noninf_values <- rnormt(nInd - nInf, 
                          range = c(-3, Inf), 
                          mu = 0, 
                          s = sqrt(intraRECCCAsd[ID, batch[1:(nInd - nInf)]]^2 + interRECCCAsd[ID]^2))
  
  # For infected individuals (rows nInd - nInf + 1 to nInd)
  inf_values <- rnormt(nInf, 
                       range = c(0, Inf), 
                       mu = 9 / (1 + exp(-multim1_REC[ID] * (Infectlvl[(nInd - nInf + 1):nInd] - multim2_REC[ID]))), 
                       s = sqrt(intraRECCCAsd[ID, batch[(nInd - nInf + 1):nInd]]^2 + intraRECCCAsd[ID]^2))
  
  # Combine the two parts into the column for Intensityday
  Intensityday <- c(noninf_values, inf_values)
  
  mIntensity1 <- Intensityday
  
  return(c(length(which(mIntensity1>=1,1,0))/nInd,
           length(which(mIntensity1>=1.5,1,0))/nInd,
           length(which(mIntensity1>=2,1,0))/nInd,
           length(which(mIntensity1>=3,1,0))/nInd))
}

modelprevKK <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- c(rep(0,nInd-nInf),rgamma(nInf,shape = sh[ID],rate = rt[ID]))
  
  Intensityday <- sapply(1, function(x) {
    c(
      rnbinom(n = (nInd - nInf), prob = K[ID]/ (0 + K[ID]), size = K[ID]),  # For non-infected
      rnbinom(n = nInf, prob =K[ID]/ (Infectlvl[(nInd-nInf+1):nInd] + K[ID]), size = K[ID])  # For infected
    )
  })
  mIntensity1 <- Intensityday
  
  return(length(which(mIntensity1>=1,1,0))/nInd)
}
