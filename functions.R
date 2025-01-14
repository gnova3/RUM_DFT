############# Functions ############
# Function to compute the probability of deciding to choose at step t

omega <- function(VA, VB, VC, VD, att1, att2, alpha, mu_1, mu_2, w_1, w_2, delta, b1, b2) {
  # Compute the probability of deciding to choose at step t
  #
  # VA, VB, VC, VD: Systematic utility of each alternative
  # att1, att2: Attribute values (length = 4 for each)
  # alpha: Memory parameter (0 <= alpha <= 1)
  # mu_1, mu_2: Scales of utilities
  # w_1, w_2: Attention weights
  # delta: threshold
  # b1, b2: Attribute parameters
  #
  # return PE float, Probability of deciding to choose at step t

  # Expected maximum of current utilities (Eq 9)
  EMU_t <- log(sum(exp(mu_1 * c(VA, VB, VC, VD)))) / mu_1
  
  # Future utilities at step t+1 for attributes 1 and 2 (Eq 8)
  u1 <- update(VA, VB, VC, VD, alpha, att1, b1)
  u2 <- update(VA, VB, VC, VD, alpha, att2, b2)
  
  # Expected maximum of future utilities (Eq 10)
  EMU_1 <- log(sum(exp(mu_2 * u1))) / mu_2
  EMU_2 <- log(sum(exp(mu_2 * u2))) / mu_2
  
  # Probability of attending to attributes (Eq 6)
  P_attr <- P_k(w_1, w_2)
  P1 <- P_attr[1]
  P2 <- P_attr[2]
  
  # Weighted expected maximum of future utilities (Eq 4)
  EMU_t_1 <- P1 * EMU_1 + P2 * EMU_2
  
  # Probability of choosing at step t (Eq 3)
  PE <- 1 / (1 + exp(-1 * (EMU_t_1 - EMU_t + delta))) -   1 / (1 + exp(-1 * (EMU_t_1 - EMU_t - delta)))
  
  return(PE)
}

# Function to compute the probability of choosing each alternative
P_i <- function(VA, VB, VC, VD) {
  # Compute the probability of choosing each alternative
  #
  # VA, VB, VC, VD: Systematic utility of each alternative
  #
  # return Probabilities for each alternative
  
  P <- exp(c(VA, VB, VC, VD)) / sum(exp(c(VA, VB, VC, VD)))
  
  return(P)
}

# Function to compute the probability of attending to each attribute
P_k <- function(w_1, w_2) {
  # Compute the probability of attending to each attribute
  #
  # w_1, w_2: Attention weights
  #
  # return Probabilities for each attribute

  PC <- exp(w_1) / (exp(w_1) + exp(w_2))
  PT <- 1 - PC
  Pk  <- c(PC, PT)
  
  return(Pk)
}

# Function to update utilities with attribute contributions
update <- function(VA, VB, VC, VD, alpha, att, b) {
  # Update utilities based on attributes and memory parameter
  #
  # VA, VB, VC, VD: Systematic utility of each alternative
  # alpha: Memory parameter
  # att: Attribute values 
  # b: Attribute parameter
  #
  # return Updated utilities

  V <- alpha * c(VA, VB, VC, VD) + (1 - alpha) * b * att
  
  return(V)
}

# Function to compute all probabilities for the log-likelihood computation
kernel <- function(VA, VB, VC, VD, mu_1, mu_2, alpha, w_1, w_2, att1, att2, delta, b1, b2) {
  # Compute all probabilities needed for log-likelihood at step t
  #
  # VA, VB, VC, VD: Systematic utility of each alternative
  # mu_1, mu_2: Scales of utilities
  # alpha: Memory parameter
  # w_1, w_2: Attention weights
  # att1, att2: Attribute values 
  # delta: Threshold
  # b1, b2: Attribute parameters
  #
  # return Probabilities [PE, PS, P1, P2, P3, P4, PC, PT]

  # Probability of deciding to choose and to continue searching
  PE <- omega(VA, VB, VC, VD, att1, att2, alpha, mu_1, mu_2, w_1, w_2, delta, b1, b2)
  PS <- 1 - PE
  
  # Probability of choosing each alternative
  Pi <- P_i(VA, VB, VC, VD)
  
  # Probability of attending each attribute
  P_attr <- P_k(w_1, w_2)
  PC <- P_attr[1]
  PT <- P_attr[2]
  
  P <- c(PE, PS, Pi, PC, PT)
  return(P)
}


###### RUM-DFT-ISP #############
LL_total <- function(betas) {
  # Log-likelihood calculation for the RUM-DFT-ISP model
  #
  # betas: Vector of model parameters.
  #   - `bt`: Beta for time
  #   - `bc`: Beta for cost
  #   - `Ft`: Memory parameter
  #   - `U_A`, `U_B`, `U_D`: Preconceived utilities for alternatives A, B, and D
  #   - `Wc`: Attention weight for the cost attribute
  #   - `d`: Threshold
  #   - `mu`: scales  
  # return sum of log-likelihood.

  
  # Extract individual parameters from betas vector
  bt  = betas[1]    # beta time
  bc  = betas[2]    # beta cost
  Ft  = betas[3]    # Memory parameter
  U_A = betas[4]    # Preconceived utility for alternative A
  U_B = betas[5]    # Preconceived utility for alternative C
  U_D = betas[6]    # Preconceived utility for alternative D
  Wc = betas[7]     # Attention weight for cost attribute
  d  = betas[8]     # threshold
  ue = betas[9]#1
  us = betas[10]#1
  
  # Extract data from global variables
  MAT     <- H
  Tmax_n  <- choicedata_RUM_DFT$tiempoxid
  TMAX    <- max(Tmax_n)
  choice  <- cbind(choicedata_RUM_DFT$c1, choicedata_RUM_DFT$c2, choicedata_RUM_DFT$c3, choicedata_RUM_DFT$c4)
  
  # Searches verifications
  count_cost <- rowSums(H == 5, na.rm = TRUE)
  count_time <- rowSums(H == 6, na.rm = TRUE)
  
  # Initialize vectors and matrices to store results
  N           <- length(choicedata_RUM_DFT[,1])
  Nc          <- vector(length = N) * 0   # Number of cost searches
  Nt          <- vector(length = N) * 0   # Number of time searches
  log_Pj      <- vector(length = N) * 0   # Log-likelihood of choosing an alternative
  log_Phi     <- vector(length = N) * 0   # Log-likelihood of attending an attribute
  log_Omega_E <- vector(length = N) * 0   # Log-likelihood of deciding to choose or search
  log_Omega_S <- vector(length = N) * 0   # Log-likelihood of deciding to continue searching
  
  # Loop through individuals
  for (n in 1:N) {
    VA   <- 0
    VB  <- 0
    VC  <- 0
    VD   <- 0
    
    # Extract attributes for each individual
    TDVn  = cbind(choicedata_RUM_DFT$TVA[n], choicedata_RUM_DFT$TVB[n], choicedata_RUM_DFT$TVC[n], choicedata_RUM_DFT$TVD[n])                               
    TDVn  = apply(TDVn, 2, as.numeric) 
    COSTn = cbind(choicedata_RUM_DFT$COSTA[n], choicedata_RUM_DFT$COSTB[n], choicedata_RUM_DFT$COSTC[n], choicedata_RUM_DFT$COSTD[n])
    COSTn = apply(COSTn, 2, as.numeric) 
    Tn    = choicedata_RUM_DFT$tiempoxid[n]
    
    # Loop through steps
    t = 1
    while (t <= Tn) {

      # Initialize the utility of each alternative with preconceived utilities
      if (t == 1) {  
        VA <- U_A
        VB <- U_B 
        VC <- 0 
        VD <- U_D  
      }
      
      # Call the kernel function to get probabilities a
      K = kernel(VA, VB, VC, VD, ue, us, Ft, Wc, 0, COSTn, TDVn, d * t^2, bc, bt)

      # Check if the individual chose an alternative
      if (MAT[n, t] == 0) {
        log_Omega_E[n] = log(K[1])
        log_Pj[n] = log(sum(K[3:6] * choice[n,]))
        t = Tn + 1

      } 
      # Check if the individual decided to continue searching     
      else {
        log_Omega_S[n] = log_Omega_S[n] + log(K[2])
        
        # Update utilities based on the chosen attribute
        if (MAT[n, t] == 5) { 
          u_next = update(VA, VB, VC, VD, Ft, COSTn, bc)
          Nc[n] = Nc[n] + 1
        } else {
          u_next = update(VA, VB, VC, VD, Ft, TDVn, bt)
          Nt[n] = Nt[n] + 1

        }
        
        VA <- u_next[1]
        VB <- u_next[2]
        VC <- u_next[3]
        VD <- u_next[4]
      }
      
      # Update the step
      t = t + 1
    }
    
    # Calculate log-likelihood of attending an attribute
    log_Phi[n] = Nc[n] * Wc - (Tn - 1) * log(1 + exp(Wc))
  }
  
  # Calculate log-likelihood considering information
  LL = log_Pj + log_Omega_E + log_Phi + log_Omega_S  
  
  # Return the negative sum of log-likelihood
  return(sum(-LL))       
}

###### RUM-DFT-SC Swissmetro ##### 
LL_Swiss <- function(beta){
  
  Bt        <- beta[1] #TDV
  Bc        <- beta[2] #Cost
  Oitrain   <- beta[3]
  Oism      <- 0
  Oicar     <- beta[4]
  desct     <- abs(beta[5])
  alpha     <- abs(beta[6])
  beta_t    <- beta[7]
  
  u     <- 1
  ue    <- 1
  us    <- 1
  
  #################  ################# Configuracion matrices
  ALT <- 2
  t <- 0:(TMAX - 1)
  H <- sum(ALT ^ t)
  T2 <- TMAX - 2
  t2 <- 0:T2
  H2 <- sum(ALT ^ t2)
  #################  ################# Matrices
  
  Vst2     <- matrix(nrow = N, ncol = H2)
  Vsc2     <- matrix(nrow = N, ncol = H2)
  Vs2      <- matrix(nrow = N, ncol = H2)
  Ve2      <- matrix(nrow = N, ncol = H2)
  CAMINO2  <- matrix(nrow = N, ncol = H)
  
  Vtrain2   <- matrix(nrow = N, ncol = H)
  Vsm2      <- matrix(nrow = N, ncol = H)
  Vcar2     <- matrix(nrow = N, ncol = H)
  
  PS2        <- matrix(nrow = N, ncol = H)
  PE2        <- matrix(nrow = N, ncol = H)
  
  Prtrainm2 <- matrix(nrow = N, ncol = H)
  Prsmm2    <- matrix(nrow = N, ncol = H)
  Prcarm2   <- matrix(nrow = N, ncol = H)
  
  for (h in 1:H) {
    if (h == 1) {
      Vtrain2[, h]<- Oitrain
      Vsm2[, h]   <- Oism
      Vcar2[, h]  <- Oicar
    }
    else{
      aux <- h %/% 2
      if (h %% 2 == 0) {
        Vtrain2[, h]   <- alpha * Vtrain2[, aux] + (1 - alpha) * Bt * train2[, 2]
        Vsm2[, h]      <- alpha * Vsm2[, aux]    + (1 - alpha) * Bt * sm2[, 2]
        Vcar2[, h]     <- alpha * Vcar2[, aux]   + (1 - alpha) * Bt * car2[, 2]
      }
      else{
        Vtrain2[, h]   <- alpha * Vtrain2[, aux] + (1 - alpha) * Bc * train2[, 1]
        Vsm2[, h]      <- alpha * Vsm2[, aux]    + (1 - alpha) * Bc * sm2[, 1]
        Vcar2[, h]     <- alpha * Vcar2[, aux]   + (1 - alpha) * Bc * car2[, 1]
      }
    }
  }
  #print(Vtrain2[1,])
  for (h in 1:H2) {
    auxt <- 2 * h
    auxc <- 2*h + 1
    Vst2[, h] <-(1 / ue) * log(exp(ue * Vtrain2[, auxt]) + exp(ue * Vsm2[, auxt]) + exp(ue *Vcar2[, auxt]) )
    Vsc2[, h] <-(1 / ue) * log(exp(ue * Vtrain2[, auxc]) + exp(ue * Vsm2[, auxc]) + exp(ue *Vcar2[, auxc]) )
    
  }
  for (h in 1:H2) {
    Vs2[,h]<- (1 / us) * log(exp(us * Vst2[,h]) + exp(us * Vsc2[,h]))
    Ve2[,h]<- (1 / ue) * log(exp(ue * Vtrain2[,h]) + exp(ue * Vsm2[,h]) + exp(ue * Vcar2[,h]) )
  }
  for (h in 1:H2) {
    t <- log2(h) %/% 1 
    desc <- desct*((t)^2)
    PS2[,h]<- (exp(u * Vs2[,h]) / (exp(u*Vs2[,h]) + exp(u * (Ve2[,h] + desc))) ) + (exp(u * Ve2[,h]) / (exp(u * (Vs2[,h]  + desc )) + exp(u * Ve2[,h])))
    
  }
  
  PS2[PS2>1]<-1
  
  PE2 <- 1 - PS2
  PE2[is.na(PE2)] <- 1
  PS2[is.na(PS2)] <- 0
  
  K_t<-beta_t
  K_c<-0
  denom_K<- exp(us*K_t)+exp(us*K_c)
  PT2  <- exp(us*K_t)/denom_K
  PC2  <- exp(us*K_c)/denom_K
  
  denom2 <- exp(ue*Vtrain2) + exp(ue*Vsm2) + exp(ue*Vcar2)
  P12 <- exp(ue*Vtrain2)  / denom2
  P22 <- exp(ue*Vsm2) / denom2
  P32 <- exp(ue*Vcar2) / denom2
  
  for (h in 1:H) {
    aux <- h %/% 2
    h2 <- h
    camino <- 1
    while (aux > 0) {
      if (h2 %% 2 == 0) {
        camino <- camino * PS2[, aux] * PT2
        h2 <- aux
        aux <- aux %/% 2
      }
      else{
        camino <- camino * PS2[, aux] * PC2
        h2 <- aux
        aux <- aux %/% 2
      }
    }
    CAMINO2[,h]  <- camino
  }
  CAMINO2<-CAMINO2*PE2
  Prtrainm2  <-  P12 * CAMINO2
  Prsmm2     <-  P22 * CAMINO2
  Prcarm2    <-  P32 * CAMINO2
  
  Prtrain2  <- rowSums(Prtrainm2)
  Prsm2     <- rowSums(Prsmm2)
  Prcar2    <- rowSums(Prcarm2)
  
  Prtrain2[is.na( Prtrain2)] <- 0
  Prsm2[is.na( Prsm2)] <- 0
  Prcar2[is.na( Prcar2)] <- 0
  
  lli2  <- log((Prtrain2 * choice[, 1]) + (Prsm2 * choice[, 2]) + (Prcar2 * choice[, 3]))
  
  
  return(-sum(lli2))  
}
LL_Swiss_fixed <- function(beta){
  
  Bt        <- beta[1] #TDV
  Bc        <- beta[2] #Cost
  Oitrain   <- beta[3]
  Oism      <- 0
  Oicar     <- beta[4]
  desct     <- abs(beta[5])
  beta_t    <- beta[6]
  
  u     <- 1
  ue    <- 1
  us    <- 1
  
  #################  ################# Configuracion matrices
  ALT <- 2
  t <- 0:(TMAX - 1)
  H <- sum(ALT ^ t)
  T2 <- TMAX - 2
  t2 <- 0:T2
  H2 <- sum(ALT ^ t2)
  #################  ################# Matrices
  
  Vst2     <- matrix(nrow = N, ncol = H2)
  Vsc2     <- matrix(nrow = N, ncol = H2)
  Vs2      <- matrix(nrow = N, ncol = H2)
  Ve2      <- matrix(nrow = N, ncol = H2)
  CAMINO2  <- matrix(nrow = N, ncol = H)
  
  Vtrain2   <- matrix(nrow = N, ncol = H)
  Vsm2      <- matrix(nrow = N, ncol = H)
  Vcar2     <- matrix(nrow = N, ncol = H)
  
  PS2        <- matrix(nrow = N, ncol = H)
  PE2        <- matrix(nrow = N, ncol = H)
  
  Prtrainm2 <- matrix(nrow = N, ncol = H)
  Prsmm2    <- matrix(nrow = N, ncol = H)
  Prcarm2   <- matrix(nrow = N, ncol = H)
  
  for (h in 1:H) {
    if (h == 1) {
      Vtrain2[, h]<- Oitrain
      Vsm2[, h]   <- Oism
      Vcar2[, h]  <- Oicar
    }
    else{
      aux <- h %/% 2
      if (h %% 2 == 0) {
        Vtrain2[, h]   <- alpha * Vtrain2[, aux] + (1 - alpha) * Bt * train2[, 2]
        Vsm2[, h]      <- alpha * Vsm2[, aux]    + (1 - alpha) * Bt * sm2[, 2]
        Vcar2[, h]     <- alpha * Vcar2[, aux]   + (1 - alpha) * Bt * car2[, 2]
      }
      else{
        Vtrain2[, h]   <- alpha * Vtrain2[, aux] + (1 - alpha) * Bc * train2[, 1]
        Vsm2[, h]      <- alpha * Vsm2[, aux]    + (1 - alpha) * Bc * sm2[, 1]
        Vcar2[, h]     <- alpha * Vcar2[, aux]   + (1 - alpha) * Bc * car2[, 1]
      }
    }
  }
  #print(Vtrain2[1,])
  for (h in 1:H2) {
    auxt <- 2 * h
    auxc <- 2*h + 1
    Vst2[, h] <-(1 / ue) * log(exp(ue * Vtrain2[, auxt]) + exp(ue * Vsm2[, auxt]) + exp(ue *Vcar2[, auxt]) )
    Vsc2[, h] <-(1 / ue) * log(exp(ue * Vtrain2[, auxc]) + exp(ue * Vsm2[, auxc]) + exp(ue *Vcar2[, auxc]) )
    
  }
  for (h in 1:H2) {
    Vs2[,h]<- (1 / us) * log(exp(us * Vst2[,h]) + exp(us * Vsc2[,h]))
    Ve2[,h]<- (1 / ue) * log(exp(ue * Vtrain2[,h]) + exp(ue * Vsm2[,h]) + exp(ue * Vcar2[,h]) )
  }
  for (h in 1:H2) {
    t <- log2(h) %/% 1 
    desc <- desct*((t)^2)
    PS2[,h]<- (exp(u * Vs2[,h]) / (exp(u*Vs2[,h]) + exp(u * (Ve2[,h] + desc))) ) + (exp(u * Ve2[,h]) / (exp(u * (Vs2[,h]  + desc )) + exp(u * Ve2[,h])))
    
  }
  
  PS2[PS2>1]<-1
  
  PE2 <- 1 - PS2
  PE2[is.na(PE2)] <- 1
  PS2[is.na(PS2)] <- 0
  
  K_t<-beta_t
  K_c<-0
  denom_K<- exp(us*K_t)+exp(us*K_c)
  PT2  <- exp(us*K_t)/denom_K
  PC2  <- exp(us*K_c)/denom_K
  
  denom2 <- exp(ue*Vtrain2) + exp(ue*Vsm2) + exp(ue*Vcar2)
  P12 <- exp(ue*Vtrain2)  / denom2
  P22 <- exp(ue*Vsm2) / denom2
  P32 <- exp(ue*Vcar2) / denom2
  
  for (h in 1:H) {
    aux <- h %/% 2
    h2 <- h
    camino <- 1
    while (aux > 0) {
      if (h2 %% 2 == 0) {
        camino <- camino * PS2[, aux] * PT2
        h2 <- aux
        aux <- aux %/% 2
      }
      else{
        camino <- camino * PS2[, aux] * PC2
        h2 <- aux
        aux <- aux %/% 2
      }
    }
    CAMINO2[,h]  <- camino
  }
  CAMINO2<-CAMINO2*PE2
  Prtrainm2  <-  P12 * CAMINO2
  Prsmm2     <-  P22 * CAMINO2
  Prcarm2    <-  P32 * CAMINO2
  
  Prtrain2  <- rowSums(Prtrainm2)
  Prsm2     <- rowSums(Prsmm2)
  Prcar2    <- rowSums(Prcarm2)
  
  Prtrain2[is.na( Prtrain2)] <- 0
  Prsm2[is.na( Prsm2)] <- 0
  Prcar2[is.na( Prcar2)] <- 0
  
  lli2  <- log((Prtrain2 * choice[, 1]) + (Prsm2 * choice[, 2]) + (Prcar2 * choice[, 3]))
  
  
  return(-sum(lli2))  
}