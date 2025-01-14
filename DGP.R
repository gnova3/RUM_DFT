rm(list = ls())
source("functions.R")

iterations = 20

for (iteration  in 1:iterations) {

    ## Define the parameters 
  DGP <- list(B_time = -2,
              B_cost = -0.5,
              UA_0   = -0.5,
              UB_0   = -0.5,
              UC_0   = 0,
              UD_0   = -1,
              Wc     = 0.3,
              Wt     = 0,
              delta  = 0.05,
              alpha  = 0.6,
              mu_1    = 0.1,
              mu_2   = 0.1)
  
  # Define the number of individuals and maximum number of searchs
  N     <- 4000
  TMAX  <- 10
  TTmax <- 30   
  
  # Preconceived utilities for each alternative
  UA_0 <- rep(DGP$UA_0, N)
  UB_0 <- rep(DGP$UB_0, N)
  UC_0 <- rep(DGP$UC_0, N)
  UD_0 <- rep(DGP$UD_0, N)
  
  # Generate random attribute values for each individual
  TVA <- abs(runif(N, 10, 20))/10 # Bici
  TVB <- abs(runif(N, 3, 22))/10 # micro
  TVC <- abs(runif(N, 4, 15))/10 #metro
  TVD <- abs(runif(N, 2, 25))/10 # auto
  
  # Generate a binary population variable and # assign costs based on population (B and C)
  poblacion <- rbinom(N, 1, 0.8)
  COSTA <- runif(N, 0, 500)/100
  COSTB <- ifelse(poblacion == 1, 700/100, 230/100)
  COSTC <- ifelse(poblacion == 1, 800/100, 230/100)
  COSTD <- runif(N, 0, 1800)/100  
  
  
  TDV  <- cbind(TVA, TVB, TVC, TVD)
  COST <- cbind(COSTA, COSTB, COSTC, COSTD)
  
  ################### ISP ######################
  # Initialize variables and matrices
  H      = vector(length=TTmax) * 0
  Choice = vector(length=N) * 0
  T_n    = vector(length=N) * 0
  
  # Loop through individuals
  for (n in 1:N) {                               
    TDVn  <- TDV[n,]   # Vector for attribute 1
    COSTn <- COST[n,]  # Vector for attribute 2
    
    UA = vector(length=TTmax) * 0  # Vector to store the utility for alternative A
    UB = vector(length=TTmax) * 0  # Vector to store the utility for alternative B
    UC = vector(length=TTmax) * 0  # Vector to store the utility for alternative C
    UD = vector(length=TTmax) * 0  # Vector to store the utility for alternative D
    
    hn = vector(length=TTmax) * 0  # Vector to store the information search process for each individual
    
    t = 1
    # Loop through steps till maximun valid value (TTMAX =30)
    while (t < TTmax) { 
      
      #cat("Timestap:", t, "\n")
      
      # Initialize the utility of each alternative with preconceived utilities
      if (t == 1) {  
        #cat("Preconcebided utilities initialised", "\n")
        VA <- UA_0[n]  
        VB <- UB_0[n] 
        VC <- UC_0[n] 
        VD <- UD_0[n]  
      }
      
      # Update the threshold
      delta = DGP$delta * t^2
      
      # Call the kernel function to get probabilities at step t
      K = kernel(VA, VB, VC, VD, DGP$mu_1, DGP$mu_2, DGP$alpha, DGP$Wc, DGP$Wt, COSTn, TDVn, delta, DGP$B_cost, DGP$B_time)
      #cat("Kernel calculated with V", t, "considering all information till", t, "\n")
      
      # Calculate probabilities
      Pac_1 <- cumsum(c(K[1], K[2]))      # Probability of choosing or searching at step t
      Pac_2 <- cumsum(K[3:6])             # Probability of choosing an alternative (A, B, C, D)
      Pac_3 <- cumsum(K[7:8])             # Probability of choosing an attribute (Att_1, Att_2)
      
      # Generate random numbers
      Random_1 <- runif(1, min = 0, max = 1)
      Random_2 <- runif(1, min = 0, max = 1)
      Random_3 <- runif(1, min = 0, max = 1)
      
      # If decides to choose, simulate the choice of an alternative
      if (Random_1 < Pac_1[1]) {  
        
        if (Random_2 < Pac_2[1]) { 
          Choice[n] <- 1
          #cat("Decide to choose alternative", 1, "at timestamp" , t, "\n")
        } else if (Random_2 > Pac_2[1] & Random_2 < Pac_2[2]) { 
          Choice[n] <- 2
          #cat("Decide to choose alternative", 2, "at timestamp" , t, "\n")
        } else if (Random_2 > Pac_2[2] & Random_2 < Pac_2[3]) {
          Choice[n] <- 3
          #cat("Decide to choose alternative", 3, "at timestamp" , t, "\n")
        } else if (Random_2 > Pac_2[3]) {
          Choice[n] <- 4
          #cat("Decide to choose alternative", 4, "at timestamp" , t, "\n")
        }
        
        # choice step is defined as the current step and the loop is stopped
        T_n[n] = t
        t = TTmax
        #cat("End of the loop", "\n")
        
      } 
      
      # If decides to attend an attribute, simulate the attention to an attribute
      else if (Random_1 >= Pac_1[1]) {  
        
        # If Att_1 is attended, the decision is stored and the utilities are updated with the attribute 1.
        if (Random_3 < Pac_3[1]) {
          #cat("Decide to search for information: Cost at timestamp" , t, "\n")
          #cat("With this new information, the utilities are updated", "\n")
          
          hn[t] = 5 
          u_next = update(VA, VB, VC, VD, DGP$alpha, COSTn, DGP$B_cost)
          #cat("The information search path is", hn, "\n")
        } 
        # If time is attended, the decision is stored and the utilities are updated with the attribute time.
        else if (Random_3 > Pac_3[1]) {
          #cat("Decide to search for information: Cost at timestamp" , t, "\n")
          #cat("With this new information, the utilities are updated", "\n")
          hn[t] = 6
          u_next = update(VA, VB, VC, VD, DGP$alpha, TDVn, DGP$B_time)
          #cat("The information search path is", hn, "\n")
        }
        
        # Update the utilities with the new values
        VA <- u_next[1]
        VB <- u_next[2]
        VC <- u_next[3]
        VD <- u_next[4]
        
      }
      
      # Update the step
      t = t + 1
    }
    
    # Store the information search path for each individual
    H = rbind(H, hn)  
  }
  
  H <- H[-1,]
  
  # Create vectors for each alternative choice
  c1 <- rep(0, N)
  c2 <- rep(0, N)
  c3 <- rep(0, N)
  c4 <- rep(0, N)
  
  c1[Choice == 1] <- 1
  c2[Choice == 2] <- 1
  c3[Choice == 3] <- 1
  c4[Choice == 4] <- 1
  
  # Combine vectors into a matrix
  choice <- cbind(c1, c2, c3, c4)
  
  # Create a data frame with simulated choice data
  choicedata_RUM_DFT <- data.frame("ID"=1:N, TDV, COST, "tiempoxid"=T_n, c1, c2, c3, c4)
  ChoiceData  = cbind(choicedata_RUM_DFT, Choice)
  
  #### Save ####
  write.csv(ChoiceData, file = paste("new_experiment/",paste(iteration,"choice_set.csv", sep = "_"), sep = ""), row.names = FALSE)
  write.csv(H,   file = paste("new_experiment/",paste(iteration,"H_set.csv", sep = "_"), sep = ""), row.names = FALSE)
}

