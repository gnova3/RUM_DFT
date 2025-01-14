# Final RUM-DFT-DT

# Clear Environment
cat("\014")
rm(list = ls())
graphics.off()
#source("functions.R")

#Boolean flags to run DEoptim:
# Set to TRUE  to re-estimate the alpha parameter for TMAX = 7/10 using DEoptim.
# Set to FALSE to use the alpha value reported in the paper.

optimize_T7 <- TRUE  
optimize_T10 <- TRUE

# Libraries
library(dplyr)   
library(DEoptim)
library(optimParallel)
library(apollo)

#### Load data and preprocess ####
data <- read.delim("swissmetro.dat", sep = "\t") %>%
  filter(PURPOSE %in% c(1, 3), CHOICE != 0, CAR_AV != 0) %>%
  # scaled variables
  mutate(
    SM_COST = SM_CO * (GA == 0),
    TRAIN_COST = TRAIN_CO * (GA == 0),
    TRAIN_TT_SCALED = TRAIN_TT / 10,
    TRAIN_COST_SCALED = TRAIN_COST / 10,
    SM_TT_SCALED = SM_TT / 10,
    SM_COST_SCALED = SM_COST / 10,
    CAR_TT_SCALED = CAR_TT / 10,
    CAR_CO_SCALED = CAR_CO / 10
  )
n		<- nrow(data)
N = length(data$ID)
ncus <- n_distinct(data$ID)
choice <- model.matrix(~ as.factor(CHOICE) - 1, data = data)
train2<-cbind(data$TRAIN_COST_SCALED,data$TRAIN_TT_SCALED)
sm2   <-cbind(data$SM_COST_SCALED,   data$SM_TT_SCALED)
car2  <-cbind(data$CAR_CO_SCALED,  data$CAR_TT_SCALED)

# Base dataframe for Apollo
base <- data.frame(
  id = data$ID,
  TRAIN_TC = data$TRAIN_COST_SCALED/10,
  TRAIN_TT = data$TRAIN_TT_SCALED/10,
  SM_TC = data$SM_COST_SCALED/10,
  SM_TT = data$SM_TT_SCALED/10,
  CAR_TC = data$CAR_CO_SCALED/10,
  CAR_TT = data$CAR_TT_SCALED/10,
  C = data$CHOICE
)

#### RUM-DFT Log-Likelihood functions  ####
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


#### Estimation ####

# 1) Deoptim - Tmax = 7 
TMAX = 7
if (optimize_T7) {
  set.seed(333)
  DEctrl <- DEoptim.control(itermax = 300,trace=TRUE,parallelType=1,packages=c(),parVar=c("N","train2","sm2","car2","choice","TMAX"))
  mymleDEa_7<-DEoptim(LL_Swiss,lower=c(-5,-5,-2,-1,-1,-1,-4) , upper=c(5,5,2,1,1,1,4),DEctrl)
  alpha = abs(mymleDEa_7$optim$bestmem[6])
  
  beta_init_7 = c(mymleDEa_7$optim$bestmem[1], 
                  mymleDEa_7$optim$bestmem[2],
                  mymleDEa_7$optim$bestmem[3], 
                  mymleDEa_7$optim$bestmem[4], 
                  mymleDEa_7$optim$bestmem[5],
                  #mymleDEa_7$optim$bestmem[6],
                  mymleDEa_7$optim$bestmem[7])
  
}else{
  alpha = 0.953
  beta_init_7 = c(-3.714,   -3.755,   -0.944,    0.001,    0.297,  0.294)
}

# 1) Optim - Tmax = 7 
mymleOP_7 <- optim(par=beta_init_7, fn=LL_Swiss_fixed, hessian=TRUE, method='BFGS',control  = list(maxit=30000, trace=TRUE, REPORT=1))

#Table
RUM_DFT_Swiss_7_par <- c(alpha ,mymleOP_7$par[1:4], abs(mymleOP_7$par[5]),mymleOP_7$par[6] )
RUM_DFT_Swiss_7_se  <- c(NaN, sqrt(diag(solve(mymleOP_7$hessian[1:6,1:6]))))
tabla_RUM_DFT_Swiss_7<-round(cbind(RUM_DFT_Swiss_7_par,RUM_DFT_Swiss_7_se, RUM_DFT_Swiss_7_par/RUM_DFT_Swiss_7_se, c(rep(0,6),mymleOP_7$value)),3)
rownames(tabla_RUM_DFT_Swiss_7)<-c("alpha", "Bt","Bc","Utrain","Ucar","delta", "Phi_t")
colnames(tabla_RUM_DFT_Swiss_7)<-c("Estimate","SE","t(0)", "LL")

#### tabla_RUM_DFT_Swiss_7 ####
#       Estimate    SE    t(0)      LL
#alpha     0.953   NaN     NaN    0.00
#Bt       -3.274 0.386  -8.491    0.00
#Bc       -4.136 0.594  -6.969    0.00
#Utrain   -0.913 0.086 -10.612    0.00
#Ucar     -0.008 0.061  -0.126    0.00
#delta     0.257 0.084   3.068    0.00
#Phi_t     0.483 0.140   3.455 4234.83

# 2) Deoptim - Tmax = 10
TMAX = 10
if (optimize_T7) {
  set.seed(333)
  DEctrl <- DEoptim.control(itermax = 300,trace=TRUE,parallelType=1,packages=c(),parVar=c("N","train2","sm2","car2","choice","TMAX"))
  mymleDEa_10<-DEoptim(LL_Swiss,lower=c(-5,-5,-2,-1,-1,-1,-4) , upper=c(5,5,2,1,1,1,4),DEctrl)
  alpha = abs(mymleDEa_10$optim$bestmem[6])
  
  beta_init_10 = c( mymleDEa_10$optim$bestmem[1], 
                    mymleDEa_10$optim$bestmem[2],
                    mymleDEa_10$optim$bestmem[3], 
                    mymleDEa_10$optim$bestmem[4], 
                    mymleDEa_10$optim$bestmem[5],
                    #mymleDEa_10$optim$bestmem[6],
                    mymleDEa_10$optim$bestmem[7])
  
}else{
  alpha = 0.960
  beta_init_10 = c(-3.545,   -4.994,   -0.953,   -0.083,    0.244,   0.543)
  
}

mymleOP_10 <- optim(par=beta_init_10, fn=LL_Swiss_fixed, hessian=TRUE, method='BFGS',control  = list(maxit=30000, trace=TRUE, REPORT=1))

#Table
RUM_DFT_Swiss_10_par <- c(alpha ,mymleOP_10$par[1:4], abs(mymleOP_10$par[5]),mymleOP_10$par[6] )
RUM_DFT_Swiss_10_se  <- c(NaN, sqrt(diag(solve(mymleOP_10$hessian[1:6,1:6]))))
tabla_RUM_DFT_Swiss_10<-round(cbind(RUM_DFT_Swiss_10_par,RUM_DFT_Swiss_10_se, RUM_DFT_Swiss_10_par/RUM_DFT_Swiss_10_se, c(rep(0,6),mymleOP_10$value)),3)
rownames(tabla_RUM_DFT_Swiss_10)<-c("alpha", "Bt","Bc","Utrain","Ucar","delta", "Phi_t")
colnames(tabla_RUM_DFT_Swiss_10)<-c("Estimate","SE","t(0)", "LL")


##### tabla_RUM_DFT_Swiss_10 #####
#Estimate    SE    t(0)       LL
#alpha     0.960   NaN     NaN    0.000
#Bt       -3.804 0.447  -8.514    0.000
#Bc       -4.830 0.692  -6.978    0.000
#Utrain   -0.898 0.085 -10.598    0.000
#Ucar     -0.008 0.060  -0.130    0.000
#delta     0.258 0.084   3.087    0.000
#Phi_t     0.486 0.139   3.489 4234.439 