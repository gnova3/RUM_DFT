# P-RRM

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



apollo_control = list(
  modelName  ="P_RRM",
  modelDescr ="P-RRM model",
  indivID    ="id"
)

database = base

###Parameters to be estimated and their starting values
apollo_beta=c(asc_train    = 0,
              asc_car      = 0,
              asc_sm       = 0,
              B_tt = 0,
              B_tc = 0)

apollo_fixed = c("asc_sm")
apollo_inputs = apollo_validateInputs()

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P_cost1 = ( - 1 / 1 ) * ( pmax( 0 , SM_TC    - TRAIN_TC ) + pmax( 0 , CAR_TC - TRAIN_TC ) )
  P_cost2 = ( - 1 / 1 ) * ( pmax( 0 , TRAIN_TC - SM_TC )    + pmax( 0 , CAR_TC - SM_TC ) )
  P_cost3 = ( - 1 / 1 ) * ( pmax( 0 , TRAIN_TC - CAR_TC )   + pmax( 0 , SM_TC  - CAR_TC )  )
  
  P_time1 = ( - 1 / 1 ) * ( pmin( 0 , SM_TT    - TRAIN_TT ) + pmin( 0 , CAR_TT - TRAIN_TT ) )
  P_time2 = ( - 1 / 1 ) * ( pmin( 0 , TRAIN_TT - SM_TT )    + pmin( 0 , CAR_TT - SM_TT )  )
  P_time3 = ( - 1 / 1 ) * ( pmin( 0 , TRAIN_TT - CAR_TT )   + pmin( 0 , SM_TT  - CAR_TT ) )

  ### List of regret functions
  R = list()
  R[['Alt1']]  =  asc_train + B_tc * P_cost1 + B_tt * P_time1 
  R[['Alt2']]  =  asc_sm    + B_tc * P_cost2 + B_tt * P_time2 
  R[['Alt3']]  =  asc_car   + B_tc * P_cost3 + B_tt * P_time3 
  
  mnl_settings = list(
    alternatives  = c(Alt1=1, Alt2=2, Alt3=3), 
    avail         = list(Alt1=1, Alt2=1, Alt3=1),
    choiceVar     = C,
    V             = R
  )
  
  P = list()
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)  
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

                    
model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)

# mu-RRM


apollo_control = list(
  modelName  ="mu-RRM",
  modelDescr ="mu-RRM model",
  indivID    ="id"
)

database = base

###Parameters to be estimated and their starting values
apollo_beta=c(asc_train    = 0,
              asc_car      = 0,
              asc_sm       = 0,
              B_tt = 0,
              B_tc = 0, 
              mu = 1)

apollo_fixed = c("asc_sm")
apollo_inputs = apollo_validateInputs()

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  P_cost2_1 = (SM_TC    - TRAIN_TC ) 
  P_cost3_1 = (CAR_TC    - TRAIN_TC ) 
  
  P_cost1_2 = (TRAIN_TC    - SM_TC ) 
  P_cost3_2 = (CAR_TC    - SM_TC ) 
  
  P_cost1_3 = (TRAIN_TC - CAR_TC) 
  P_cost2_3 = (SM_TC - CAR_TC) 
  
  
  
  P_time2_1 = (SM_TT    - TRAIN_TT ) 
  P_time3_1 = (CAR_TT    - TRAIN_TT ) 
  
  P_time1_2 = (TRAIN_TT    - SM_TT ) 
  P_time3_2 = (CAR_TT    - SM_TT ) 
  
  P_time1_3 = (TRAIN_TT - CAR_TT) 
  P_time2_3 = (SM_TT - CAR_TT) 

  
  ### List of regret functions
  R = list()

  R[['Alt1']]  = asc_train + mu * ( - log( 1 + exp( ( B_tt / mu ) * P_time2_1 ) ) - log( 1 + exp( ( B_tt / mu ) * P_time3_1 ) ) -  log( 1 + exp( ( B_tc / mu ) * P_cost2_1 ) ) - log( 1 + exp( ( B_tt / mu ) * P_cost3_1 ) ) )
  R[['Alt2']]  = asc_sm    + mu * ( - log( 1 + exp( ( B_tt / mu ) * P_time1_2 ) ) - log( 1 + exp( ( B_tt / mu ) * P_time3_2 ) ) -  log( 1 + exp( ( B_tc / mu ) * P_cost1_2 ) ) - log( 1 + exp( ( B_tt / mu ) * P_cost3_2 ) ) )
  R[['Alt3']]  = asc_car   + mu * ( - log( 1 + exp( ( B_tt / mu ) * P_time1_3 ) ) - log( 1 + exp( ( B_tt / mu ) * P_time2_3 ) ) -  log( 1 + exp( ( B_tc / mu ) * P_cost1_3 ) ) - log( 1 + exp( ( B_tt / mu ) * P_cost2_3 ) ) )
  
  
  mnl_settings = list(
    alternatives  = c(Alt1=1, Alt2=2, Alt3=3), 
    avail         = list(Alt1=1, Alt2=1, Alt3=1),
    choiceVar     = C,
    V             = R
  )
  
  P = list()
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)  
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}


model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)
