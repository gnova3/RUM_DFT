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


#### MNL ####
library(apollo)
apollo_initialise()

apollo_control = list(
  modelName       = "MNL",
  modelDescr      = "swiss",
  indivID         = "id",
  outputDirectory = "output"
)

database = base
apollo_beta=c(asc_train    = 0,
              asc_car      = 0,
              asc_sm       = 0,
              b_tt         = 0,
              b_cost       = 0)

apollo_fixed = c("asc_car")
apollo_inputs = apollo_validateInputs()

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  V = list()
  V[["train"]]  = asc_train  + b_tt  * TRAIN_TT  + b_cost * TRAIN_TC
  V[["car"]]    = asc_car    + b_tt  * CAR_TT    + b_cost * CAR_TC 
  V[["sm"]]     = asc_sm     + b_tt  * SM_TT     + b_cost * SM_TC    
  
  mnl_settings = list(
    alternatives  = c(train=1, car=3, sm=2), 
    avail         = list(train=1, car=1, sm=1), 
    choiceVar     = C,
    utilities     = V
  )
  
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)