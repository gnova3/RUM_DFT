rm(list = ls())
source("functions.R")

iterations = 20

for (iteration  in 1:iterations) {
  set.seed(123 + iteration) 
  
  ##### read database ######
  D_path = paste("Repo/",paste(iteration,"choice_set.csv", sep = "_"), sep = "")
  H_path = paste("Repo/",paste(iteration,"H_set.csv", sep = "_"), sep = "")
  
  H = read.csv2(H_path, sep = ",")
  ChoiceData = read.csv2(D_path, sep = ",")
  choicedata_RUM_DFT = read.csv2(D_path, sep = ",")
  
  #### Estimation #####
  init  = c (-1,-.1,.1,-.1,-.1,-.5,0.1,0.01, 0.1, 0.1)
  result <- optim(par =init , fn = LL_total, hessian = TRUE, method = "BFGS", control = list(maxit = 3000, trace = TRUE, REPORT = 1)  )
  
  ### Results ####
  True  = c (-2,-.5,.6,-.5,-.5,-1,0.3,0.05, 0.1, 0.1)
  parameter_estimated = result$par
  standard_error = sqrt(diag(solve(result$hessian)))
  t_statistic_b = (parameter_estimated-True)/standard_error
  LL = result$value
  estimates <- c(parameter_estimated, standard_error, t_statistic_b, rep(0, 9), LL)
  result_table <- matrix(estimates, ncol = 4)
  result_table <- cbind(True, result_table)
  colnames(result_table) <- c("True", "Estimate", "S.E.", "t(b)", "LL")
  rownames(result_table) <- c("b_t", "b_c", "alpha", 'ASC_A', 'ASC_B', "ASC_D", "Wc", "delta", "mu1", "mu2")
  result_table
  
  write.csv2(result_table, paste("Repo/",paste(iteration,"FINAL_test_seed_repo_final_RUM_DFT_ISP_results", sep = "_"), sep = ""))
  
}
  
#### Outcomes ####

install.packages("tidyverse")
library(tidyverse)


files_path <- "Repo"  
csv_files <- list.files(files_path, pattern = "FINAL_test_seed_repo_final_RUM_DFT_ISP_results\\.csv$", full.names = TRUE)

estimates_list <- list()

for (i in seq_along(csv_files)) {
  current_file <- read.csv(csv_files[i], sep = ";")
  parameter_names <- current_file$X
  estimates <- current_file$Estimate
  
  estimates <- as.numeric(gsub(",", ".", estimates))
  
  ll_value <- as.numeric(gsub(",", ".", current_file$LL[nrow(current_file)]))
  
  estimates_named <- setNames(estimates, parameter_names)
  estimates_named["LL"] <- ll_value
  
  estimates_list[[i]] <- estimates_named
}

estimates_df <- do.call(rbind, estimates_list)
summary_df <- as.data.frame(estimates_df)

summary_stats <- summary_df %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    std = ~ sd(.x, na.rm = TRUE)
  )))

summary_stats_renamed <- as.data.frame(t(summary_stats))  
colnames(summary_stats_renamed) <- c("Value")  

summary_stats_renamed$Metric <- ifelse(grepl("_mean$", rownames(summary_stats_renamed)), "_mean", "_std")
summary_stats_renamed$Parameter <- gsub("_mean$|_std$", "", rownames(summary_stats_renamed))

final_summary <- pivot_wider(
  summary_stats_renamed,
  names_from = Metric,
  values_from = Value
)

print(final_summary)
write.csv(final_summary, file = "Repo/summary.csv", row.names = FALSE)


  
  