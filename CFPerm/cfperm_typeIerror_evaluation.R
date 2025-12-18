#cfperm type-I error evaluation: Linear Model & Mars Model:
source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)

#snr represent the signal-to-noise ratio
snr_seq <- c(0.14, 0.25, 0.42, 0.71, 1.22, 2.07, 3.52, 6)
corr_seq <- seq(0, 0.7, 0.1)
B <- 500
for(k in 1:8){
  power_df <- data.frame(
      snr_seq = snr_seq,
      power = numeric(8)
    )
  for(r in 1:8){
    power_list <- numeric(B)
    for(i in 1:B){
      df1 <- correlation_model_df(n = 1200, p = 20, snr = snr_seq[r], rho = corr_seq[k])
      df2 <- correlation_model_df(n = 800, p = 20, snr = snr_seq[r], rho = corr_seq[k])
      n_col <- ncol(df1)
      power_list[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    }
    power_df$power[r] <- mean(power_list)
  }
  write.csv(power_df, paste(k, "SNR_LM_cfperm_typeIerror.csv", sep = "_"))
}
