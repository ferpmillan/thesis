## Purpose: Replicates Chapter 5 copula-based multivariate modeling and backtesting.
## Inputs: Equity and bond return series with Chapter 3/4 helper functions.
## Dependencies: Helper functions sourced in scripts/00_run_all.R.

##### Chapter 5. Copula Multivariate Model. ##### 

# Prepare bond returns
bond_returns <- xts(final_bond_df$Log_Return, order.by = as.Date(final_bond_df$Date))
bond_returns <- na.omit(bond_returns)

# Define train and test periods
train_start <- as.Date("2014-01-01")
train_end   <- as.Date("2024-12-31")
test_start  <- as.Date("2025-01-01")
test_end    <- as.Date("2025-04-30")

# Split data
bond_train <- bond_returns[paste0(train_start, "/", train_end)]
bond_test  <- bond_returns[paste0(test_start, "/", test_end)]

# Combine train + test
returns_full_bond <- c(bond_train, bond_test)

# Number of observations
n_total_bond <- length(returns_full_bond)
n_train_bond <- length(bond_train)
n_test_bond  <- length(bond_test)

# GARCH(1,1).
garch_spec_bond <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
garch_fit_bond_full <- ugarchfit(
  spec = garch_spec_bond,
  data = returns_full_bond,
  out.sample = n_test_bond
)

# GJR-GARCH(1,1). 
gjr_spec_bond <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
gjr_fit_bond_full <- ugarchfit(
  spec = gjr_spec_bond,
  data = returns_full_bond,
  out.sample = n_test_bond,
  solver = "hybrid",
  fit.control = list(trace = TRUE)
)


# Extract GARCH(1,1) parameters
garch_params_bond <- coef(garch_fit_bond_full)
cat("GARCH(1,1) Parameters (Bond):\n")
print(garch_params_bond)


# Extract GJR(1,1) parameters
gjr_params_bond <- coef(gjr_fit_bond_full)
cat("\nGJR-GARCH(1,1) Parameters (Bond):\n")
print(gjr_params_bond)

# Extract standardized retrurns 
z_garch_bond <- residuals(garch_fit_bond_full, standardize = TRUE)
z_gjr_bond   <- residuals(gjr_fit_bond_full, standardize = TRUE)

# Asymmetric t via MLE
fit_garch_asyt_bond <- optim(
  par = c(6, 0.5),
  fn = loglik_asyt,
  z = as.numeric(z_garch_bond),
  method = "L-BFGS-B",
  lower = c(2.01, -1.01),
  upper = c(100, 0.99)
)

fit_gjr_asyt_bond <- optim(
  par = c(6, 0.5),
  fn = loglik_asyt,
  z = as.numeric(z_gjr_bond),
  method = "L-BFGS-B",
  lower = c(2.01, -1.01),
  upper = c(100, 0.99)
)

# Bond parameters
d1_garch_bond <- fit_garch_asyt_bond$par[1]
d2_garch_bond <- fit_garch_asyt_bond$par[2]

d1_gjr_bond <- fit_gjr_asyt_bond$par[1]
d2_gjr_bond <- fit_gjr_asyt_bond$par[2]

# Copula inputs (PIT's)

# Figure 5.2: Threshold Correlation of GARCH(1,1) standardized returns (training period)

z1 <- na.omit(z_garch)           # S&P 500 standardized residuals
z2 <- na.omit(z_garch_bond)          # Bond standardized residuals

aligned <- na.omit(merge(
  xts(z1, order.by = index(sp500_log_returns)),
  xts(z2, order.by = index(bond_train))
))

colnames(aligned) <- c("z1", "z2")

quantiles <- seq(0.1, 0.9, by = 0.01)

threshold_corr <- sapply(quantiles, function(p) {
  q1 <- quantile(aligned$z1, p)
  q2 <- quantile(aligned$z2, p)
  
  if (p <= 0.5) {
    mask <- aligned$z1 <= q1 & aligned$z2 <= q2
  } else {
    mask <- aligned$z1 > q1 & aligned$z2 > q2
  }
  
  if (sum(mask) > 25) {
    cor(aligned$z1[mask], aligned$z2[mask])
  } else {
    NA  # NA if not enough observations
  }
})

df_threshold <- data.frame(Quantile = quantiles, Correlation = threshold_corr)

ggplot(df_threshold, aes(x = Quantile, y = Correlation)) +
  geom_line(color = "deepskyblue3", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Threshold correlation", x = "Quantile") +
  theme_minimal(base_family = "Times")+
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),  
        text = element_text(size = 18)) 

# Figure 5.1: 25-Day rolling correlation of standardized returns (training period). 
roll_cor <- rollapply(
  aligned, width = 25, FUN = function(x) cor(x[,1], x[,2]),
  by.column = FALSE, align = "right", fill = NA
)

df_cor <- data.frame(Date = index(roll_cor), Covariance = as.numeric(roll_cor))

ggplot(df_cor, aes(x = Date, y = Covariance)) +
  geom_line(color = "deepskyblue3", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs( y = "Correlation", x = "Date") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),  
        text = element_text(size = 18)) 




# Copula modeling 

pit_garch_sp500 <- xts(pasyt(as.numeric(z_garch), d1 = d1_garch_asyt, d2 = d2_garch_asyt),
                       order.by = index(z_garch))
pit_gjr_sp500 <- xts(pasyt(as.numeric(z_gjr), d1 = d1_gjr_asyt, d2 = d2_gjr_asyt),
                     order.by = index(z_gjr))

pit_garch_bond <- xts(pasyt(as.numeric(z_garch_bond), d1 = d1_garch_bond, d2 = d2_garch_bond),
                      order.by = index(z_garch_bond))
pit_gjr_bond <- xts(pasyt(as.numeric(z_gjr_bond), d1 = d1_gjr_bond, d2 = d2_gjr_bond),
                    order.by = index(z_gjr_bond))

pit_pair_garch <- merge(pit_garch_sp500, pit_garch_bond, join = "inner")

pit_clean <- pmax(pmin(pit_pair_garch, 1 - 1e-10), 1e-10)

u_data <- as.matrix(pit_clean)

t_cop <- tCopula(dim = 2, df = 4, df.fixed = FALSE, dispstr = "un")
fit_t <- fitCopula(t_cop, u_data, method = "ml")

summary(fit_t)

pit_pair_gjr <- merge(pit_gjr_sp500, pit_gjr_bond, join = "inner")

pit_clean_gjr <- pmax(pmin(pit_pair_gjr, 1 - 1e-10), 1e-10)

u_data_gjr <- as.matrix(pit_clean_gjr)

t_cop_gjr <- tCopula(dim = 2, df = 4, df.fixed = FALSE, dispstr = "un")
fit_t_gjr <- fitCopula(t_cop_gjr, u_data_gjr, method = "ml")

summary(fit_t_gjr)

copula_fit <- fit_t@copula

set.seed(191023)
u_sim <- rCopula(1000000, copula_fit)
set.seed(191023)

z_sim_sp500 <- qasyt(u_sim[, 1], d1 = d1_garch_asyt, d2 = d2_garch_asyt)
z_sim_bond  <- qasyt(u_sim[, 2], d1 = d1_garch_bond,  d2 = d2_garch_bond)

sim_df <- data.frame(SP500 = z_sim_sp500, Bond = z_sim_bond)

thresh_corr_sim <- sapply(quantiles, function(p) {
  q1 <- quantile(sim_df$SP500, p)
  q2 <- quantile(sim_df$Bond,  p)
  
  if (p <= 0.5) {
    mask <- sim_df$SP500 <= q1 & sim_df$Bond <= q2
  } else {
    mask <- sim_df$SP500 > q1 & sim_df$Bond > q2
  }
  
  if (sum(mask) > 25) {
    cor(sim_df$SP500[mask], sim_df$Bond[mask])
  } else {
    NA
  }
})

df_sim_threshold <- data.frame(
  Quantile = quantiles,
  Correlation = thresh_corr_sim,
  Type = "Simulated"
)

threshold_plot_df <- rbind(
  cbind(df_threshold, Type = "Empirical"),
  df_sim_threshold
)


copula_fit_gjr <- fit_t_gjr@copula

set.seed(191023)
u_sim_gjr <- rCopula(1000000, copula_fit_gjr)
set.seed(191023)

z_sim_sp500_gjr <- qasyt(u_sim_gjr[, 1], d1 = d1_gjr_asyt, d2 = d2_gjr_asyt)
z_sim_bond_gjr  <- qasyt(u_sim_gjr[, 2], d1 = d1_gjr_bond,  d2 = d2_gjr_bond)


sim_df_gjr <- data.frame(SP500 = z_sim_sp500_gjr, Bond = z_sim_bond_gjr)

thresh_corr_sim_gjr <- sapply(quantiles, function(p) {
  q1 <- quantile(sim_df_gjr$SP500, p)
  q2 <- quantile(sim_df_gjr$Bond,  p)
  
  if (p <= 0.5) {
    mask <- sim_df_gjr$SP500 <= q1 & sim_df_gjr$Bond <= q2
  } else {
    mask <- sim_df_gjr$SP500 > q1 & sim_df_gjr$Bond > q2
  }
  
  if (sum(mask) > 25) {
    cor(sim_df_gjr$SP500[mask], sim_df_gjr$Bond[mask])
  } else {
    NA
  }
})

df_sim_gjr_threshold <- data.frame(
  Quantile = quantiles,
  Correlation = thresh_corr_sim_gjr,
  Type = "Simulated"
)

z1 <- na.omit(z_gjr)        # S&P 500
z2 <- na.omit(z_gjr_bond)   # Bond

aligned_gjr <- na.omit(merge(
  xts(z1, order.by = index(sp500_log_returns)),
  xts(z2, order.by = index(bond_train))
))

colnames(aligned_gjr) <- c("z1", "z2")

threshold_corr_gjr <- sapply(quantiles, function(p) {
  q1 <- quantile(aligned_gjr$z1, p)
  q2 <- quantile(aligned_gjr$z2, p)
  
  if (p <= 0.5) {
    mask <- aligned_gjr$z1 <= q1 & aligned_gjr$z2 <= q2
  } else {
    mask <- aligned_gjr$z1 > q1 & aligned_gjr$z2 > q2
  }
  
  if (sum(mask) > 25) {
    cor(aligned_gjr$z1[mask], aligned_gjr$z2[mask])
  } else {
    NA
  }
})

df_threshold_gjr <- data.frame(
  Quantile = quantiles,
  Correlation = threshold_corr_gjr
)


threshold_plot_df_gjr <- rbind(
  cbind(df_threshold_gjr, Type = "Empirical"),
  df_sim_gjr_threshold
)

# Figure 5.3: Empirical and simulated threshold correlations of GARCH-GARCH copula.
ggplot(threshold_plot_df, aes(x = Quantile, y = Correlation, color = Type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Empirical" = "deepskyblue3", "Simulated" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Threshold correlation", x = "Quantile") +
  theme_minimal(base_family = "Times") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = "right")

# Figure 5.4: Empirical and simulated threshold correlations of GJR-GJR copula.
ggplot(threshold_plot_df_gjr, aes(x = Quantile, y = Correlation, color = Type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Empirical" = "deepskyblue3", "Simulated" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Threshold correlation", x = "Quantile") +
  theme_minimal(base_family = "Times") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = "right")

# Forecasting using the Copula Multivariate Model

garch_forecast_bond <- ugarchforecast(garch_fit_bond_full, n.ahead = 1, n.roll = n_test-1)
gjr_forecast_bond   <- ugarchforecast(gjr_fit_bond_full,   n.ahead = 1, n.roll = n_test - 1)

M <- 10000 #MC simulations per day 
n_test <- length(sigma(garch_forecast))  # S&P test length
dates_test <- as.Date(index(sp500_test_returns))
VaR_025_garch <- numeric(n_test)
VaR_025_gjr   <- numeric(n_test)
PIT_garch <- numeric(n_test)
PIT_gjr   <- numeric(n_test)
ret_sp500 <- as.numeric(sp500_test_returns)
ret_bond  <- as.numeric(bond_test[dates_test])

# Simulation
for (i in 1:n_test) {
  
  # Step 1: Conditional volatilities
  sigma_sp_garch  <- sigma(garch_forecast)[i]
  sigma_bond_garch<- sigma(garch_forecast_bond)[i]
  
  sigma_sp_gjr    <- sigma(gjr_forecast)[i]
  sigma_bond_gjr  <- sigma(gjr_forecast_bond)[i]
  
  # Step 2: Simulated PITs
  set.seed(191023)
  u_sim_garch <- rCopula(M, fit_t@copula)
  u_sim_gjr   <- rCopula(M, fit_t_gjr@copula)
  
  # Step 3: Transform to shocks
  z_sp_garch  <- qasyt(u_sim_garch[,1], d1_garch_asyt, d2_garch_asyt)
  z_bond_garch<- qasyt(u_sim_garch[,2], d1_garch_bond, d2_garch_bond)
  
  z_sp_gjr    <- qasyt(u_sim_gjr[,1], d1_gjr_asyt, d2_gjr_asyt)
  z_bond_gjr  <- qasyt(u_sim_gjr[,2], d1_gjr_bond, d2_gjr_bond)
  
  # Step 4: Scale to returns
  r_sp_garch   <- sigma_sp_garch * z_sp_garch
  r_bond_garch <- sigma_bond_garch * z_bond_garch
  
  r_sp_gjr     <- sigma_sp_gjr * z_sp_gjr
  r_bond_gjr   <- sigma_bond_gjr * z_bond_gjr
  
  # Step 5: Simulated portfolio returns
  port_ret_garch <- 0.5 * r_sp_garch + 0.5 * r_bond_garch
  port_ret_gjr   <- 0.5 * r_sp_gjr   + 0.5 * r_bond_gjr
  
  # Step 6: Compute 2.5% VaR
  VaR_025_garch[i] <- quantile(port_ret_garch, 0.025)
  VaR_025_gjr[i]   <- quantile(port_ret_gjr, 0.025)
  
  # Compute empirical PIT (fraction of sims <= actual return)
  actual_return <- 0.5 * ret_sp500[i] + 0.5 * ret_bond[i]
  PIT_garch[i] <- mean(port_ret_garch <= actual_return)
  PIT_gjr[i]   <- mean(port_ret_gjr   <= actual_return)
  
  # Progress
  if (i %% 10 == 0) cat("Day", i, "of", n_test, "done\n")
  set.seed(191023)
}


df_VaR <- data.frame(
  Date             = dates_test,
  VaR_GARCH        = VaR_025_garch,
  VaR_GJR          = VaR_025_gjr,
  PIT_GARCH        = PIT_garch,
  PIT_GJR          = PIT_gjr,
  Return_SP500     = ret_sp500,
  Return_Bond      = ret_bond
)

df_VaR$Portfolio_Return <- 0.5 * df_VaR$Return_SP500 + 0.5 * df_VaR$Return_Bond

df_VaR$Breach_GARCH <- as.integer(df_VaR$Portfolio_Return < df_VaR$VaR_GARCH)
df_VaR$Breach_GJR   <- as.integer(df_VaR$Portfolio_Return < df_VaR$VaR_GJR)
df_VaR <- na.omit(df_VaR)

# Figure 5.6: 2.5% VaR breaches, GARCH-GACRH Copula. 
ggplot(df_VaR, aes(x = Date)) +
  geom_line(aes(y = Portfolio_Return), color = "black") +
  geom_line(aes(y = VaR_GARCH), color = "deepskyblue3", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = subset(df_VaR, Breach_GARCH == 1),
    aes(y = Portfolio_Return),
    color = "blue",
    fill = "white",
    size = 2,
    shape = 21,
    stroke = 1
  ) +
  labs(
    x = "Date",
    y = "Portfolio Return"
  ) +
  theme_tufte(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    text = element_text(size = 18)
  )

# Figure 5.7: 2.5% VaR breaches, GJR-GJR Copula. 
ggplot(df_VaR, aes(x = Date)) +
  geom_line(aes(y = Portfolio_Return), color = "black") +
  geom_line(aes(y = VaR_GJR), color = "seagreen4", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = subset(df_VaR, Breach_GJR == 1),
    aes(y = Portfolio_Return),
    color = "green4",
    fill = "white",
    size = 2,
    shape = 21,
    stroke = 1
  ) +
  labs(
    x = "Date",
    y = "Portfolio Return"
  ) +
  
  theme_tufte(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    text = element_text(size = 18)
  )

# Figure 5.8: Histogram of PIT values, GARCH-GARCH copula.
ggplot(df_VaR, aes(x = PIT_GARCH)) +
  geom_histogram(aes(y = ..density..), bins = 20, 
                 fill = "deepskyblue3", color = "black", alpha = 1) +
  labs(
    x = expression(u[t]),
    y = "Density") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 5.9: Histogram of PIT values, GJR-GJR copula.
ggplot(df_VaR, aes(x = PIT_GJR)) +
  geom_histogram(aes(y = ..density..), bins = 20, 
                 fill = "deepskyblue3", color = "black", alpha = 1) +
  labs(
    x = expression(u[t]),
    y = "Density") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Transform PITs to standard normal scale
z_pit_garch <- qnorm(df_VaR$PIT_GARCH)
z_pit_gjr   <- qnorm(df_VaR$PIT_GJR)

# Kolmogorov-Smirnov tests
ks_garch <- ks.test(z_pit_garch, "pnorm")
ks_gjr   <- ks.test(z_pit_gjr, "pnorm")

print(ks_garch)
print(ks_gjr)

# Function to prepare data for plotting

df_garch_ks <- get_ks_plot_data(z_pit_garch)
df_gjr_ks   <- get_ks_plot_data(z_pit_gjr)

# Figure 5.10: Kolmogorov-Smirnov test on ˆzt+1, GARCH-GARCH copula.
ggplot(df_garch_ks, aes(x = z)) +
  geom_step(aes(y = ECDF, color = "Empirical CDF"), linewidth = 1) +
  geom_line(aes(y = NormalCDF, color = "Standard Normal CDF"), 
            linewidth = 1, linetype = "dashed") +
  geom_segment(aes(x = max_z, xend = max_z,
                   y = max_ecdf, yend = max_norm),
               color = "black", linewidth = 1.2) +
  geom_point(aes(x = max_z, y = max_ecdf), shape = 21, size = 2,
             fill = "deepskyblue3", color = "black", stroke = 0.6) +
  geom_point(aes(x = max_z, y = max_norm), shape = 21, size = 2,
             fill = "red", color = "black", stroke = 0.6) +
  scale_color_manual(values = c(
    "Empirical CDF" = "deepskyblue3",
    "Standard Normal CDF" = "red"
  )) +
  labs(
    x = expression(hat(z)[t]),
    y = "CDF",
    color = "Legend") +
  theme_tufte(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    text = element_text(size = 18),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

# Figure 5.11: Kolmogorov-Smirnov test on ˆzt+1, GJR-GJR copula.
ggplot(df_gjr_ks, aes(x = z)) +
  geom_step(aes(y = ECDF, color = "Empirical CDF"), linewidth = 1) +
  geom_line(aes(y = NormalCDF, color = "Standard Normal CDF"), 
            linewidth = 1, linetype = "dashed") +
  geom_segment(aes(x = max_z, xend = max_z,
                   y = max_ecdf, yend = max_norm),
               color = "black", linewidth = 1.2) +
  geom_point(aes(x = max_z, y = max_ecdf), shape = 21, size = 2,
             fill = "deepskyblue3", color = "black", stroke = 0.6) +
  geom_point(aes(x = max_z, y = max_norm), shape = 21, size = 2,
             fill = "red", color = "black", stroke = 0.6) +
  scale_color_manual(values = c(
    "Empirical CDF" = "deepskyblue3",
    "Standard Normal CDF" = "red"
  )) +
  labs(
    x = expression(hat(z)[t]),
    y = "CDF",
    color = "Legend") +
  theme_tufte(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    text = element_text(size = 18),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

# VaR tests
cc_copula_garch <- christoffersen_test(df_VaR$Breach_GARCH, alpha = 0.025)
cc_copula_gjr <- christoffersen_test(df_VaR$Breach_GJR, alpha = 0.025)

print(cc_copula_garch)
print(cc_copula_gjr)
