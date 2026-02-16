## Purpose: Replicates Chapter 4 forecasting, VaR/ES evaluation, and stress tests.
## Inputs: Fitted models, distribution parameters, and Chapter 3 residual objects.
## Dependencies: Helper functions sourced in scripts/00_run_all.R.

##### Chapter 4. Forecasting and Evaluation. ######

# Define date range for the test set
start_test <- as.Date("2024-12-31")
end_test   <- as.Date("2025-05-01")

# Download data (SP 500 from January to April 2025)
getSymbols("^GSPC", from = start_test, to = end_test, src = "yahoo")

# Extract closing prices
sp500_test <- Cl(GSPC)

# Compute log returns
sp500_test_returns <- na.omit(diff(log(sp500_test)))

# Combine full dataset forforecasting
returns_full <- c(sp500_log_returns, sp500_test_returns)

# Total number of observations
n_total <- length(returns_full)
n_train <- length(sp500_log_returns)
n_test  <- length(sp500_test_returns)

# Re-fit with training set
garch_fit <- ugarchfit(spec = garch_spec, data = returns_full, out.sample = n_test)
gjr_fit   <- ugarchfit(spec = gjr_spec,   data = returns_full, out.sample = n_test)

# Proper rolling forecast
garch_forecast <- ugarchforecast(garch_fit, n.ahead = 1, n.roll = n_test - 1)
gjr_forecast   <- ugarchforecast(gjr_fit,   n.ahead = 1, n.roll = n_test - 1)

# Extract 1-step-ahead variances
garch_sigma_forecast <- sigma(garch_forecast)^2
gjr_sigma_forecast   <- sigma(gjr_forecast)^2


## ACF of squared standardized returns ##

# Squared standardized returns
z2_garch_train <- as.numeric(z_garch)^2
z2_gjr_train   <- as.numeric(z_gjr)^2

# Compute ACF values of the squared standardized returns
acf_garch <- acf(z2_garch_train, lag.max = 100, plot = FALSE)
acf_gjr   <- acf(z2_gjr_train,   lag.max = 100, plot = FALSE)

# Build dataframes for plotting
acf_df_garch <- data.frame(
  Lag = acf_garch$lag[-1],
  ACF = acf_garch$acf[-1]
)

acf_df_gjr <- data.frame(
  Lag = acf_gjr$lag[-1],
  ACF = acf_gjr$acf[-1]
)

# 95% confidence bounds of the autocorrelations (approx.)
n_garch <- length(z2_garch_train)
n_gjr   <- length(z2_gjr_train)
conf_bound_garch <- 1.96 / sqrt(n_garch)
conf_bound_gjr   <- 1.96 / sqrt(n_gjr)

# Add to DataFrames
acf_df_garch$Upper <- conf_bound_garch
acf_df_garch$Lower <- -conf_bound_garch

acf_df_gjr$Upper <- conf_bound_gjr
acf_df_gjr$Lower <- -conf_bound_gjr


# Figure 4.1: Autocorrelation of squared returns over variance (2014-2024), GARCH(1,1).
ggplot(acf_df_garch, aes(x = Lag, y = ACF)) +
  geom_line(color = "deepskyblue3", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(aes(yintercept = Upper), linetype = "dashed", color = "gray40") +
  geom_hline(aes(yintercept = Lower), linetype = "dashed", color = "gray40") +
  labs(
    x = "Lag Order",
    y = "Autocorrelation of squared shocks"
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.2: Autocorrelation of squared returns over variance (2014-2024), GJR-GARCH(1,1).
ggplot(acf_df_gjr, aes(x = Lag, y = ACF)) +
  geom_line(color = "deepskyblue3", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(aes(yintercept = Upper), linetype = "dashed", color = "gray40") +
  geom_hline(aes(yintercept = Lower), linetype = "dashed", color = "gray40") +
  labs(
    x = "Lag Order",
    y = "Autocorrelation of squared shocks"
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))


## QLIKE loss function ##
## One of the returns in the training set is 0, 
## to calculate the loss function we replace the standardized return with a small epsilon.


# Mean QLIKE
qlike_garch_mean <- mean(qlike(z2_garch_train), na.rm = TRUE)
qlike_gjr_mean   <- mean(qlike(z2_gjr_train), na.rm = TRUE)

# Report
cat("QLIKE Loss:\n")
cat(sprintf("  GARCH(1,1):      %.6f\n", qlike_garch_mean))
cat(sprintf("  GJR-GARCH(1,1):  %.6f\n", qlike_gjr_mean))


## Likelihood ratio test (goodness of fit)

# Log-likelihood for GARCH(1,1)
ll_garch_qmle <- sum(log(fasyt(z_garch, d1_garch_asyt, d2_garch_asyt)))

# Log-likelihood for GJR-GARCH(1,1)
ll_gjr_qmle   <- sum(log(fasyt(z_gjr, d1_gjr_asyt, d2_gjr_asyt)))

# LR test statistic
LR_stat_qmle <- -2 * (ll_garch_qmle - ll_gjr_qmle)
p_val_qmle   <- 1 - pchisq(LR_stat_qmle, df = 1)

cat("LR Test:\n")
cat("LR statistic =", round(LR_stat_qmle, 4), "\n")
cat("p-value =", round(p_val_qmle, 4), "\n")


## Asymmetric-t distribution evaluation (VaR and density forecast) ##

# Quantile levels
p_lower <- 0.025
p_upper <- 0.975

# Forecasted standard deviation (from sigma forecast)
garch_sd <- sigma(garch_forecast)
gjr_sd   <- sigma(gjr_forecast)

# Compute VaR for test period
VaR_025_garch <- garch_sd * qasyt(p_lower, d1_garch_asyt, d2_garch_asyt)
VaR_975_garch <- garch_sd * qasyt(p_upper, d1_garch_asyt, d2_garch_asyt)

VaR_025_gjr <- gjr_sd * qasyt(p_lower, d1_gjr_asyt, d2_gjr_asyt)
VaR_975_gjr <- gjr_sd * qasyt(p_upper, d1_gjr_asyt, d2_gjr_asyt)

# Expected shortfall of the asymmetric t


# Expected shortfall in the test period under the asymmetric t-distribution
ES_025_garch <- -garch_sd * es_asyt(p_lower, d1_garch_asyt, d2_garch_asyt)
ES_025_gjr   <- -gjr_sd   * es_asyt(p_lower, d1_gjr_asyt,   d2_gjr_asyt)

df_var <- data.frame(
  Date = index(sp500_test_returns),
  Returns = as.numeric(sp500_test_returns),
  VaR_025_GARCH = as.numeric(VaR_025_garch),
  VaR_975_GARCH = as.numeric(VaR_975_garch),
  VaR_025_GJR   = as.numeric(VaR_025_gjr),
  VaR_975_GJR   = as.numeric(VaR_975_gjr)
)

df_var$ES_025_GARCH <- as.numeric(ES_025_garch)
df_var$ES_025_GJR   <- as.numeric(ES_025_gjr)


# 
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black",linewidth=0.5) +
  geom_line(aes(y = VaR_025_GARCH), color = "deepskyblue3", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = VaR_975_GARCH), color = "deepskyblue3", linetype = "dashed", linewidth = 1) +
  labs(
       x = "Date", y = "Log return") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

#
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black",linewidth=0.5) +
  geom_line(aes(y = VaR_025_GJR), color = "seagreen4", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = VaR_975_GJR), color = "seagreen4", linetype = "dashed", linewidth = 1) +
  labs(
       x = "Date", y = "Log return") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Identify VaR breaches (lower tail only)
df_var$Breach_GARCH <- ifelse(df_var$Returns < df_var$VaR_025_GARCH, 1, 0)
df_var$Breach_GJR   <- ifelse(df_var$Returns < df_var$VaR_025_GJR,   1, 0)

# Figure 4.3: 2.5% VaR forecast and breaches in the test period, GARCH-Asyt.
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black") +
  geom_line(aes(y = VaR_025_GARCH), color = "deepskyblue3", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = subset(df_var, Breach_GARCH == 1),
    aes(y = Returns), color = "blue", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  labs(
       x = "Date", y = "Log return") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.4: 2.5% VaR forecast and breaches in the test period,  GJR-Asyt.
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black") +
  geom_line(aes(y = VaR_025_GJR), color = "seagreen4", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = subset(df_var, Breach_GJR == 1),
    aes(y = Returns), color = "green4", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  labs(
       x = "Date", y = "Log return") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.5: 2.5% Expected Shortfall forecast in the test period,  GARCH-Asyt.
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black",linewidth=0.5) +
  geom_line(aes(y = VaR_025_GARCH), color = "deepskyblue3", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = ES_025_GARCH), color = "deepskyblue3", linetype = "solid", linewidth = 0.7) +
  geom_point(
    data = subset(df_var, Breach_GARCH == 1),
    aes(y = Returns),
    color = "blue", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  
  labs(
    x = "Date", y = "Log return"
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.6: 2.5% Expected Shortfall forecast in the test period,  GJR-Asyt).
 ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black",linewidth=0.5) +
  geom_line(aes(y = VaR_025_GJR), color = "seagreen4", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = ES_025_GJR), color = "seagreen4", linetype = "solid", linewidth = 0.7) +
  geom_point(
    data = subset(df_var, Breach_GJR == 1),
    aes(y = Returns),
    color = "green4", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  
  labs(
    x = "Date", y = "Log return"
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))


# Extract forecasted sigma (test)
garch_sd_test <- as.numeric(sigma(garch_forecast))
gjr_sd_test   <- as.numeric(sigma(gjr_forecast))

# Test returns
r_test <- as.numeric(sp500_test_returns)

# Compute standardized residuals manually
z_garch_test <- r_test / garch_sd_test
z_gjr_test   <- r_test / gjr_sd_test

pit_garch <- pasyt(z_garch_test, d1 = d1_garch_asyt, d2 = d2_garch_asyt)
pit_gjr   <- pasyt(z_gjr_test,   d1 = d1_gjr_asyt,   d2 = d2_gjr_asyt)

# Figure 4.7: Histogram of PIT values, GARCH- Asyt.
ggplot(data.frame(PIT = pit_garch), aes(x = PIT)) +
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

#Figure 4.7: Histogram of PIT values, GJR- Asyt.
ggplot(data.frame(PIT = pit_gjr), aes(x = PIT)) +
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

# Apply inverse standard normal (probit) to PIT values
z_pit_garch <- qnorm(pit_garch)
z_pit_gjr   <- qnorm(pit_gjr)


# Kolmogorov-Smirnov tests on transformed pit values
ks.test(z_pit_garch, "pnorm")
ks.test(z_pit_gjr, "pnorm")

# For GARCH(1,1)
BerkowitzTest(z_pit_garch)

# For GJR-GARCH(1,1)
BerkowitzTest(z_pit_gjr)


# Diebold-Mariano Test
loss_garch <- -log(fasyt(z_garch_test, d1_garch_asyt, d2_garch_asyt))
loss_gjr   <- -log(fasyt(z_gjr_test,   d1_gjr_asyt,   d2_gjr_asyt))
dm.test(
  e1 = loss_garch,
  e2 = loss_gjr,
  alternative = "two.sided",
  h = 1, power = 2
)



df_garch_ks <- get_ks_plot_data(z_pit_garch)
df_gjr_ks   <- get_ks_plot_data(z_pit_gjr)

# Figure 4.9: Kolmogorov-Smirnov test on ˆzt+1, GARCH-Asyt.
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
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18),
        legend.position = "right",
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 14))

# Figure 4.10: Kolmogorov-Smirnov test on ˆzt+1, GJR- Asyt.
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
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18),
        legend.position = "right",
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 14))


cc_garch <- christoffersen_test(df_var$Breach_GARCH, alpha = 0.025)
cc_gjr   <- christoffersen_test(df_var$Breach_GJR,   alpha = 0.025)

print(cc_garch)
print(cc_gjr)

## EVT approach evaluation (VaR) ##

# Total obs in training set
T_train <- length(y_garch_vec)

# Number of threshold exceedances
N_garch_exceed <- length(excess_garch)
N_gjr_exceed   <- length(excess_gjr)

# Confidence level
q <- 0.975
u <- threshold  # = 1.8

# EVT-based VaR from GPD (in standardized residual space)
VaR_gpd_garch_std <- u - (beta_garch / xi_garch) * (1 - ((T_train / N_garch_exceed) * (1 - q))^(-xi_garch))
VaR_gpd_gjr_std   <- u - (beta_gjr / xi_gjr) * (1 - ((T_train / N_gjr_exceed) * (1 - q))^(-xi_gjr))

# Transform to return space
VaR_gpd_garch <- -garch_sd * VaR_gpd_garch_std
VaR_gpd_gjr   <- -gjr_sd * VaR_gpd_gjr_std

df_var$VaR_EVT_GARCH <- as.numeric(VaR_gpd_garch)
df_var$VaR_EVT_GJR   <- as.numeric(VaR_gpd_gjr)

df_var$Breach_EVT_GARCH <- ifelse(df_var$Returns < df_var$VaR_EVT_GARCH, 1, 0)
df_var$Breach_EVT_GJR   <- ifelse(df_var$Returns < df_var$VaR_EVT_GJR, 1, 0)

# GPD Expected Shortfall (standardized residual space) 
ES_gpd_garch_std <- VaR_gpd_garch_std / (1 - xi_garch) +
  (beta_garch - xi_garch * u) / (1 - xi_garch)

ES_gpd_gjr_std <- VaR_gpd_gjr_std / (1 - xi_gjr) +
  (beta_gjr - xi_gjr * u) / (1 - xi_gjr)

# Transform to return space
ES_gpd_garch <- -garch_sd * ES_gpd_garch_std
ES_gpd_gjr   <- -gjr_sd   * ES_gpd_gjr_std

df_var$ES_EVT_GARCH <- as.numeric(ES_gpd_garch)
df_var$ES_EVT_GJR   <- as.numeric(ES_gpd_gjr)

# Figure 4.13
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black") +
  geom_line(aes(y = VaR_EVT_GARCH), color = "red4", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = subset(df_var, Breach_EVT_GARCH == 1),
    aes(y = Returns), color = "red", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  labs(
       x = "Date", y = "Log Return") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.13
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black") +
  geom_line(aes(y = VaR_EVT_GJR), color = "purple4", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = subset(df_var, Breach_EVT_GJR == 1),
    aes(y = Returns), color = "purple", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  labs(
       x = "Date", y = "Log Return") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.11: 2.5% VaR ES, and VaR breaches in the test period, GARCH-EVT
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black") +
  geom_line(aes(y = VaR_EVT_GARCH), color = "red4", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = ES_EVT_GARCH), color = "red4", linetype = "solid", linewidth = 1) +
  geom_point(
    data = subset(df_var, Breach_EVT_GARCH == 1),
    aes(y = Returns),
    color = "red", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  
  labs(
    x = "Date", y = "Log Return"
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 4.12: 2.5% VaR ES, and VaR breaches in the test period, GJR-EVT
ggplot(df_var, aes(x = Date)) +
  geom_line(aes(y = Returns), color = "black") +
  geom_line(aes(y = VaR_EVT_GJR), color = "purple4", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = ES_EVT_GJR), color = "purple4", linetype = "solid", linewidth = 1) +
  geom_point(
    data = subset(df_var, Breach_EVT_GJR == 1),
    aes(y = Returns),
    color = "purple", fill = "white", size = 2, shape = 21, stroke = 1
  ) +
  
  labs(
    x = "Date", y = "Log Return"
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

cc_evt_garch <- christoffersen_test(df_var$Breach_EVT_GARCH, alpha = 0.025)
cc_evt_gjr   <- christoffersen_test(df_var$Breach_EVT_GJR,   alpha = 0.025)

print(cc_evt_garch)
print(cc_evt_gjr)


## Backtest 2008 financial crisis ##

# Get 2008 returns
getSymbols("^GSPC", from = "2008-01-01", to = "2008-12-31")
sp500_2008 <- na.omit(Cl(GSPC))
sp500_2008_log_returns <- diff(log(sp500_2008))[-1]

garch_spec_fixed <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm",
  fixed.pars = as.list(coef(garch_fit))
)


gjr_spec_fixed <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm",
  fixed.pars = as.list(coef(gjr_fit))
)

# Forecast for 2008
garch_filtered_2008 <- ugarchfilter(spec = garch_spec_fixed, data = sp500_2008_log_returns)
gjr_filtered_2008 <- ugarchfilter(spec = gjr_spec_fixed, data = sp500_2008_log_returns)

# Extract conditional std deviation forecasts
garch_sd_2008 <- as.numeric(sigma(garch_filtered_2008))
gjr_sd_2008   <- as.numeric(sigma(gjr_filtered_2008))

# Asymmetric-t
VaR_garch_asyt <- garch_sd_2008 * qasyt(0.025, d1_garch_asyt, d2_garch_asyt)
VaR_gjr_asyt   <- gjr_sd_2008   * qasyt(0.025, d1_gjr_asyt, d2_gjr_asyt)

# EVT-GPD
VaR_std_evt_garch <- u - (beta_garch / xi_garch) * (1 - ((length(y_garch_vec) / length(excess_garch)) * 0.025)^(-xi_garch))
VaR_std_evt_gjr   <- u - (beta_gjr   / xi_gjr)   * (1 - ((length(y_gjr_vec)   / length(excess_gjr))   * 0.025)^(-xi_gjr))

VaR_garch_evt <- -garch_sd_2008 * VaR_std_evt_garch
VaR_gjr_evt   <- -gjr_sd_2008   * VaR_std_evt_gjr

df_2008 <- data.frame(
  Date           = index(sp500_2008_log_returns),
  Returns        = as.numeric(sp500_2008_log_returns),
  VaR_GARCH_ASYT = VaR_garch_asyt,
  VaR_GJR_ASYT   = VaR_gjr_asyt,
  VaR_GARCH_EVT  = VaR_garch_evt,
  VaR_GJR_EVT    = VaR_gjr_evt
)

df_2008$Breach_GARCH_ASYT <- df_2008$Returns < df_2008$VaR_GARCH_ASYT
df_2008$Breach_GJR_ASYT   <- df_2008$Returns < df_2008$VaR_GJR_ASYT
df_2008$Breach_GARCH_EVT  <- df_2008$Returns < df_2008$VaR_GARCH_EVT
df_2008$Breach_GJR_EVT    <- df_2008$Returns < df_2008$VaR_GJR_EVT

# Count
colSums(df_2008[, grep("Breach", colnames(df_2008))])



p1 <- plot_var(df_2008, "VaR_GARCH_ASYT", "Breach_GARCH_ASYT", "GARCH(1,1) + Asymmetric-t", "deepskyblue3")
p2 <- plot_var(df_2008, "VaR_GJR_ASYT",   "Breach_GJR_ASYT",   "GJR-GARCH(1,1) + Asymmetric-t", "seagreen4")
p3 <- plot_var(df_2008, "VaR_GARCH_EVT",  "Breach_GARCH_EVT",  "GARCH(1,1) + EVT (GPD)", "red4")
p4 <- plot_var(df_2008, "VaR_GJR_EVT",    "Breach_GJR_EVT",    "GJR-GARCH(1,1) + EVT (GPD)", "purple4")


#Figure 4.13: Stress test of the models
grid.arrange(p1, p2, p3, p4, ncol = 2)
