## Purpose: Replicates Chapter 3 shock distribution modeling and EVT analysis.
## Inputs: Fitted GARCH/GJR objects from Chapter 2.
## Dependencies: Helper functions sourced in scripts/00_run_all.R.

##### Chapter 3. Shock Modeling. #####

# Extract standardized residuals
z_garch <- residuals(garch_fit, standardize = TRUE)
z_gjr <- residuals(gjr_fit, standardize = TRUE)

# Convert to numeric vectors
z_garch_vec <- as.numeric(z_garch)
z_gjr_vec <- as.numeric(z_gjr)

# KS test against standard normal
ks_garch <- ks.test(z_garch_vec, "pnorm", mean = 0, sd = 1)
ks_gjr   <- ks.test(z_gjr_vec, "pnorm", mean = 0, sd = 1)

# Print results
print("KS Test - GARCH(1,1) Residuals:")
print(ks_garch)

print("KS Test - GJR-GARCH(1,1) Residuals:")
print(ks_gjr)

# Figure 3.1: Histogram of standardized returns, GARCH(1,1).
ggplot(data.frame(z = z_garch_vec), aes(x = z)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "deepskyblue3", color = "black", alpha = 1) +
  stat_function(fun = dnorm, 
                args = list(mean = 0, sd = 1),
                color = "red", linewidth = 1.2) +
  labs(
       x = "Standardized return",
       y = "Density") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 3.2: Histogram of standardized returns, GJR-GARCH(1,1).
ggplot(data.frame(z = z_gjr_vec), aes(x = z)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "deepskyblue3", color = "black", alpha = 1) +
  stat_function(fun = dnorm, 
                args = list(mean = 0, sd = 1),
                color = "red", linewidth = 1.2) +
  labs(
       x = "Standardized return",
       y = "Density") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 3.3: Normal QQ plot of daily standardized S&P 500 returns, GARCH(1,1).
ggplot(data.frame(z = z_garch_vec), aes(sample = z)) +
  stat_qq(shape = 21, fill = "deepskyblue3", color = "black", size = 2.5, stroke = 0.4) +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(
       x = "Normal quantile",
       y = "Standardized return quantile") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 3.4: Normal QQ plot of daily standardized S&P 500 returns, GJR-GARCH(1,1).
ggplot(data.frame(z = z_gjr_vec), aes(sample = z)) +
  stat_qq(shape = 21, fill = "deepskyblue3", color = "black",size = 2.5, stroke = 0.3) +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(
    x = "Normal quantile",
    y = "Standardized return quantile") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size =20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))


## Asymmetric Student's t-distribution ##
# The asymmetric t-distribution is not included in any package, we implement it manually here.

# Optimize Christoffersen asymmetric t on GARCH standardized residuals
fit_garch_asyt <- optim(par = c(6, 0.5), fn = loglik_asyt, z = z_garch_vec,
                        method = "L-BFGS-B", lower = c(2.01, -1.01), upper = c(100, .99))

fit_gjr_asyt <- optim(par = c(6, 0.5), fn = loglik_asyt, z = z_gjr_vec,
                      method = "L-BFGS-B", lower = c(2.01, -1.01), upper = c(100, .99))

# Extract estimates
d1_garch_asyt <- fit_garch_asyt$par[1]
d2_garch_asyt <- fit_garch_asyt$par[2]

d1_gjr_asyt <- fit_gjr_asyt$par[1]
d2_gjr_asyt <- fit_gjr_asyt$par[2]

cat("GARCH(1,1) fitted Asyt:\n")
cat(sprintf("  d1 = %.4f\n  d2 = %.4f\n", d1_garch_asyt, d2_garch_asyt))

cat("\nGJR-GARCH(1,1) fitted Asyt:\n")
cat(sprintf("  d1 = %.4f\n  d2 = %.4f\n", d1_gjr_asyt, d2_gjr_asyt))

# Figure 3.5: Histogram of standardized returns with the fitted asymmetric t-distribution, GARCH(1,1).
plot_fitted_asyt(z_garch_vec, d1_garch_asyt, d2_garch_asyt)

# Figure 3.6: Histogram of standardized returns with the fitted asymmetric t-distribution, GJR-GARCH(1,1).
plot_fitted_asyt(z_gjr_vec, d1_gjr_asyt, d2_gjr_asyt)

# Figure 3.7: QQ plot of standardized returns against the asymmetric t-distribution, GARCH(1,1).
plot_qq_asyt(z_garch_vec, d1_garch_asyt, d2_garch_asyt)

# Figure 3.8: QQ plot of standardized returns against the asymmetric t-distribution, GJR- GARCH(1,1).
plot_qq_asyt(z_gjr_vec, d1_gjr_asyt, d2_gjr_asyt)

# Figure 3.9: Emprical CDF and theoretical asymmetric Student’s t CDF of standardized returns, GARCH(1,1).
plot_cdf_asyt(z_garch_vec, d1_garch_asyt, d2_garch_asyt)

# Figure 3.10: Emprical CDF and theoretical asymmetric Student’s t CDF of standardized returns, GJR-GARCH(1,1).
plot_cdf_asyt(z_gjr_vec, d1_gjr_asyt, d2_gjr_asyt)


## Generalized Pareto Distribution ##

# Standardized negative residuals
y_garch_vec <- as.numeric(-z_garch)
y_gjr_vec   <- as.numeric(-z_gjr)

# Compute and filter
mef_garch <- get_mef_data(y_garch_vec)
mef_gjr   <- get_mef_data(y_gjr_vec)
mef_garch_clean <- subset(mef_garch, n >= 50)
mef_gjr_clean   <- subset(mef_gjr, n >= 50)

# Figure 3.11: Mean excess function e(u), GARCH(1,1).
ggplot(mef_garch_clean, aes(x = u, y = e)) +
  geom_point(shape = 21, fill = "deepskyblue3", color = "black", size = 2.5, stroke = 0.3) +
  labs(
    x = expression("Threshold " * italic(u)),
    y = expression("Mean excess " * italic(e)(italic(u)))
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 3.12: Mean excess function e(u), GJR-GARCH(1,1).
ggplot(mef_gjr_clean, aes(x = u, y = e)) +
  geom_point(shape = 21, fill = "deepskyblue3", color = "black", size = 2.5, stroke = 0.3) +
  labs(
    x = expression("Threshold " * italic(u)),
    y = expression("Mean excess " * italic(e)(italic(u)))
  ) +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))


# Set threshold
threshold <- 1.8

# Fit GPS
fit_gpd_garch <- gpd(y_garch_vec, threshold = threshold)
fit_gpd_gjr   <- gpd(y_gjr_vec, threshold = threshold)

# Extract exceedances (for histogram)
excess_garch <- y_garch_vec[y_garch_vec > threshold] - threshold
excess_gjr   <- y_gjr_vec[y_gjr_vec > threshold] - threshold

# GPD parameters
xi_garch   <- fit_gpd_garch$par.ests["xi"]
beta_garch <- fit_gpd_garch$par.ests["beta"]

xi_gjr   <- fit_gpd_gjr$par.ests["xi"]
beta_gjr <- fit_gpd_gjr$par.ests["beta"]

# Theoretical quantiles from GPD
q_theoretical_garch <- qgpd(ppoints(length(excess_garch)), xi = xi_garch, beta = beta_garch)
q_theoretical_gjr   <- qgpd(ppoints(length(excess_gjr)),   xi = xi_gjr,   beta = beta_gjr)

# Create dataframes
df_qq_garch <- data.frame(
  Theoretical = q_theoretical_garch,
  Sample = sort(excess_garch)
)

df_qq_gjr <- data.frame(
  Theoretical = q_theoretical_gjr,
  Sample = sort(excess_gjr)
)

# Figure 3.13: QQ plot of standardized returns against the generalized Pareto distribution, GARCH(1,1).
ggplot(df_qq_garch, aes(x = -Theoretical, y = -Sample)) +
  geom_point(shape = 21, fill = "deepskyblue3", color = "black", size = 2, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(
       x = "EVT quantile",
       y = "Empirical quantile") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

#Figure 3.14: QQ plot of standardized returns against the generalized Pareto distribution, GJR-GARCH(1,1).
ggplot(df_qq_gjr, aes(x = -Theoretical, y = -Sample)) +
  geom_point(shape = 21, fill = "deepskyblue3", color = "black", size = 2, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(
       x = "EVT quantile",
       y = "Empirical quantile") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Figure 3.15: Empirical CDF and theoretical EVT CDF of standardized returns, GARCH(1,1).
plot_gpd_cdf(
  excess = excess_garch,
  xi = xi_garch,
  beta = beta_garch
)

# Figure 3.16: Empirical CDF and theoretical EVT CDF of standardized returns, GJR-GARCH(1,1).
plot_gpd_cdf(
  excess = excess_gjr,
  xi = xi_gjr,
  beta = beta_gjr
)

