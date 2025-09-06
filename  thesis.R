##### THESIS. 
##### Fernando Pérez Millán :)

#### Libraries ####
library(readxl)
library(quantmod)
library(ggplot2)
library(ggthemes)
library(forecast)
library(MASS)
library(rugarch)
library(rmgarch)
library(reshape2)
library(fGarch)
library(evir)
library(xts)
library(tseries)
library(gridExtra)
library(zoo)
library(copula)
library(matrixcalc)
library(pracma)
library(MCMCpack)
library(viridis)


##### Chapter 1. Empirical Characteristics of Returns. ##### 

set.seed(191023)

# Define the date range (training set)
start_date <- "2014-01-01"
end_date <- "2025-01-01"

# Get S&P 500 data from Yahoo Finance
getSymbols("^GSPC", src = "yahoo", from = start_date, to = end_date)

# Extract the closing prices
sp500_close <- Cl(GSPC)

# Calculate daily log returns
sp500_log_returns <- diff(log(sp500_close))
sp500_log_returns <- na.omit(sp500_log_returns)

# Calculate autocorrelation for lags 1 to 100
acf_values <- Acf(sp500_log_returns, lag.max = 100, plot = FALSE)

# Convert to data frame and remove lag 0
acf_df <- data.frame(Lag = acf_values$lag[-1], ACF = acf_values$acf[-1])

# 
n_returns_train <- length(sp500_log_returns)
conf_bound <- 1.96/sqrt(n_returns_train)
acf_df$Upper <- conf_bound
acf_df$Lower <- -conf_bound

# Figure 1.1: Autocorrelation of daily S&P 500 log returns from January 1, 2014, to December 31, 2024.
ggplot(acf_df, aes(x = Lag, y = ACF)) +
  geom_line(stat = "identity", color = "deepskyblue3", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  geom_hline(aes(yintercept = Upper), linetype = "dashed", color = "gray40") +
  geom_hline(aes(yintercept = Lower), linetype = "dashed", color = "gray40") +
  labs(
    x = "Lag order", 
    y = "Autocorrelation of daily log returns") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),  
        text = element_text(size = 18)) 


# Fit a normal distribution using Maximum Likelihood Estimation (MLE)
log_returns_vector <- as.numeric(sp500_log_returns)
mle_fit <- fitdistr(log_returns_vector, "normal")
mu_hat <- mle_fit$estimate["mean"]
sigma_hat <- mle_fit$estimate["sd"]

# Figure 1.2: Histogram of daily S&P 500 log returns and the normal distribution January 1, 2014 - December 31, 2024.
ggplot(data.frame(Returns = log_returns_vector), aes(x = Returns)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "deepskyblue3", color = "black", alpha =1) +
  stat_function(fun = dnorm, args = list(mean = mu_hat, sd = sigma_hat), 
                color = "red", linetype = "solid", linewidth = 1.2) +
  labs(
    x = "Log Return",
    y = "Density") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 26),    
        axis.title = element_text(size = 28),
        text = element_text(size = 24))

# Squared returns
sp500_squared_returns <- sp500_log_returns^2

# Calculate autocorrelation for squared returns (lag 1 to 100)
acf_values_squared <- Acf(sp500_squared_returns, lag.max = 100, plot = FALSE)
acf_df_squared <- data.frame(Lag = acf_values_squared $lag[-1], ACF = acf_values_squared $acf[-1])

# Figure 1.3: Autocorrelation of squared daily S&P 500 log returns from January 1, 2014, to December 31, 2024.
ggplot(acf_df_squared, aes(x = Lag, y = ACF)) +
  geom_line(stat = "identity", color = "deepskyblue3",linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  labs( 
       x = "Lag order", 
       y = "Autocorrelation of squared returns") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Load Bond Data
file_path <- "/Users/fernandoperezmillan/Desktop/final_bond_data.xlsx"
bond_data <- read_excel(file_path)

# Convert Date column to Date type
bond_data$Dates <- as.Date(bond_data$Dates)

# Convert to xts
bond_xts <- xts(bond_data$PX_LAST, order.by = bond_data$Dates)

# Calculate log returns
bond_log_returns <- diff(log(bond_xts))
bond_log_returns <- na.omit(bond_log_returns)

# Convert to data frame if needed
final_bond_df <- data.frame(
  Date = index(bond_log_returns),
  Log_Return = as.numeric(bond_log_returns)
)

# Merge with S&P 500 log returns (assuming sp500_log_returns already exists)
merged_returns <- merge(
  sp500_log_returns,
  bond_log_returns,
  join = "inner"
)

# 7. Rename columns for clarity
colnames(merged_returns) <- c("SP500", "Bond")

# Compute rolling 25-day covariance
rolling_corr <- rollapply(
  data = merged_returns,
  width = 25, 
  FUN = function(x) cov(x[,1], x[,2], use = "complete.obs"),
  by.column = FALSE,
  align = "right"
)

# Convert to data frame for plotting
rolling_cov_df <- data.frame(
  Date = index(rolling_corr),
  Correlation = coredata(rolling_corr)
)

# Figure 1.4: Rolling covariance between S&P 500 and 10-year treasury note index.
ggplot(rolling_cov_df, aes(x = Date, y = Correlation)) +
  geom_line(color = "deepskyblue3", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  labs(
       x = "Date",
       y = "Covariance") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

##### Chapter 2. Volatility Modeling. #####

# Convert to data frame for plotting
squared_returns_df <- data.frame(
  Date = index(sp500_squared_returns),
  Squared_Returns = coredata(sp500_squared_returns)
)

# Figure 2.1: Squared S&P 500 returns, 2014-2024.
ggplot(squared_returns_df, aes(x = Date, y = GSPC.Close)) +
  geom_line(color = "deepskyblue3", linewidth = 0.5) +
  labs(
       x = "Date",
       y = "Squared returns") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),
        text = element_text(size = 18))


# Define the number of lags for Ljun-Box test as m = ln(T)
T <- length(sp500_squared_returns)
m <- floor(log(T))  # Take the natural log and round down

# Perform Ljung-Box test
ljung_box_result <- Box.test(sp500_squared_returns, lag = m, type = "Ljung-Box")
print(ljung_box_result)

# Specify GARCH(1,1) Model
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "norm" #QMLE
)

# Fit GARCH(1,1) Model
garch_fit <- ugarchfit(spec = garch_spec, data = sp500_log_returns)

# Specify GJR-GARCH(1,1) Model
gjr_spec <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "norm"
)

# Fit GJR-GARCH(1,1) Model
gjr_fit <- ugarchfit(spec = gjr_spec, data = sp500_log_returns)

# Extract conditional variance (forecasted volatility)
garch_vol <- sigma(garch_fit)^2  # GARCH(1,1) variance
gjr_vol <- sigma(gjr_fit)^2  # GJR-GARCH(1,1) variance

# Convert to data frame for plotting
garch_df <- data.frame(
  Date = index(sp500_squared_returns),
  Squared_Returns = coredata(sp500_squared_returns),
  GARCH_Variance = coredata(garch_vol)
)

gjr_df <- data.frame(
  Date = index(sp500_squared_returns),
  Squared_Returns = coredata(sp500_squared_returns),
  GJR_GARCH_Variance = coredata(gjr_vol)
)

# Figure 2.2: Squared S&P 500 returns with GARCH(1,1) variance parameters estimated using QMLE.
ggplot(garch_df, aes(x = Date)) +
  geom_line(aes(y = GSPC.Close, color = "Squared returns"), linewidth = 0.5) +
  geom_line(aes(y = GARCH_Variance, color = "GARCH(1,1)"), linewidth = 0.5) +
  scale_color_manual(values = c("red", "deepskyblue3")) +
  labs(
       x = "Date",
       y = "Variance") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),
        text = element_text(size = 18)) +
  guides(color = guide_legend(""))

# Figure 2.3: Squared S&P 500 returns with GJR-GARCH(1,1) variance parameters estimated using QMLE.
ggplot(gjr_df, aes(x = Date)) +
  geom_line(aes(y = GSPC.Close, color = "Squared returns"), linewidth = 0.5) +
  geom_line(aes(y = GJR_GARCH_Variance, color = "GJR-GARCH(1,1)"), linewidth = 0.5) +
  scale_color_manual(values = c("red", "deepskyblue3")) +
  labs(
       x = "Date",
       y = "Variance") +
  theme_tufte(base_family = "Times") + 
  theme(panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 20),    
        axis.title = element_text(size = 22),
        text = element_text(size = 18)) +
  guides(color = guide_legend(title = ""))


# Extract GARCH(1,1) parameters
garch_params <- coef(garch_fit)
print("GARCH(1,1) Parameters:")
print(garch_params)

# Extract GJR-GARCH(1,1) parameters
gjr_params <- coef(gjr_fit)
print("GJR-GARCH(1,1) Parameters:")
print(gjr_params)

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

# Density function (Christoffersen, 2012)
fasyt <- function(z, d1, d2) {
  # Constants
  C <- gamma((d1 + 1) / 2) / (gamma(d1 / 2) * sqrt(pi * (d1 - 2)))
  A <- 4 * d2 * C * (d1 - 2) / (d1 - 1)
  B <- sqrt(1 + 3 * d2^2 - A^2)
  
  # Piecewise PDF
  left  <- z < -A / B
  right <- !left
  
  dens <- numeric(length(z))
  denom_left  <- (1 - d2)^2 * (d1 - 2)
  denom_right <- (1 + d2)^2 * (d1 - 2)
  
  # Density formulas
  dens[left] <- B * C * (1 + ((B * z[left] + A)^2) / denom_left)^(-(1 + d1) / 2)
  dens[right] <- B * C * (1 + ((B * z[right] + A)^2) / denom_right)^(-(1 + d1) / 2)
  
  return(dens)
}

# Loglikelihood to optimize
loglik_asyt <- function(params, z) {
  d1 <- params[1]
  d2 <- params[2]
  
  # Constraints to ensure valid pdf
  if (d1 <= 2 || d2 <= 0) return(1e6)
  
  fz <- fasyt(z, d1, d2)
  
  # Avoid log(0) / invalid
  if (any(fz <= 0 | is.nan(fz))) return(1e6)
  
  return(-sum(log(fz)))
}

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

# PDF for plotting
plot_fitted_asyt <- function(z_vec, d1, d2) {
  x_grid <- seq(min(z_vec), max(z_vec), length.out = 1000)
  fitted_density <- fasyt(x_grid, d1, d2)
  
  ggplot(data.frame(z = z_vec), aes(x = z)) +
    geom_histogram(aes(y = ..density..), bins = 50, fill = "deepskyblue3", color = "black") +
    geom_line(data = data.frame(x = x_grid, y = fitted_density),
              aes(x = x, y = y), color = "red", linewidth = 1.2) +
    labs(
      x = "Standardized return",
      y = "Density") +
    theme_tufte(base_family = "Times") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          text = element_text(size = 18))
}

# Figure 3.5: Histogram of standardized returns with the fitted asymmetric t-distribution, GARCH(1,1).
plot_fitted_asyt(z_garch_vec, d1_garch_asyt, d2_garch_asyt)

# Figure 3.6: Histogram of standardized returns with the fitted asymmetric t-distribution, GJR-GARCH(1,1).
plot_fitted_asyt(z_gjr_vec, d1_gjr_asyt, d2_gjr_asyt)

# Quantiles of the asymmetric t-distribution (Christoffersen, 2012)
qasyt <- function(p, d1, d2) {
  # Constants
  C <- gamma((d1 + 1) / 2) / (gamma(d1 / 2) * sqrt(pi * (d1 - 2)))
  A <- 4 * d2 * C * (d1 - 2) / (d1 - 1)
  B <- sqrt(1 + 3 * d2^2 - A^2)
  
  q <- numeric(length(p))
  cutoff <- (1 - d2) / 2
  
  # Left tail
  idx_left <- p < cutoff
  t_quant_left <- qt(p[idx_left] / (1 - d2), df = d1)
  q[idx_left] <- (1 / B) * ((1 - d2) * sqrt((d1 - 2) / d1) * t_quant_left - A)
  
  # Right tail
  idx_right <- !idx_left
  t_quant_right <- qt((p[idx_right] + d2) / (1 + d2), df = d1)
  q[idx_right] <- (1 / B) * ((1 + d2) * sqrt((d1 - 2) / d1) * t_quant_right - A)
  
  return(q)
}

plot_qq_asyt <- function(z_vec, d1, d2) {
  n <- length(z_vec)
  p <- ppoints(n)
  theoretical <- qasyt(p, d1, d2)
  sample <- sort(z_vec)
  
  df <- data.frame(Theoretical = theoretical, Sample = sample)
  
  ggplot(df, aes(x = Theoretical, y = Sample)) +
    geom_point(shape = 21, fill = "deepskyblue3", color = "black", size = 2, stroke = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
    labs(
      x = "Asymmetric Student's t quantile",
      y = "Empirical quantile") +
    theme_tufte(base_family = "Times") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          text = element_text(size = 18))
}

# Figure 3.7: QQ plot of standardized returns against the asymmetric t-distribution, GARCH(1,1).
plot_qq_asyt(z_garch_vec, d1_garch_asyt, d2_garch_asyt)

# Figure 3.8: QQ plot of standardized returns against the asymmetric t-distribution, GJR- GARCH(1,1).
plot_qq_asyt(z_gjr_vec, d1_gjr_asyt, d2_gjr_asyt)

# CDF (Christoffersen, 2012)
pasyt <- function(z, d1, d2) {
  # Constants
  C <- gamma((d1 + 1) / 2) / (gamma(d1 / 2) * sqrt(pi * (d1 - 2)))
  A <- 4 * d2 * C * (d1 - 2) / (d1 - 1)
  B <- sqrt(1 + 3 * d2^2 - A^2)
  
  cutoff <- -A / B
  p <- numeric(length(z))
  
  # Left region
  idx_left <- z < cutoff
  arg_left <- ((B * z[idx_left] + A) / ((1 - d2) * sqrt((d1 - 2) / d1)))
  p[idx_left] <- (1 - d2) * pt(arg_left, df = d1)
  
  # Right region
  idx_right <- !idx_left
  arg_right <- ((B * z[idx_right] + A) / ((1 + d2) * sqrt((d1 - 2) / d1)))
  p[idx_right] <- -d2 + (1 + d2) * pt(arg_right, df = d1)
  
  return(p)
}

plot_cdf_asyt <- function(z_vec, d1, d2) {
  z_sorted <- sort(z_vec)
  emp_cdf <- ecdf(z_vec)(z_sorted)
  theo_cdf <- pasyt(z_sorted, d1, d2)
  
  df <- data.frame(z = z_sorted,
                   Empirical = emp_cdf,
                   Theoretical = theo_cdf)
  
  ggplot(df, aes(x = z)) +
    geom_line(aes(y = Empirical, color = "Empirical CDF"), linewidth = 1) +
    geom_line(aes(y = Theoretical, color = "Theoretical CDF"), linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("deepskyblue3", "red")) +
    labs(
      x ="Standardized return",
      y = "CDF",
      color = "") +
    theme_tufte(base_family = "Times") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          text = element_text(size = 18),
          legend.position = "right")
}

# Figure 3.9: Emprical CDF and theoretical asymmetric Student’s t CDF of standardized returns, GARCH(1,1).
plot_cdf_asyt(z_garch_vec, d1_garch_asyt, d2_garch_asyt)

# Figure 3.10: Emprical CDF and theoretical asymmetric Student’s t CDF of standardized returns, GJR-GARCH(1,1).
plot_cdf_asyt(z_gjr_vec, d1_gjr_asyt, d2_gjr_asyt)


## Generalized Pareto Distribution ##

# Standardized negative residuals
y_garch_vec <- as.numeric(-z_garch)
y_gjr_vec   <- as.numeric(-z_gjr)

# Obtain the mean excess function
get_mef_data <- function(x) {
  x <- sort(x)
  n <- length(x)
  u <- x[1:(n - 1)]
  e <- sapply(1:(n - 1), function(i) {
    excess <- x[(i + 1):n] - x[i]
    mean(excess)
  })
  count <- sapply(1:(n - 1), function(i) n - i)
  data.frame(u = u, e = e, n = count)
}

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
  labs(,
       x = "EVT quantile",
       y = "Empirical quantile") +
  theme_tufte(base_family = "Times") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        text = element_text(size = 18))

# Plot function for empirical vs theoretical CDF of GPD
plot_gpd_cdf <- function(excess, xi, beta, title) {
  x_sorted <- sort(excess)
  ecdf_vals <- ecdf(excess)(x_sorted)
  theoretical_vals <- pgpd(x_sorted, xi = xi, beta = beta)
  
  df <- data.frame(
    x = x_sorted,
    Empirical = ecdf_vals,
    Theoretical = theoretical_vals
  )
  
  ggplot(df, aes(x = x)) +
    geom_line(aes(y = Empirical, color = "Empirical CDF"), linewidth = 1) +
    geom_line(aes(y = Theoretical, color = "Theoretical CDF"), linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("deepskyblue3", "red")) +
    labs(
         x = "Excess over threshold",
         y = "CDF",
         color = "") +
    theme_tufte(base_family = "Times") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          text = element_text(size = 18),
          legend.position = "right")
}

# Figure 3.15: Empirical CDF and theoretical EVT CDF of standardized returns, GARCH(1,1).
plot_gpd_cdf(
  excess = excess_garch,
  xi = xi_garch,
  beta = beta_garch,
)

# Figure 3.16: Empirical CDF and theoretical EVT CDF of standardized returns, GJR-GARCH(1,1).
plot_gpd_cdf(
  excess = excess_gjr,
  xi = xi_gjr,
  beta = beta_gjr,
)


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

qlike <- function(z2, eps = 1e-8) {
  z2_adj <- pmax(z2, eps)
  z2_adj - log(z2_adj) - 1
}

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

es_asyt <- function(p, d1, d2) {
  A <- 4 * d2 * (gamma((d1 + 1) / 2) / (gamma(d1 / 2) * sqrt(pi * (d1 - 2)))) * ((d1 - 2) / d1)
  B <- sqrt(1 + 3 * d2^2 - A^2)
  C <- gamma((d1 + 1) / 2) / (gamma(d1 / 2) * sqrt(pi * (d1 - 2)))
  
  Q <- qasyt(p, d1, d2)
  
  term1 <- (C * (1 - d2)^2) / (B * p)
  term2 <- (1 + (1 / (d1 - 2)) * ((B * Q + A) / (1 - d2))^2) ^ ((1 - d1) / 2)
  
  term3 <- (A * C * (1 - d2) * sqrt(pi * (d1 - 2)) * gamma(d1 / 2)) /
    (B * p * gamma((d1 + 1) / 2))
  
  td1_val <- dt(
    sqrt(d1 / (d1 - 2)) * (B * Q + A) / (1 - d2),
    df = d1
  )
  
  ES <- term1 * term2 - term3 * td1_val
  return(ES)
}

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


get_ks_plot_data <- function(z_vec) {
  ecdf_func <- ecdf(z_vec)
  z_sorted <- sort(z_vec)
  ecdf_vals <- ecdf_func(z_sorted)
  theo_cdf <- pnorm(z_sorted)
  ks_dist <- abs(ecdf_vals - theo_cdf)
  max_idx <- which.max(ks_dist)
  
  data.frame(
    z = z_sorted,
    ECDF = ecdf_vals,
    NormalCDF = theo_cdf,
    KS = ks_dist,
    max_z = z_sorted[max_idx],
    max_ecdf = ecdf_vals[max_idx],
    max_norm = theo_cdf[max_idx]
  )
}

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

christoffersen_test <- function(breach_vec, alpha = 0.025, B = 999) {
  n <- length(breach_vec)
  x <- sum(breach_vec)
  pi_hat <- x / n
  
  # UC: Theoretical 
  L0 <- ((1 - alpha)^(n - x)) * (alpha^x)
  L1 <- ((1 - pi_hat)^(n - x)) * (pi_hat^x)
  LR_uc <- -2 * log(L0 / L1)
  p_uc <- 1 - pchisq(LR_uc, df = 1)
  
  # IND: Theoretical
  v <- breach_vec
  n00 <- sum(v[-n] == 0 & v[-1] == 0)
  n01 <- sum(v[-n] == 0 & v[-1] == 1)
  n10 <- sum(v[-n] == 1 & v[-1] == 0)
  n11 <- sum(v[-n] == 1 & v[-1] == 1)
  
  pi01 <- n01 / (n00 + n01)
  pi11 <- n11 / (n10 + n11)
  pi_hat2 <- (n01 + n11) / (n00 + n01 + n10 + n11)
  
  L_indep <- ((1 - pi01)^n00) * (pi01^n01) * ((1 - pi11)^n10) * (pi11^n11)
  L_dep <- ((1 - pi_hat2)^(n00 + n10)) * (pi_hat2^(n01 + n11))
  LR_ind <- -2 * log(L_dep / L_indep)
  p_ind <- 1 - pchisq(LR_ind, df = 1)
  
  # CC
  LR_cc <- LR_uc + LR_ind
  p_cc <- 1 - pchisq(LR_cc, df = 2)
  
  # Simulation for empirical p-values
  sim_LR_uc <- numeric(B)
  sim_LR_ind <- numeric(B)
  sim_LR_cc <- numeric(B)
  
  for (b in 1:B) {
    sim <- rbinom(n, 1, alpha)
    
    # UC
    xb <- sum(sim)
    pib <- xb / n
    L0b <- ((1 - alpha)^(n - xb)) * (alpha^xb)
    L1b <- ((1 - pib)^(n - xb)) * (pib^xb)
    sim_LR_uc[b] <- -2 * log(L0b / L1b)
    
    # IND
    n00b <- sum(sim[-n] == 0 & sim[-1] == 0)
    n01b <- sum(sim[-n] == 0 & sim[-1] == 1)
    n10b <- sum(sim[-n] == 1 & sim[-1] == 0)
    n11b <- sum(sim[-n] == 1 & sim[-1] == 1)
    
    pi01b <- ifelse(n00b + n01b > 0, n01b / (n00b + n01b), 0)
    pi11b <- ifelse(n10b + n11b > 0, n11b / (n10b + n11b), 0)
    pi_hat2b <- ifelse((n00b + n01b + n10b + n11b) > 0,
                       (n01b + n11b) / (n00b + n01b + n10b + n11b), 0)
    
    L_indepb <- ((1 - pi01b)^n00b) * (pi01b^n01b) *
      ((1 - pi11b)^n10b) * (pi11b^n11b)
    L_depb <- ((1 - pi_hat2b)^(n00b + n10b)) * (pi_hat2b^(n01b + n11b))
    sim_LR_ind[b] <- -2 * log(L_depb / L_indepb)
    
    sim_LR_cc[b] <- sim_LR_uc[b] + sim_LR_ind[b]
  }
  
  # Empirical (simulated) p-values
  p_uc_emp  <- (1 + sum(sim_LR_uc  > LR_uc))  / (B + 1)
  p_ind_emp <- (1 + sum(sim_LR_ind > LR_ind)) / (B + 1)
  p_cc_emp  <- (1 + sum(sim_LR_cc  > LR_cc))  / (B + 1)
  
  return(list(
    Total_Obs = n,
    Violations = x,
    
    # Theoretical
    LR_uc = round(LR_uc, 4), p_uc_theo = round(p_uc, 4),
    LR_ind = round(LR_ind, 4), p_ind_theo = round(p_ind, 4),
    LR_cc = round(LR_cc, 4), p_cc_theo = round(p_cc, 4),
    
    # Empirical
    p_uc_emp  = round(p_uc_emp, 4),
    p_ind_emp = round(p_ind_emp, 4),
    p_cc_emp  = round(p_cc_emp, 4)
  ))
}

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


plot_var <- function(df, var_col, breach_col, title, color) {
  ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Returns), color = "black") +
    geom_line(aes_string(y = var_col), color = color, linetype = "dashed", linewidth = 1) +
    geom_point(
      data = subset(df, df[[breach_col]]),
      aes(y = Returns),
      shape = 21, fill = "white", color = color, stroke = 1.2, size = 2.5
    ) +
    scale_x_date(date_labels = "%b") +  # show only month abbrev (Jan, Feb, ...)
    labs(x = "Month", y = "Log Return") +
    theme_tufte(base_family = "Times") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          text = element_text(size = 18))
}

p1 <- plot_var(df_2008, "VaR_GARCH_ASYT", "Breach_GARCH_ASYT", "GARCH(1,1) + Asymmetric-t", "deepskyblue3")
p2 <- plot_var(df_2008, "VaR_GJR_ASYT",   "Breach_GJR_ASYT",   "GJR-GARCH(1,1) + Asymmetric-t", "seagreen4")
p3 <- plot_var(df_2008, "VaR_GARCH_EVT",  "Breach_GARCH_EVT",  "GARCH(1,1) + EVT (GPD)", "red4")
p4 <- plot_var(df_2008, "VaR_GJR_EVT",    "Breach_GJR_EVT",    "GJR-GARCH(1,1) + EVT (GPD)", "purple4")


#Figure 4.13: Stress test of the models
grid.arrange(p1, p2, p3, p4, ncol = 2)

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

u_sim <- rCopula(1000000, copula_fit)

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

u_sim_gjr <- rCopula(1000000, copula_fit_gjr)

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

# 3. Function to prepare data for plotting
get_ks_plot_data <- function(z_vec) {
  ecdf_func <- ecdf(z_vec)
  z_sorted <- sort(z_vec)
  ecdf_vals <- ecdf_func(z_sorted)
  theo_cdf <- pnorm(z_sorted)
  ks_dist <- abs(ecdf_vals - theo_cdf)
  max_idx <- which.max(ks_dist)
  
  data.frame(
    z = z_sorted,
    ECDF = ecdf_vals,
    NormalCDF = theo_cdf,
    KS = ks_dist,
    max_z = z_sorted[max_idx],
    max_ecdf = ecdf_vals[max_idx],
    max_norm = theo_cdf[max_idx]
  )
}

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