## Purpose: Replicates Chapter 1 empirical return characteristics analysis.
## Inputs: Yahoo S&P 500 data and local bond data.
## Dependencies: R/00_setup.R (packages, seed, and shared paths).

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
  geom_col(fill = "deepskyblue3", color = "black", width = 0.7) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_hline(aes(yintercept = Upper), linetype = "dashed", color = "gray40") +
  geom_hline(aes(yintercept = Lower), linetype = "dashed", color = "gray40") +
  labs(
    x = "Lag order", 
    y = "Autocorrelation of daily log returns"
  ) +
  theme_tufte(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    text = element_text(size = 18)
  )


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
  geom_col(fill = "deepskyblue3", color = "black", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  labs(
    x = "Lag order", 
    y = "Autocorrelation of squared returns"
  ) +
  theme_tufte(base_family = "Times") + 
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.text = element_text(size = 20),    
    axis.title = element_text(size = 22),
    text = element_text(size = 18)
  )

# Load Bond Data
file_path <- bond_data_path
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
