## Purpose: Replicates Chapter 2 conditional volatility modeling and diagnostics.
## Inputs: Chapter 1 returns objects.
## Dependencies: GARCH specs and returns produced in previous scripts.

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
