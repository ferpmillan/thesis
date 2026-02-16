# Forecast evaluation and backtesting utility functions.

qlike <- function(z2, eps = 1e-8) {
  z2_adj <- pmax(z2, eps)
  z2_adj - log(z2_adj) - 1
}

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
  
  set.seed(191023)
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
  set.seed(191023)
}

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
