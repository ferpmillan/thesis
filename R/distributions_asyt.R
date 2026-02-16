# Asymmetric Student's t distribution and diagnostics helpers.
# Source: Christoffersen (2012) implementation used in the thesis.

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
