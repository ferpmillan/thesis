# Extreme value theory helper functions used in Chapter 3 and Chapter 4.

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

# NOTE: `title` is kept optional for backward compatibility with original call sites.
plot_gpd_cdf <- function(excess, xi, beta, title = NULL) {
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
