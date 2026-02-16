# Goodness-of-fit diagnostic plotting helpers.

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
