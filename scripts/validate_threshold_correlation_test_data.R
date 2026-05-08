## Purpose: Validate the chapter 5 threshold-correlation workflow on a new dataset.
## Inputs: Excel file with Date, sp500, and bond columns.
## Outputs: Threshold-correlation CSV and PNG in outputs/.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg))
  repo_root <- normalizePath(file.path(dirname(script_path), ".."))
} else {
  repo_root <- getwd()
}

source(file.path(repo_root, "R", "00_setup.R"))

input_path <- "/Users/fernandoperezmillan/Desktop/test_data_corr.xlsx"
output_dir <- file.path(repo_root, "outputs")

# The provided workbook appears to contain price levels, not precomputed log returns.
input_is_log_returns <- FALSE

normalize_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}

pick_column <- function(df_names, candidates, label) {
  normalized <- normalize_name(df_names)
  candidate_idx <- match(normalize_name(candidates), normalized)
  candidate_idx <- candidate_idx[!is.na(candidate_idx)]
  if (length(candidate_idx) == 0) {
    stop(sprintf("Could not find a %s column in %s.", label, input_path))
  }
  df_names[candidate_idx[1]]
}

sheet_name <- excel_sheets(input_path)[1]
raw_df <- read_excel(input_path, sheet = sheet_name)

date_col <- pick_column(names(raw_df), c("Date"), "date")
sp_col <- pick_column(names(raw_df), c("sp500", "sp_500", "sp"), "S&P 500")
bond_col <- pick_column(names(raw_df), c("bond", "bon"), "bond")

df <- data.frame(
  Date = as.Date(raw_df[[date_col]]),
  sp500 = as.numeric(raw_df[[sp_col]]),
  bond = as.numeric(raw_df[[bond_col]])
)

df <- df[complete.cases(df), ]
df <- df[order(df$Date), ]

sp500_xts <- xts(df$sp500, order.by = df$Date)
bond_xts <- xts(df$bond, order.by = df$Date)

if (input_is_log_returns) {
  sp500_returns <- na.omit(sp500_xts)
  bond_returns <- na.omit(bond_xts)
} else {
  sp500_returns <- na.omit(diff(log(sp500_xts)))
  bond_returns <- na.omit(diff(log(bond_xts)))
}

garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)

garch_fit_sp500 <- ugarchfit(spec = garch_spec, data = sp500_returns)
garch_fit_bond <- ugarchfit(spec = garch_spec, data = bond_returns)

z_garch_sp500 <- residuals(garch_fit_sp500, standardize = TRUE)
z_garch_bond <- residuals(garch_fit_bond, standardize = TRUE)

aligned <- na.omit(merge(z_garch_sp500, z_garch_bond, join = "inner"))
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
    NA_real_
  }
})

df_threshold <- data.frame(
  Quantile = quantiles,
  Correlation = as.numeric(threshold_corr)
)

threshold_plot <- ggplot(df_threshold, aes(x = Quantile, y = Correlation)) +
  geom_line(color = "deepskyblue3", linewidth = 1, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Quantile",
    y = "Threshold correlation"
  ) +
  theme_minimal(base_family = "Times") +
  theme_tufte(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    text = element_text(size = 18)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(output_dir, "test_data_threshold_correlation.csv")
png_path <- file.path(output_dir, "test_data_threshold_correlation.png")

write.csv(df_threshold, csv_path, row.names = FALSE)
ggsave(
  filename = png_path,
  plot = threshold_plot,
  width = 8,
  height = 5,
  dpi = 300
)

print(threshold_plot)

cat("Sheet used:", sheet_name, "\n")
cat("Input treated as", if (input_is_log_returns) "log returns" else "price levels", "\n")
cat("Aligned observations:", nrow(aligned), "\n")
cat("Non-missing threshold-correlation points:", sum(!is.na(df_threshold$Correlation)), "\n")
cat("Threshold-correlation CSV:", csv_path, "\n")
cat("Threshold-correlation PNG:", png_path, "\n")
