# Master script: sources setup, helper functions, and chapter scripts in order.
# Run from repository root:
# Rscript scripts/00_run_all.R

source('R/00_setup.R')
source('R/distributions_asyt.R')
source('R/evt_utils.R')
source('R/backtest_utils.R')
source('R/diagnostics_utils.R')

source('scripts/01_chapter1_empirical.R')
source('scripts/02_chapter2_volatility.R')
source('scripts/03_chapter3_shock_modeling.R')
source('scripts/04_chapter4_forecasting_evaluation.R')
source('scripts/05_chapter5_copula_multivariate.R')
