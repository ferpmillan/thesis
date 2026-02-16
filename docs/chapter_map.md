# Chapter to Script Map

## Original source

- `legacy/thesis_original.R` (copy of the monolithic script)

## Modular mapping

- Chapter 1 -> `scripts/01_chapter1_empirical.R`
- Chapter 2 -> `scripts/02_chapter2_volatility.R`
- Chapter 3 -> `scripts/03_chapter3_shock_modeling.R`
- Chapter 4 -> `scripts/04_chapter4_forecasting_evaluation.R`
- Chapter 5 -> `scripts/05_chapter5_copula_multivariate.R`

## Shared function modules

- Asymmetric t distribution helpers: `R/distributions_asyt.R`
- EVT helpers: `R/evt_utils.R`
- Backtesting helpers: `R/backtest_utils.R`
- Diagnostics helpers: `R/diagnostics_utils.R`
