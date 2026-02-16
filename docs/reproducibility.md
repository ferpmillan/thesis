# Reproducibility Notes

## Environment

- Language: R
- Main execution script: `scripts/00_run_all.R`
- Seed: `191023` (set in `R/00_setup.R` ans in various parts of the code)

## Required packages

`readxl`, `quantmod`, `ggplot2`, `ggthemes`, `forecast`, `MASS`, `rugarch`, `rmgarch`, `reshape2`, `fGarch`, `evir`, `xts`, `tseries`, `gridExtra`, `zoo`, `copula`, `matrixcalc`, `pracma`, `MCMCpack`, `viridis`

## Data sources

- Local bond data: `data/raw/final_bond_data.xlsx`
- Market data downloaded via Yahoo Finance through `quantmod::getSymbols`

## Execution order

1. `R/00_setup.R`
2. Function modules under `R/`
3. Chapter scripts `scripts/01_...` to `scripts/05_...`

