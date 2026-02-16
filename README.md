# Thesis Repository

# Modeling Financial Returns: Conditional Volatility, Heavy Tails, and Multivariate Dependence

## Abstract

This thesis develops a comprehensive, sequential framework for modeling financial returns, designed as an applied "recipe" that starts with univariate models and culminates in a multivariate, copula-based risk approach. The study begins by documenting stylized facts of returns and then models conditional volatility using GARCH(1,1) and GJR-GARCH(1,1) specifications. To capture heavy tails and asymmetry in standardized shocks, it compares two strategies: a skewed t distribution for the entire density and an extreme value theory (POT/EVT) approach focused on the left tail. The evaluation combines statistical and risk management criteria: out-of-sample forecasts and metrics such as VaR and Expected Shortfall, along with calibration diagnostics and tests based on VaR breaches. Results show that models with a skewed t distribution perform poorly in terms of VaR coverage and breach independence, while EVT-based models improve coverage but still exhibit independence issues, suggesting potential for hybrid specifications that adapt asymmetry and kurtosis more quickly to volatility shocks. Finally, to model dependence across assets, a copula approach (static over standardized returns) is proposed, enabling flexible estimation of marginal models by asset and capturing nonlinear dependence relevant for risk, offering a practical alternative to the linear emphasis of approaches such as DCC.


## Structure

- `scripts/00_run_all.R`: master script that runs the full workflow in order.
- `scripts/01_chapter1_empirical.R` to `scripts/05_chapter5_copula_multivariate.R`: chapter-level analysis scripts.
- `R/00_setup.R`: shared package loading, seed, and shared paths.
- `R/distributions_asyt.R`: asymmetric Student's t helper functions.
- `R/evt_utils.R`: EVT helper functions.
- `R/backtest_utils.R`: backtesting/evaluation helper functions.
- `R/diagnostics_utils.R`: diagnostics helper functions.
- `data/raw/final_bond_data.xlsx`: bond input data.
- `docs/thesis_draft.pdf`: thesis draft.
- `legacy/thesis_original.R`: unchanged monolithic source script.
- `outputs/figures/`: placeholder for saved figures.

## Run

From the repository root:

```bash
Rscript scripts/00_run_all.R
```

## Documentation

- `docs/reproducibility.md`
- `docs/data_dictionary.md`
- `docs/chapter_map.md`
