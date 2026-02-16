# Data Dictionary

## `data/raw/final_bond_data.xlsx`

Expected fields used by the code:

- `Dates`: observation date, converted to `Date`.
- `PX_LAST`: bond index level used to compute log returns.

## Downloaded market series

- `^GSPC` from Yahoo Finance.
- Closing price extracted with `Cl()`.
- Daily log returns computed as `diff(log(price))`.
