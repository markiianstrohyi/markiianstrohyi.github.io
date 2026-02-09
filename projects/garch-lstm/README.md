# Forecasting Volatility and Risk with a Dynamically Weighted GARCH-LSTM Ensemble Model

Bachelor thesis project investigating volatility forecasting and risk estimation using S&P 500 daily return data (2004-2024).

## Overview

This repository contains the implementation and analysis of three volatility forecasting approaches:

- **GARCH(1,1)**: Classical econometric model with rolling window estimation
- **LSTM**: Deep learning model with walk-forward validation  
- **GARCH-LSTM Ensemble**: Dynamically weighted hybrid model with inverse-error weighting

## Repository Structure

```
├── data/                    # S&P 500 return data and volatility series
├── GARCH.R                  # GARCH model implementation and rolling forecasts
├── GARCH-LSTM.ipynb         # LSTM tuning, training, ensemble forecasting, VaR/ES estimation
└── Paper.pdf                # Complete thesis document
```

## Methodology

### GARCH Model

- Rolling window approach (252-day window, 20-day refit frequency)
- One-step-ahead conditional volatility forecasts
- Model selection via AIC/BIC and ARCH-LM diagnostic tests

### LSTM Network

- Walk-forward validation (1008-day training with 33% validation split, 504-day test windows)
- Hyperparameter tuning via Random Search
- Sequential input structure for capturing temporal dependencies

### Ensemble Model

- Dynamic weighting based on inverse Mean Absolute Error (MAE)
- Adaptive re-weighting over rolling windows
- Combines GARCH interpretability with LSTM flexibility

### Risk Estimation

- Value-at-Risk (VaR) and Expected Shortfall (ES) at 95% confidence level
- Backtesting via Kupiec POF, Christoffersen Conditional Coverage, and McNeil-Frey tests

## Key Findings

- **LSTM** achieved the lowest point forecast errors (MSE, RMSE, MAE, HMSE) but produced the highest number of VaR violations and failed key backtests

- **Ensemble** model outperformed both individual models in point forecast accuracy metrics but did not fully satisfy backtesting criteria

- **GARCH(1,1)** demonstrated the most reliable tail risk estimation, coming closest to meeting statistical requirements for VaR coverage despite higher point forecast errors

### Conclusion

Improved point forecast accuracy for realized volatility does not necessarily translate into superior risk estimation. The ensemble offers a promising middle ground, though GARCH remains most reliable for regulatory risk applications.

## Data

Daily S&P 500 closing prices from Yahoo Finance (2004-01-01 to 2024-01-01), converted to logarithmic returns. Realized volatility proxied by squared returns.

## Evaluation Metrics

**Forecast Accuracy:**
- Mean Squared Error (MSE)
- Root Mean Squared Error (RMSE)
- Mean Absolute Error (MAE)
- Heteroscedastic MSE (HMSE)

**Risk Validation:**
- Kupiec Proportion of Failures (POF) test
- Christoffersen Conditional Coverage test
- McNeil-Frey test for Expected Shortfall

## Reference

Full methodology, literature review, and empirical results available in [`Paper.pdf`](Paper.pdf).

---
