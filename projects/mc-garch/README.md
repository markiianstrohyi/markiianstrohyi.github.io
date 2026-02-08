# Monte Carlo Value-at-Risk with GARCH Volatility

## Overview
This project implements a Monte Carlo Value-at-Risk (VaR) model for a multi-asset portfolio, incorporating **time-varying volatility using GARCH(1,1)** with t-student distribution.

---

## Libraries
- **Python**: numpy, pandas, yfinance, arch, matplotlib
- **Notebook**: Jupyter Notebook

---

## What it does
- Constructs 5-asset equaly weighted portfolio
- Builds one day ahead covariance matrix with predicted GARCH volatility 
- Simulates multi-asset portfolio outcomes using Monte Carlo with sampling from t-student distribution
- Estimates next day portfolio VaR for 95% and 99% confidence level
- Performes Kupiec and Christoffersen tests
- Compares with other VaR approaches, namely parametric and EWMA (RiskMetrics)

---

## How to Run
1. Clone the repo:
```bash
git clone https://github.com/markiianstrohyi/markiianstrohyi.github.io.git
cd projects/mc-garch
