
# GARCH Analysis of S&P 500 Returns

# Install and load required packages
required_packages <- c(
  "quantmod", 
  "rugarch",   
  "FinTS",     
  "ggplot2",    
  "zoo"         
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(123)


# Data Download & Returns

cat("Downloading S&P 500 data...\n")
getSymbols("^GSPC", from = "2004-01-01", to = "2024-01-01", auto.assign = TRUE)

returns <- dailyReturn(Cl(GSPC), type = "log")[-1]
returns_df <- data.frame(
  Date = index(returns),
  Returns = coredata(returns)
)

cat("Data downloaded successfully. Total observations:", length(returns), "\n\n")


# Model Selection - Grid Search

cat("Starting model selection grid search...\n")

# Initialize results data frame
results <- data.frame(
  ar = integer(),
  ma = integer(),
  p = integer(),
  q = integer(),
  aic = numeric(),
  bic = numeric(),
  model = character(),
  stringsAsFactors = FALSE
)

# Loop through all combinations between 0 and 2 for arma and garch parameters
for (p in 0:2) {   
  for (q in 0:2) {  
    for (ar in 0:2) { 
      for (ma in 0:2) { 
        
        # Skip invalid GARCH(0,0)
        if (p == 0 && q == 0) next
        
        # Specify the ARMA-GARCH model
        spec <- ugarchspec(
          variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
          mean.model = list(armaOrder = c(ar, ma), include.mean = TRUE),
          distribution.model = "norm"
        )
        
        # Attempt to fit the model, handling errors
        fit <- tryCatch(
          ugarchfit(spec = spec, data = returns, solver = "hybrid", 
                    solver.control = list(tol = 1e-3)),
          error = function(e) NULL
        )
        
        # If the fit is successful, extract AIC and BIC
        if (!is.null(fit) && inherits(fit, "uGARCHfit")) {
          ic <- infocriteria(fit)
          aic <- ic[1]  # AIC
          bic <- ic[2]  # BIC
          model_str <- paste0("ARMA(", ar, ",", ma, ")-GARCH(", p, ",", q, ")")
          
          # Add the results to the data frame
          results <- rbind(results, data.frame(
            ar = ar,
            ma = ma,
            p = p,
            q = q,
            aic = aic,
            bic = bic,
            model = model_str,
            stringsAsFactors = FALSE
          ))
          
          cat("Fitted:", model_str, "\n")
        }
      }
    }
  }
}

cat("\nGrid search completed. Total models fitted:", nrow(results), "\n\n")


# Model Ranking

# Sort and select the top 5 models based on AIC
top5_aic <- results[order(results$aic), ][1:5, ]
cat("Top 5 models based on AIC:\n")
print(top5_aic[, c("model", "aic")])
cat("\n")

# Sort and select the top 5 models based on BIC
top5_bic <- results[order(results$bic), ][1:5, ]
cat("Top 5 models based on BIC:\n")
print(top5_bic[, c("model", "bic")])
cat("\n")

# Calculate ranks for AIC and BIC
results$aic_rank <- rank(results$aic, ties.method = "min")
results$bic_rank <- rank(results$bic, ties.method = "min")

# Calculate the sum of ranks
results$sum_rank <- results$aic_rank + results$bic_rank

# Sort by sum_rank and select the top 5 models
top5_joint <- results[order(results$sum_rank), ][1:5, ]

# Display the top 5 models
cat("Top 5 models based on joint AIC and BIC ranks:\n")
print(top5_joint[, c("model", "aic", "bic", "aic_rank", "bic_rank", "sum_rank")])
cat("\n")


# ARCH-LM Tests

cat("Performing ARCH-LM tests...\n\n")

# Define GARCH(1,1) specification
spec11 <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "norm"
)

# Fit GARCH(1,1) model
fit11 <- ugarchfit(spec = spec11, data = returns, solver = "hybrid")
std_res11 <- residuals(fit11, standardize = TRUE)

# Define GARCH(1,2) specification
spec12 <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 2)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "norm"
)

# Fit GARCH(1,2) model
fit12 <- ugarchfit(spec = spec12, data = returns, solver = "hybrid")
std_res12 <- residuals(fit12, standardize = TRUE)

# Set number of lags for the test
lags <- 10

# ARCH-LM test for GARCH(1,1)
arch_test11 <- ArchTest(std_res11^2, lags = lags)

# ARCH-LM test for GARCH(1,2)
arch_test12 <- ArchTest(std_res12^2, lags = lags)

# Display results
cat("ARCH-LM Test Results:\n")
cat(" - GARCH(1,1) p-value:", arch_test11$p.value, "\n")
cat(" - GARCH(1,2) p-value:", arch_test12$p.value, "\n\n")

# Model chosen - ARMA(1,1)-GARCH(1,1) 


# Rolling Window Forecast

cat("Running rolling window forecast...\n")

# Define GARCH(1,1) specification
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1))
)

# Calculate number of forecasts possible
n <- length(returns)
window_size <- 252
n_forecasts <- n - window_size

# Rolling forecast
roll <- ugarchroll(
  spec, 
  data = returns, 
  n.ahead = 1,
  forecast.length = n_forecasts,
  refit.every = 20,
  refit.window = "moving",
  window.size = window_size,
  calculate.VaR = FALSE
)

# Extract forecasts
scond_volatility <- roll@forecast$density[, "Sigma"]^2

# Real values are the observations after the initial window
realized_volatility <- returns[253:length(returns)]^2

# Get dates
dates <- index(returns)[253:length(returns)]

cat("Rolling forecast completed. Forecasts generated:", length(scond_volatility), "\n\n")


# Visualization

cat("Creating visualization...\n")

# Create a data frame with the date, realized volatility, and forecasted volatility
vol_df_combined <- data.frame(
  Date = dates,
  Realized = as.numeric(realized_volatility),
  Forecast = as.numeric(scond_volatility)
)

# Create the plot
p <- ggplot(vol_df_combined, aes(x = Date)) +
  geom_line(aes(y = Realized, color = "Realized Volatility")) +
  geom_line(aes(y = Forecast, color = "Conditional Volatility")) +
  labs(
    title = "GARCH(1,1) Realized vs Conditional Volatility",
    x = "Date",
    y = "Values",
    color = "Legend"
  ) +
  scale_color_manual(
    values = c("Realized Volatility" = "blue",
               "Conditional Volatility" = "red")
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(p)


# Export Data

cat("Exporting data to CSV files...\n")

# Prepare data frames for export
returns_export <- data.frame(
  Date = index(returns), 
  daily.returns = coredata(returns)
)

vol_df <- data.frame(
  Date = dates,
  Conditional_Volatility = as.numeric(scond_volatility)
)

vol_df_real <- data.frame(
  Date = dates,
  Realized_Volatility = as.numeric(realized_volatility)
)

# Export to CSV
write.csv(returns_export, file = "returns_data.csv", row.names = FALSE)
write.csv(vol_df, file = "conditional_volatility.csv", row.names = FALSE)
write.csv(vol_df_real, file = "realized_volatility.csv", row.names = FALSE)

cat("Data exported successfully!\n")
cat(" - returns_data.csv\n")
cat(" - conditional_volatility.csv\n")
cat(" - realized_volatility.csv\n\n")


# Summary Statistics

cat("Summary Statistics:\n")
cat("Returns - Mean:", mean(returns), "SD:", sd(returns), "\n")
cat("Realized Volatility - Mean:", mean(realized_volatility), "\n")
cat("Conditional Volatility - Mean:", mean(scond_volatility), "\n")
