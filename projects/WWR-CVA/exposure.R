if (!exists("exposure")) {
  library(ggplot2)
  source("market_sim.R")

  # Option Parameters 
  K <- 100     # strike price (at-the-money at inception)
  r <- 0.05    # risk-free rate (annual)

  # Black-Scholes Call Price Function 
  bs_call <- function(S, K, r, sigma, tau) {
    # tau = time remaining to maturity
    # returns 0 if expired
    if (tau <= 0) return(rep(0, length(S)))
    
    d1 <- (log(S / K) + (r + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma * sqrt(tau)
    
    price <- S * pnorm(d1) - K * exp(-r * tau) * pnorm(d2)
    return(pmax(price, 0))
  }

  # Compute Exposure at Every Time Step
  # Exposure = MTM value of the call = max(V_t, 0)

  # Pre-allocate exposure matrix: same dims as S_paths
  exposure <- matrix(0, nrow = N, ncol = steps + 1)

  for (j in 1:(steps + 1)) {
    tau_remaining <- T - time_grid[j]
    exposure[, j] <- bs_call(S_paths[, j], K, r, sigma, tau_remaining)
  }


  # Expected Positive Exposure (EPE) 
  # EPE(t) = average exposure across all paths at time t
  # This is what feeds directly into CVA calculation
  EPE <- colMeans(exposure)

  # Potential Future Exposure (PFE) 
  # PFE(t) = 95th percentile exposure at time t
  # Used by credit risk desk for counterparty limit setting
  PFE <- apply(exposure, 2, quantile, probs = 0.95)

  # Summary numbers 
  peak_PFE <- max(PFE)
  cat("Peak PFE (95th percentile):", round(peak_PFE, 4), "\n")
  cat("EPE at t=0:                ", round(EPE[1], 4), "\n")

  # Plot: EPE and PFE profiles together 
  exposure_long_df <- rbind(
    data.frame(time = time_grid, value = EPE, metric = "EPE (Mean Exposure)"),
    data.frame(time = time_grid, value = PFE, metric = "PFE (95th Percentile)")
  )

  ggplot(exposure_long_df, aes(x = time, y = value, colour = metric)) +
    geom_line(linewidth = 1.2) +
    geom_area(data = subset(exposure_long_df, metric == "EPE (Mean Exposure)"),
              aes(y = value), fill = "steelblue", alpha = 0.10) +
    scale_colour_manual(values = c(
      "EPE (Mean Exposure)"   = "steelblue",
      "PFE (95th Percentile)" = "firebrick"
    )) +
    annotate("text", x = 0.6, y = peak_PFE * 1.04,
            label = paste("Peak PFE =", round(peak_PFE, 2)),
            colour = "firebrick", size = 3.5) +
    labs(
      title    = "Exposure Profiles: EPE and PFE",
      subtitle = "European call | ATM | EPE = mean | PFE = 95th percentile",
      x        = "Time (years)",
      y        = "Exposure",
      colour   = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
}
