if (!exists("S_paths")) {
  library(ggplot2)

  #Parameters 
  S0    <- 100      # initial asset price
  mu    <- 0.05     # drift (annual)
  sigma <- 0.20     # volatility (annual)
  T     <- 1        # time horizon in years
  N     <- 100000    # number of simulation paths
  steps <- 252      # daily steps (trading days in a year)

  dt    <- T / steps

  set.seed(42)

  # Matrix of standard normal shocks: N paths x steps
  Z <- matrix(rnorm(N * steps), nrow = N, ncol = steps)

  # Log returns at each step
  log_returns <- (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z

  # Cumulative sum across columns to build log price path
  log_paths <- t(apply(log_returns, 1, cumsum))

  # Convert to price paths and prepend S0 as column 1
  S_paths <- S0 * exp(cbind(0, log_paths))

  # Plots
  time_grid <- seq(0, T, length.out = steps + 1)

  # Sample 100 paths for plotting 
  plot_idx   <- sample(1:N, 100)
  plot_df    <- data.frame(
    time  = rep(time_grid, each = 100),
    price = as.vector(t(S_paths[plot_idx, ])),
    path  = rep(1:100, times = steps + 1)
  )

  ggplot(plot_df, aes(x = time, y = price, group = path)) +
    geom_line(alpha = 0.15, colour = "steelblue", linewidth = 0.3) +
    stat_summary(aes(group = 1),
                fun = mean, geom = "line",
                colour = "firebrick", linewidth = 1.2) +
    labs(
      title    = "GBM Simulated Asset Paths",
      subtitle = paste0("N = ", N, " paths | σ = ", sigma, " | μ = ", mu),
      x        = "Time (years)",
      y        = "Asset Price"
    ) +
    theme_minimal()

}
