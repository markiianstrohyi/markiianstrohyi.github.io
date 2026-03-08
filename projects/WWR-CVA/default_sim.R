if (!exists("tau")) {
  library(ggplot2)
  source("market_sim.R")  

  # Credit Parameters
  lambda <- 0.05   # hazard rate (5% annual default intensity)
                  # implies ~5% probability of default in year 1
  R_rec  <- 0.40   # recovery rate (40% as standard market convention)
                  # LGD = 1 - R_rec = 0.60

  # Simulate Default Times 
  # From hazard rate model: tau = -ln(U) / lambda
  # U ~ Uniform(0,1)

  set.seed(123)
  U      <- runif(N)
  tau    <- -log(U) / lambda    # vector of N default times (in years)


  # Plot: Default Time Distribution 
  tau_df <- data.frame(tau = tau)

  ggplot(tau_df, aes(x = tau)) +
    geom_histogram(aes(y = after_stat(density)),
                  bins = 80,
                  fill = "steelblue", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = T, colour = "firebrick",
              linetype = "dashed", linewidth = 1) +
    annotate("text", x = T + 0.5, y = 0.03,
            label = paste("T =", T, "yr"), colour = "firebrick", size = 4) +
    labs(
      title    = "Simulated Default Time Distribution",
      subtitle = paste0("Exponential | λ = ", lambda,
                        " | E[τ] = ", round(1/lambda, 0), " years"),
      x        = "Default Time τ (years)",
      y        = "Density"
    ) +
    theme_minimal()

}
