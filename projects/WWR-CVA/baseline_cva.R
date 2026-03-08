if (!exists("CVA_baseline")) {
  library(ggplot2)
  source("exposure.R")    
  source("default_sim.R") 

  # CVA Calculation 
  # CVA = (1 - R) * (1/N) * sum over all paths of:
  #       exposure at default time * discount factor

  # For each path i:
  #   - if tau[i] > T : counterparty survives, no loss
  #   - if tau[i] <= T: look up exposure at the time step closest to tau[i]


  LGD <- 1 - R_rec   

  # Map each default time to the nearest column index in exposure matrix
  # time_grid goes from 0 to T in (steps+1) points
  get_exposure_at_default <- function(tau_i) {
    if (tau_i > T) return(0)   # survived past horizon, no CVA contribution
    
    # Find closest time step index
    idx <- which.min(abs(time_grid - tau_i))
    return(exposure[, idx])    # returns exposure for ALL paths at that time
  }

  # Vectorised approach: for each path, get its own exposure at its own tau
  losses <- numeric(N)

  for (i in 1:N) {
    if (tau[i] > T) {
      losses[i] <- 0
    } else {
      idx        <- which.min(abs(time_grid - tau[i]))
      losses[i]  <- LGD * exposure[i, idx] * exp(-r * tau[i])
    }
  }

  CVA_baseline <- mean(losses)

  # ── Results ───────────────────────────────────────────────────
  cat("==========================================\n")
  cat("  BASELINE CVA (No Wrong-Way Risk)\n")
  cat("==========================================\n")
  cat("CVA:                ", round(CVA_baseline, 4), "\n")
  cat("LGD:                ", LGD, "\n")
  cat("Lambda:             ", lambda, "\n")
  cat("Recovery Rate:      ", R_rec, "\n")
  cat("Paths defaulting:   ", sum(tau <= T), "/", N, "\n")
  cat("==========================================\n")

  # ── Plot: Loss Distribution (defaulted paths only) ────────────
  defaulted_losses <- losses[losses > 0]

  loss_df <- data.frame(loss = defaulted_losses)

  ggplot(loss_df, aes(x = loss)) +
    geom_histogram(aes(y = after_stat(density)),
                  bins  = 60,
                  fill  = "steelblue",
                  colour = "white",
                  alpha = 0.8) +
    geom_vline(xintercept = CVA_baseline, colour = "firebrick",
              linetype = "dashed", linewidth = 1) +
    annotate("text",
            x     = CVA_baseline * 1.15,
            y     = Inf,
            vjust = 1.5,
            label = paste("CVA =", round(CVA_baseline, 4)),
            colour = "firebrick",
            size  = 4) +
    labs(
      title    = "Distribution of Discounted Losses at Default",
      subtitle = "Baseline model — independence assumption | Defaulted paths only",
      x        = "Discounted Loss",
      y        = "Density"
    ) +
    theme_minimal()

}
