if (exists("results_df")) return()
library(ggplot2)
source("exposure.R")
source("default_sim.R")

# Pre-generate base random draws (Gaussian copula) 
set.seed(999)
Z1_base <- rnorm(N)
Z_indep <- rnorm(N)

# Gaussian Copula CVA 
cva_wwr <- function(rho) {

  Z2       <- rho * Z1_base + sqrt(1 - rho^2) * Z_indep
  U_credit <- pnorm(Z2)
  tau_wwr  <- -log(U_credit) / lambda

  losses <- numeric(N)

  for (i in 1:N) {
    if (tau_wwr[i] > T) {
      losses[i] <- 0
    } else {
      tau_i         <- tau_wwr[i]
      S_at_default  <- S0 * exp((mu - 0.5 * sigma^2) * tau_i +
                                   sigma * sqrt(tau_i) * Z1_base[i])
      tau_remaining <- T - tau_i
      exp_i         <- bs_call(S_at_default, K, r, sigma, tau_remaining)
      losses[i]     <- LGD * exp_i * exp(-r * tau_i)
    }
  }

  mean(losses)
}

# t-Copula CVA 
# Extends Gaussian copula by introducing tail dependence
# nu = degrees of freedom: lower = fatter tails = stronger tail dependence
# nu -> infinity converges to Gaussian copula

cva_wwr_t <- function(rho, nu = 4) {

  # Common chi-squared scaling factor - shared across both variables
  # When W is large (tail event), both market and credit are stressed together
  W  <- sqrt(nu / rchisq(N, df = nu))

  # Correlated standard normals
  Z1_t <- rnorm(N)
  Z2_t <- rho * Z1_t + sqrt(1 - rho^2) * rnorm(N)

  # Scale to get correlated t variables
  T1 <- Z1_t * W    # market t-factor
  T2 <- Z2_t * W    # credit t-factor

  # Transform to uniforms using t-CDF
  U_credit <- pt(T2, df = nu)

  # Map to default times
  tau_wwr_t <- -log(U_credit) / lambda

  losses <- numeric(N)

  for (i in 1:N) {
    if (tau_wwr_t[i] > T) {
      losses[i] <- 0
    } else {
      tau_i         <- tau_wwr_t[i]
      S_at_default  <- S0 * exp((mu - 0.5 * sigma^2) * tau_i +
                                   sigma * sqrt(tau_i) * T1[i])
      tau_remaining <- T - tau_i
      exp_i         <- bs_call(S_at_default, K, r, sigma, tau_remaining)
      losses[i]     <- LGD * exp_i * exp(-r * tau_i)
    }
  }

  mean(losses)
}

# Sensitivity: Gaussian copula 
rho_grid   <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 0.90)
CVA_gaussian <- sapply(rho_grid, cva_wwr)

results_df <- data.frame(rho = rho_grid, CVA = CVA_gaussian)

# Sensitivity: t-copula (nu=4 and nu=10) 
set.seed(999)
CVA_t4  <- sapply(rho_grid, function(r) cva_wwr_t(rho = r, nu = 4))

set.seed(999)
CVA_t10 <- sapply(rho_grid, function(r) cva_wwr_t(rho = r, nu = 10))

#  Combined results 
results_copula_df <- rbind(
  data.frame(rho = rho_grid, CVA = CVA_gaussian, copula = "Gaussian"),
  data.frame(rho = rho_grid, CVA = CVA_t4,       copula = "t (ν = 4)"),
  data.frame(rho = rho_grid, CVA = CVA_t10,      copula = "t (ν = 10)")
)

#  Console output 
cat("\n=====================================================\n")
cat("  CVA COMPARISON: GAUSSIAN vs t-COPULA\n")
cat("=====================================================\n")
print(round(
  data.frame(
    rho      = rho_grid,
    Gaussian = CVA_gaussian,
    t_nu4    = CVA_t4,
    t_nu10   = CVA_t10
  ), 4))
cat("=====================================================\n")

# Gaussian sensitivity plot 
ggplot(results_df, aes(x = rho, y = CVA)) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_point(colour = "steelblue", size = 3) +
  geom_hline(yintercept = CVA_baseline, colour = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = -0.55, y = CVA_baseline * 1.08,
           label = paste("Baseline CVA =", round(CVA_baseline, 4)),
           colour = "firebrick", size = 3.5) +
  labs(
    title    = "CVA Sensitivity to Wrong-Way Risk Correlation",
    subtitle = "Gaussian copula | ρ > 0 → Wrong-Way Risk | ρ < 0 → Right-Way Risk",
    x        = "Copula Correlation ρ",
    y        = "CVA"
  ) +
  theme_minimal()

# Save RDS for visuals 
saveRDS(results_df,        "results_df.rds")
saveRDS(results_copula_df, "results_copula_df.rds")