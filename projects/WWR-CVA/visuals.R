library(ggplot2)

# Plot 1 (p1): GBM Simulated Paths 
set.seed(42)
plot_idx <- sample(1:N, 100)
plot_df  <- data.frame(
  time  = rep(time_grid, each = 100),
  price = as.vector(t(S_paths[plot_idx, ])),
  path  = rep(1:100, times = steps + 1)
)

p1 <- ggplot(plot_df, aes(x = time, y = price, group = path)) +
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
  theme_minimal(base_size = 13)

# Plot 2 (p2): Default Time Distribution 
tau_df <- data.frame(tau = tau)

p2 <- ggplot(tau_df, aes(x = tau)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80,
                 fill = "steelblue", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = T, colour = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = T + 1.5, y = 0.025,
           label = paste("T =", T, "yr"), colour = "firebrick", size = 4) +
  labs(
    title    = "Simulated Default Time Distribution",
    subtitle = paste0("Exponential | λ = ", lambda,
                      " | E[τ] = ", round(1/lambda, 0), " years"),
    x        = "Default Time τ (years)",
    y        = "Density"
  ) +
  theme_minimal(base_size = 13)

# Plot 3 (p3): EPE and PFE profiles 
exposure_long_df <- rbind(
  data.frame(time = time_grid, value = EPE, metric = "EPE (Mean Exposure)"),
  data.frame(time = time_grid, value = PFE, metric = "PFE (95th Percentile)")
)

peak_PFE <- max(PFE)

p3 <- ggplot(exposure_long_df, aes(x = time, y = value, colour = metric)) +
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

# Plot 4 (p4): CVA Sensitivity to Correlation 
p4 <- ggplot(results_df, aes(x = rho, y = CVA)) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_point(colour = "steelblue", size = 3) +
  geom_hline(yintercept = CVA_baseline, colour = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = -0.55, y = CVA_baseline * 1.12,
           label = paste("Baseline CVA =", round(CVA_baseline, 4)),
           colour = "firebrick", size = 3.8) +
  annotate("text", x = 0.55, y = 0.15,
           label = "Right-Way Risk", colour = "darkgreen", size = 3.5) +
  annotate("text", x = 0.55, y = 0.72,
           label = "Wrong-Way Risk", colour = "firebrick", size = 3.5) +
  labs(
    title    = "CVA Sensitivity to Wrong-Way Risk",
    subtitle = "Gaussian copula | ρ drives exposure–default dependence",
    x        = "Copula Correlation ρ",
    y        = "CVA"
  ) +
  theme_minimal(base_size = 13)

# Plot 5 (p5): Loss Distribution rho=0 vs rho=0.75 
get_losses <- function(rho) {
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
  losses
}

losses_rho0  <- get_losses(0)
losses_rho75 <- get_losses(0.75)

loss_compare_df <- rbind(
  data.frame(loss = losses_rho0[losses_rho0 > 0],   rho = "ρ = 0.00 (No WWR)"),
  data.frame(loss = losses_rho75[losses_rho75 > 0], rho = "ρ = 0.75 (WWR)")
)

p5 <- ggplot(loss_compare_df, aes(x = loss, fill = rho)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, alpha = 0.6,
                 position = "identity", colour = "white") +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  labs(
    title    = "Loss Distribution: No WWR vs Wrong-Way Risk",
    subtitle = "Defaulted paths only | ρ = 0.75 shifts mass to higher losses",
    x        = "Discounted Loss at Default",
    y        = "Density",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# Plot 6 (p6): Gaussian vs t-copula comparison 
results_copula_df <- readRDS("results_copula_df.rds")

p6 <- ggplot(results_copula_df, aes(x = rho, y = CVA,
                                     colour   = copula,
                                     linetype = copula)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = CVA_baseline, colour = "grey40",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = -0.55, y = CVA_baseline * 1.10,
           label = paste("Baseline CVA =", round(CVA_baseline, 4)),
           colour = "grey40", size = 3.5) +
  scale_colour_manual(values = c(
    "Gaussian"   = "steelblue",
    "t (ν = 4)"  = "firebrick",
    "t (ν = 10)" = "darkorange"
  )) +
  labs(
    title    = "CVA Comparison: Gaussian vs t-Copula",
    subtitle = "t-copula assigns more weight to joint tail events | lower ν = fatter tails",
    x        = "Copula Correlation ρ",
    y        = "CVA",
    colour   = "Copula",
    linetype = "Copula"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# Print all plots 
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

# Console summary 
cat("\nAll plots rendered.\n")
cat("Baseline CVA:                   ", round(CVA_baseline, 4), "\n")
cat("WWR CVA at rho=0.75 (Gaussian): ",
    round(results_df$CVA[results_df$rho == 0.75], 4),
    "—", round(results_df$CVA[results_df$rho == 0.75] / CVA_baseline, 2),
    "x baseline\n")
cat("WWR CVA at rho=0.75 (t nu=4):   ",
    round(results_copula_df$CVA[results_copula_df$rho == 0.75 &
                                  results_copula_df$copula == "t (ν = 4)"], 4),
    "—", round(results_copula_df$CVA[results_copula_df$rho == 0.75 &
                                       results_copula_df$copula == "t (ν = 4)"] / CVA_baseline, 2),
    "x baseline\n")

# ── Save all plots ────────────────────────────────────────────
ggsave("outputs/plots/p1_gbm_paths.png",         p1, width = 8, height = 5, dpi = 150, bg = "white")
ggsave("outputs/plots/p2_default_distribution.png", p2, width = 8, height = 5, dpi = 150, bg = "white")
ggsave("outputs/plots/p3_epe_pfe.png",           p3, width = 8, height = 5, dpi = 150, bg = "white")
ggsave("outputs/plots/p4_cva_sensitivity.png",   p4, width = 8, height = 5, dpi = 150, bg = "white")
ggsave("outputs/plots/p5_loss_distribution.png", p5, width = 8, height = 5, dpi = 150, bg = "white")
ggsave("outputs/plots/p6_copula_comparison.png", p6, width = 8, height = 5, dpi = 150, bg = "white")