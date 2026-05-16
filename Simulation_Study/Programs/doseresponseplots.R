# Load necessary library
library(ggplot2)

# 1. Define the continuous dose range from 0 to 10
d <- seq(0, 8, length.out = 500)

# 2. Calculate the true mean response (mu) for each scenario
mu_linear    <- -4.59 + 7.33 * d
mu_quadratic <- 0 + 28 * d - 4 * d^2
mu_emax      <- -10 + (100 * d) / (0.5 + d)

# 3. Combine into a single data frame for ggplot
plot_data <- data.frame(
  Dose = rep(d, 3),
  Response = c(mu_linear, mu_quadratic, mu_emax),
  Scenario = factor(rep(c("Linear", "Quadratic", "Emax"), each = length(d)),
                    levels = c("Linear", "Quadratic", "Emax")) # Sets the order of the panels
)

# 4. Generate the 3-panel plot
p_curves <- ggplot(plot_data, aes(x = Dose, y = Response)) +
  geom_line(color = "#2B435C", linewidth = 1.2) +
  facet_wrap(~ Scenario, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    title = "True Dose-Response Data-Generating Scenarios",
    x = "Dose (x1000 IU)",
    y = expression(paste("True Mean Response (", mu, ")"))
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# View the plot
print(p_curves)

# Save the plot (optional)
# ggsave("True_Dose_Response_Curves.png", p_curves, width = 10, height = 4, dpi = 300)