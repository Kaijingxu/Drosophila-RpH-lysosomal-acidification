library(readxl)
library(dplyr)
library(ggplot2)
library(emmeans)
library(sandwich)
library(ggpubr)

# Import data
eye_baf_raw <- readxl::read_excel("data/GMR_RpH/BafA1.xlsx")

# Settings
scale_factor <- 100
eps <- 1e-12

condition_levels <- c("DMSO",
                      "1 uM BafA1",
                      "2.5 uM BafA1",
                      "5 uM BafA1")

condition_colours <- c("DMSO"       = "black",
                       "1 uM BafA1" = "pink",
                       "2.5 uM BafA1" = "pink2",
                       "5 uM BafA1" = "red")

fmt_p <- function(p) {vapply(p,
                             FUN.VALUE = character(1),
                             FUN = function(x) {if (is.na(x)) {return(NA_character_)}
                               if (x < 1e-4) {
                                 paste0("p = ", format(signif(x, 3), scientific = TRUE))} 
                               else if (x < 0.01) {
                                 paste0("p = ", formatC(x, format = "f", digits = 4))} 
                               else if (x < 0.1) {
                                 paste0("p = ", formatC(x, format = "f", digits = 3))} 
                               else {paste0("p = ", formatC(x, format = "f", digits = 2))}})}

# Prepare analysis dataset
eye_baf <- eye_baf_raw %>%
  mutate(Condition = factor(Condition, 
                            levels = condition_levels)) %>%
  filter(puncta_area > 0,
         IntDen_corr_mC > 0,
         IntDen_corr_pH > 0) %>%
  mutate(RpH_ratio = IntDen_corr_pH / IntDen_corr_mC,
         ratio_plot = scale_factor * RpH_ratio,
         log_RpH_ratio = log(RpH_ratio + eps)) %>%
  filter(is.finite(RpH_ratio),
         is.finite(log_RpH_ratio))

# Model and robust inference
fit_eye_baf <- lm(log_RpH_ratio ~ Condition, data = eye_baf)
vcov_eye_baf <- vcovHC(fit_eye_baf,
                       type = "HC3")

emm_eye_baf <- emmeans(fit_eye_baf,
                       ~ Condition,
                       vcov. = vcov_eye_baf)

baf_contrasts <- contrast(emm_eye_baf,
                          method = "trt.vs.ctrl",
                          ref = "DMSO",
                          adjust = "holm") %>%
  summary() %>%
  as.data.frame()

# Plot annotations
observed_max <- max(eye_baf$ratio_plot, na.rm = TRUE)

bracket_positions <- observed_max + observed_max * c(0.18, 0.34, 0.50)
y_limit <- max(bracket_positions) * 1.06

stats_eye_baf <- baf_contrasts %>%
  mutate(group1 = "DMSO",
         group2 = sub(" - DMSO$", "", contrast),
         label = fmt_p(p.value),
         y.position = bracket_positions)

sample_sizes <- eye_baf %>%
  count(Condition, name = "n") %>%
  mutate(y = 0.08 * y_limit,
         label = paste0("n = ", n))

# Panel B
set.seed(1)

figure_4b <- ggplot(eye_baf,
                    aes(x = Condition,
                        y = ratio_plot,
                        fill = Condition)) +
  geom_boxplot(width = 0.48,
               outlier.shape = NA,
               colour = "black",
               coef = Inf) +
  geom_point(position = position_jitter(width = 0.12, height = 0),
             shape = 21,
             size = 3,
             fill = "black",
             colour = "black",
             stroke = 0) +
  geom_text(data = sample_sizes,
            aes(x = Condition,
                y = y,
                label = label),
            inherit.aes = FALSE,
            size = 4) +
  stat_pvalue_manual(stats_eye_baf,
                     label = "label",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  scale_fill_manual(values = condition_colours) +
  scale_x_discrete(labels = c("DMSO",
                              "1 µM\nBafA1",
                              "2.5 µM\nBafA1",
                              "5 µM\nBafA1")) +
  scale_y_continuous(limits = c(0, y_limit),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Condition",
       y = expression("RpH ratio (" * CTF[pHluorin] / CTF[mCherry] * ") × 100"),
       title = "B: Validation with BafA1 titration") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_4b

# Supplementary Figure S9A: Background variability and mask coverage: BafA1 titration
cv_threshold       <- 0.30
coverage_threshold <- 25
eps                 <- 1e-12

# Prepare QC dataset
eye_baf_qc <- eye_baf %>%
  mutate(CV_mCherry = if_else(is.finite(back_mean_mC) & back_mean_mC > eps,
                              SD_mC / back_mean_mC,
                              NA_real_),
         CV_pHluorin = if_else(is.finite(back_mean_pH) & back_mean_pH > eps,
                               SD_pH / back_mean_pH,
                               NA_real_),
         mask_coverage = if_else(is.finite(whole_eye_area) & whole_eye_area > 0,
                                 100 * puncta_area / whole_eye_area,
                                 NA_real_)) %>%
  filter(is.finite(CV_mCherry),
         is.finite(CV_pHluorin),
         is.finite(mask_coverage))

# QC summaries
cv_pass_rate <- mean(eye_baf_qc$CV_mCherry <= cv_threshold & eye_baf_qc$CV_pHluorin <= cv_threshold)

coverage_pass_rate <- mean(eye_baf_qc$mask_coverage >= coverage_threshold)

coverage_median <- median(eye_baf_qc$mask_coverage,
                          na.rm = TRUE)

coverage_iqr <- quantile(eye_baf_qc$mask_coverage,
                         probs = c(0.25, 0.75),
                         na.rm = TRUE,
                         names = FALSE)

qc_subtitle <- sprintf(paste0("Background CV thresholds: mCherry ≤ %.2f and pHluorin ≤ %.2f; %d%% passed\n",
                              "Mask coverage: median %.2f%% (IQR %.2f–%.2f%%); %d%% met the ≥%d%% threshold"),
                       cv_threshold,
                       cv_threshold,
                       round(100 * cv_pass_rate),
                       coverage_median,
                       coverage_iqr[1],
                       coverage_iqr[2],
                       round(100 * coverage_pass_rate),
                       coverage_threshold)

# Axis limits
x_upper <- max(eye_baf_qc$CV_mCherry,
               cv_threshold,
               na.rm = TRUE) * 1.15

y_upper <- max(eye_baf_qc$CV_pHluorin,
               cv_threshold,
               na.rm = TRUE) * 1.15

# Supplementary Figure S9A plot
figure_s9a <- ggplot(eye_baf_qc,
                     aes(x = CV_mCherry,
                         y = CV_pHluorin,
                         colour = Condition)) +
  annotate("rect",
           xmin = 0,
           xmax = cv_threshold,
           ymin = 0,
           ymax = cv_threshold,
           alpha = 0.08) +
  geom_vline(xintercept = cv_threshold,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_hline(yintercept = cv_threshold,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_point(size = 3,
             alpha = 1) +
  annotate("text",
           x = 0.52 * cv_threshold,
           y = 1.06 * cv_threshold,
           label = "Reference thresholds",
           hjust = 0,
           size = 4.5) +
  scale_colour_manual(values = condition_colours,
                      breaks = condition_levels,
                      name = "Condition") +
  scale_x_continuous(limits = c(0, x_upper),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(limits = c(0, y_upper),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(title = "SA: Background CV and mask coverage (BafA1 titration)",
       subtitle = qc_subtitle,
       x = "Background CV, mCherry",
       y = "Background CV, pHluorin") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_s9a