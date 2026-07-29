library(dplyr)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(purrr)
library(readxl)
library(sandwich)
library(tidyr)

# Import and prepare vehicle-control replicate data
df_vehicle_raw <- read_excel("data/cg_RpH/R1-R3.xlsx")

scale_factor <- 1e4
df_vehicle <- df_vehicle_raw %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("+/+", "s64w/+", "s64w/s64w")),
         Replicate = factor(Replicate,
                            levels = c("1", "2", "3"),
                            labels = c("R1", "R2", "R3")),
         CTF_mCherry = IntDen_mC - Area * mean_back_mC,
         CTF_pHluorin = IntDen_pH - Area * mean_back_pH) %>%
  filter(Area > 0,
         CTF_mCherry > 0,
         CTF_pHluorin > 0) %>%
  mutate(RpH_ratio = CTF_pHluorin / CTF_mCherry,
         RpH_ratio_plot = scale_factor * RpH_ratio,
         log_RpH_ratio = log(RpH_ratio))

# Statistical analysis within each replicate
vehicle_contrasts <- df_vehicle %>%
  split(.$Replicate) %>%
  imap_dfr(function(data_rep, replicate_name) {
    fit <- lm(log_RpH_ratio ~ Genotype,
              data = data_rep)
  vcov_hc3 <- sandwich::vcovHC(fit,
                               type = "HC3")
  emmeans(fit,
          ~ Genotype,
          vcov. = vcov_hc3) %>%
    contrast(method = "pairwise",
             adjust = "holm") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(Replicate = replicate_name)}) %>%
  mutate(Replicate = factor(Replicate,
                            levels = c("R1", "R2", "R3")))

format_p <- function(p) {
  case_when(is.na(p) ~ NA_character_,
            p < 1e-4 ~ paste0("p = ",
                              format(signif(p, 3), scientific = TRUE)),
            p < 0.01 ~ paste0("p = ",
                              formatC(p, format = "f", digits = 4)),
            p < 0.10 ~ paste0("p = ",
                              formatC(p, format = "f", digits = 3)),
            TRUE ~ paste0("p = ",
                          formatC(p, format = "f", digits = 2)))}

# Prepare statistical annotations
genotype_levels <- levels(df_vehicle$Genotype)

x_centres <- setNames(seq_along(genotype_levels),
                      genotype_levels)

replicate_ymax <- df_vehicle %>%
  group_by(Replicate) %>%
  summarise(ymax = max(RpH_ratio_plot, na.rm = TRUE),
            .groups = "drop")

stats_vehicle <- vehicle_contrasts %>%
  separate(contrast,
           into = c("group1", "group2"),
           sep = " - ",
           remove = FALSE) %>%
  mutate(group1 = trimws(gsub("[`()]", "", group1)),
         group2 = trimws(gsub("[`()]", "", group2))) %>%
  group_by(Replicate) %>%
  mutate(comparison_order = row_number()) %>%
  ungroup() %>%
  left_join(replicate_ymax,
            by = "Replicate") %>%
  transmute(Replicate,
            group1,
            group2,
            xmin = as.numeric(x_centres[group1]),
            xmax = as.numeric(x_centres[group2]),
            y.position = ymax * (1.12 + 0.12 * (comparison_order - 1)),
            label = format_p(p.value))

# Samples size
sample_sizes_vehicle <- df_vehicle %>%
  count(Replicate,
        Genotype,
        name = "n") %>%
  left_join(replicate_ymax,
            by = "Replicate") %>%
  mutate(x = as.numeric(x_centres[as.character(Genotype)]),
         y = 0.05 * ymax,
         label = paste0("n = ", n))

# Generate vehicle-control replicate Figure 3C
genotype_colours <- c("+/+" = "blue",
                      "s64w/+" = "green",
                      "s64w/s64w" = "red")

plot_ymax <- max(stats_vehicle$y.position,
                 df_vehicle$RpH_ratio_plot,
                 na.rm = TRUE)

set.seed(1)

figure_3c <- ggplot(df_vehicle,
                    aes(x = Genotype,
                        y = RpH_ratio_plot,
                        fill = Genotype)) +
  facet_wrap(~ Replicate,
             nrow = 1) +
  geom_boxplot(width = 0.45,
               outlier.shape = NA,
               colour = "black") +
  geom_jitter(width = 0.15,
              height = 0,
              size = 2,
              alpha = 0.70) +
  geom_text(data = sample_sizes_vehicle,
            aes(x = x,
                y = y,
                label = label),
            inherit.aes = FALSE,
            size = 4) +
  stat_pvalue_manual(stats_vehicle,
                     label = "label",
                     xmin = "xmin",
                     xmax = "xmax",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  scale_fill_manual(values = genotype_colours,
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(0, 1.01 * plot_ymax),
                  clip = "off") +
  labs(x = "Genotype",
       y = expression("RpH ratio (" *  CTF[pHluorin] / CTF[mCherry] * ") " %*%  10^4),
       title = "C: Reproducibility across vehicle-control replicates") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        strip.background = element_rect(fill = "grey85",
                                        colour = "black"),
        strip.text = element_text(face = "plain"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_3c

# Supplementary Figure S8C: cg-RpH fat-body imaging quality control
cbr_threshold <- 5
snr_threshold <- 5

df_qc_vehicle <- df_vehicle %>%
  mutate(net_mCherry = mean_mC - mean_back_mC,
         CBR_mCherry = if_else(mean_back_mC > 0,
                               mean_mC / mean_back_mC,
                               NA_real_),
         SNR_mCherry = if_else(SD_back_mC > 0,
                               net_mCherry / SD_back_mC,
                               NA_real_),
         delta_background = if_else(net_mCherry > 0,
                                    100 * abs(mean_back_g5_mC - mean_back_mC) / net_mCherry,
                                    NA_real_),
         qc_pass = (CBR_mCherry >= cbr_threshold & SNR_mCherry >= snr_threshold)) %>%
  filter(is.finite(CBR_mCherry),
         is.finite(SNR_mCherry),
         is.finite(delta_background))

qc_summary_vehicle <- df_qc_vehicle %>%
  group_by(Replicate) %>%
  summarise(n = n(),
            pass_rate = mean(qc_pass),
            delta_background_median = median(delta_background),
            delta_background_q1 = quantile(delta_background, 0.25),
            delta_background_q3 = quantile(delta_background, 0.75),
            .groups = "drop") %>%
  mutate(qc_label = sprintf("%d%% passed; \u0394B %.2f%% (%.2f\u2013%.2f%%)",
                            round(100 * pass_rate),
                            delta_background_median,
                            delta_background_q1,
                            delta_background_q3))

replicate_labels <- setNames(paste0(qc_summary_vehicle$Replicate,
                                    "\n",
                                    qc_summary_vehicle$qc_label),
                             qc_summary_vehicle$Replicate)

x_range_vehicle <- range(df_qc_vehicle$CBR_mCherry,
                         na.rm = TRUE)

y_range_vehicle <- range(df_qc_vehicle$SNR_mCherry,
                         na.rm = TRUE)

x_limits_vehicle <- c(min(floor(x_range_vehicle[1]) - 1, cbr_threshold - 1),
                      ceiling(x_range_vehicle[2]) + 1)

y_limits_vehicle <- c(max(0,
                          min(floor(y_range_vehicle[1]) - 1, snr_threshold - 2)),
                      ceiling(y_range_vehicle[2]) + 1)

# Supplementary Figure S8C plot
figure_s7c <- ggplot(df_qc_vehicle,
                     aes(x = CBR_mCherry,
                         y = SNR_mCherry,
                         colour = Genotype)) +
  annotate("rect",
           xmin = cbr_threshold,
           xmax = Inf,
           ymin = snr_threshold,
           ymax = Inf,
           fill = "grey50",
           alpha = 0.08) +
  geom_vline(xintercept = cbr_threshold,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_hline(yintercept = snr_threshold,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_point(size = 2.6,
             alpha = 0.85) +
  annotate("text",
           x = cbr_threshold + 0.4,
           y = snr_threshold - 0.8,
           label = "Reference thresholds",
           hjust = 0,
           size = 4) +
  facet_wrap(~ Replicate,
             nrow = 1,
             labeller = as_labeller(replicate_labels)) +
  scale_colour_manual(values = genotype_colours,
                      breaks = names(genotype_colours),
                      name = "Genotype") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = x_limits_vehicle,
                  ylim = y_limits_vehicle,
                  clip = "off") +
  labs(x = "Contrast-to-background ratio, mCherry",
       y = "Signal-to-noise ratio, mCherry",
       title = "SC: Signal and background adequacy (vehicle-control reproducibility)",
       subtitle = "CBR \u2265 5 and SNR \u2265 5") +
  guides(colour = guide_legend(order = 1)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_rect(fill = "grey85",
                                        colour = "black"),
        strip.text = element_text(face = "plain",
                                  size = 9),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_s7c

# Supplementary Figure S8D: Sensitivity checks for vehicle-control replicates
background_threshold <- 0.5
mask_threshold <- 7

df_sensitivity_vehicle <- df_vehicle %>%
  mutate(CTF_mCherry_bg5 = IntDen_mC - Area * mean_back_g5_mC,
         CTF_pHluorin_bg5 = IntDen_pH - Area * mean_back_g5_pH,   
         CTF_mCherry_mask1 = IntDen_mC_tight - Area_tight * mean_back_mC,
         CTF_pHluorin_mask1 = IntDen_pH_tight - Area_tight * mean_back_pH) %>%
  mutate(RpH_ratio_bg5 = if_else(CTF_mCherry_bg5 > 0 & CTF_pHluorin_bg5 > 0,
                                 CTF_pHluorin_bg5 / CTF_mCherry_bg5,
                                 NA_real_),
         RpH_ratio_mask1 = if_else(CTF_mCherry_mask1 > 0 & CTF_pHluorin_mask1 > 0,
                                   CTF_pHluorin_mask1 / CTF_mCherry_mask1,
                                   NA_real_)) %>%
  mutate(background_change = 100 * abs(RpH_ratio_bg5 - RpH_ratio) / RpH_ratio,
         mask_change =  100 * abs(RpH_ratio_mask1 - RpH_ratio) / RpH_ratio,
         sensitivity_pass = background_change <= background_threshold & mask_change <= mask_threshold) %>%
  filter(RpH_ratio > 0,
         is.finite(background_change),
         is.finite(mask_change))

sensitivity_summary_vehicle <- df_sensitivity_vehicle %>%
  group_by(Replicate) %>%
  summarise(n = n(),
            pass_rate = mean(sensitivity_pass),
            background_median = median(background_change),
            background_q1 = quantile(background_change, 0.25),
            background_q3 = quantile(background_change, 0.75),
            mask_median = median(mask_change),
            mask_q1 = quantile(mask_change, 0.25),
            mask_q3 = quantile(mask_change, 0.75),
            .groups = "drop") %>%
  mutate(sensitivity_label = sprintf("%d%% passed",
                                     round(100 * pass_rate)))

replicate_labels_sensitivity <- setNames(paste0(sensitivity_summary_vehicle$Replicate,
                                                "\n",
                                                sensitivity_summary_vehicle$sensitivity_label),
                                         sensitivity_summary_vehicle$Replicate)

x_max_sensitivity_vehicle <- 1.10 * max(mask_threshold,
                                        df_sensitivity_vehicle$mask_change,
                                        na.rm = TRUE)
y_max_sensitivity_vehicle <- 1.10 * max(background_threshold,
                                        df_sensitivity_vehicle$background_change,
                                        na.rm = TRUE)

# Supplementary Figure S8D plot
figure_s7d <- ggplot(df_sensitivity_vehicle,
                     aes(x = mask_change,
                         y = background_change,
                         colour = Genotype)) +
  annotate("rect",
           xmin = 0,
           xmax = mask_threshold,
           ymin = 0,
           ymax = background_threshold,
           fill = "grey50",
           alpha = 0.08) +
  geom_vline(xintercept = mask_threshold,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_hline(yintercept = background_threshold,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_point(size = 2.6,
             alpha = 0.90) +
  annotate("text",
           x = 0.55 * mask_threshold,
           y = 1.04 * background_threshold,
           label = "Reference thresholds",
           hjust = 0.5,
           vjust = 0,
           size = 4) +
  facet_wrap(~ Replicate,
             nrow = 1,
             labeller = as_labeller(replicate_labels_sensitivity)) +
  scale_colour_manual(values = genotype_colours,
                      breaks = names(genotype_colours),
                      name = "Genotype") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +  
  coord_cartesian(xlim = c(0, x_max_sensitivity_vehicle),
                  ylim = c(0, y_max_sensitivity_vehicle),
                  clip = "off") +
  labs(x = "Mask \u22121 px, absolute % change",
       y = "Background +5, absolute % change",
       title = "SD: Sensitivity checks (vehicle-control reproducibility)",
       subtitle = paste0("Background +5 \u2264 ",
                         background_threshold,
                         "% and Mask \u22121 px \u2264 ",
                         mask_threshold,
                         "%")) +
  guides(colour = guide_legend(order = 1)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.background = element_rect(fill = "grey85",
                                        colour = "black"),
        strip.text = element_text(face = "plain",
                                  size = 9),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_s7d