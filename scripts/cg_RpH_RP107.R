library(dplyr)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(readxl)
library(sandwich)
library(tidyr)

# Import and prepare RP107 data
df_rp107_raw <- read_excel("data/cg_RpH/RP107.xlsx")

scale_factor <- 1e4
dodge_width <- 0.60

df_rp107 <- df_rp107_raw %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("+/+", "s64w/+", "s64w/s64w")),
         Condition = factor(Condition,
                            levels = c("DMSO", "RP107")),
         CTF_mCherry = IntDen_mC - Area * mean_back_mC,
         CTF_pHluorin = IntDen_pH - Area * mean_back_pH) %>%
  filter(Area > 0,
         CTF_mCherry > 0,
         CTF_pHluorin > 0) %>%
  mutate(RpH_ratio = CTF_pHluorin / CTF_mCherry,
         RpH_ratio_plot = scale_factor * RpH_ratio,
         log_RpH_ratio = log(RpH_ratio))

# Statistical analysis
fit_rp107 <- lm(log_RpH_ratio ~ Genotype * Condition,
                data = df_rp107)

vcov_rp107_hc3 <- sandwich::vcovHC(fit_rp107,
                                   type = "HC3")

# RP107 versus DMSO within each genotype
emm_condition_rp107 <- emmeans(fit_rp107,
                               ~ Condition | Genotype,
                               vcov. = vcov_rp107_hc3)

condition_contrasts_rp107 <- contrast(emm_condition_rp107,
                                      method = "revpairwise",
                                      adjust = "holm") %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# Genotype comparisons under DMSO
emm_genotype_rp107 <- emmeans(fit_rp107,
                              ~ Genotype | Condition,
                              vcov. = vcov_rp107_hc3)

dmso_genotype_contrasts_rp107 <- contrast(emm_genotype_rp107,
                                          method = "pairwise",
                                          adjust = "holm") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  filter(Condition == "DMSO")

format_p <- function(p) {case_when(is.na(p) ~ NA_character_,
                                   p < 1e-4 ~ paste0("p = ",
                                                     format(signif(p, 3), scientific = TRUE)),
                                   p < 0.01 ~ paste0("p = ",
                                                     formatC(p, format = "f", digits = 4)),
                                   p < 0.10 ~ paste0("p = ",
                                                     formatC(p, format = "f", digits = 3)),
                                   TRUE ~ paste0("p = ",
                                                 formatC(p, format = "f", digits = 2)))}

# Prepare statistical annotations
genotype_levels <- levels(df_rp107$Genotype)

x_centres <- setNames(seq_along(genotype_levels),
                      genotype_levels)

x_offsets <- c(DMSO = -dodge_width / 4,
               RP107 = dodge_width / 4)

y_max_rp107 <- max(df_rp107$RpH_ratio_plot,
                   na.rm = TRUE)

genotype_ymax_rp107 <- df_rp107 %>%
  group_by(Genotype) %>%
  summarise(ymax = max(RpH_ratio_plot, na.rm = TRUE),
            .groups = "drop")

stats_condition_rp107 <- condition_contrasts_rp107 %>%
  transmute(Genotype = as.character(Genotype),
            group1 = "DMSO",
            group2 = "RP107",
            label = format_p(p.value),
            xmin = as.numeric(x_centres[Genotype] + x_offsets["DMSO"]),
            xmax = as.numeric(x_centres[Genotype] + x_offsets["RP107"])) %>%
  left_join(genotype_ymax_rp107,
            by = "Genotype") %>%
  mutate(y.position = 1.12 * ymax)

dmso_pairs_rp107 <- tibble::tibble(group1 = c("+/+", "+/+", "s64w/+"),
                                   group2 = c("s64w/+", "s64w/s64w", "s64w/s64w"),
                                   comparison_order = 1:3)

dmso_contrast_pairs_rp107 <- dmso_genotype_contrasts_rp107 %>%
  tidyr::separate(contrast,
                  into = c("group1", "group2"),
                  sep = " - ",
                  remove = FALSE) %>%
  mutate(group1 = trimws(gsub("[`()]", "", group1)),
         group2 = trimws(gsub("[`()]", "", group2)))

stats_dmso_genotype_rp107 <- dmso_pairs_rp107 %>%
  left_join(dmso_contrast_pairs_rp107 %>%
              select(group1, group2, p.value),
            by = c("group1", "group2")) %>%
  mutate(label = format_p(p.value),
         xmin = as.numeric(x_centres[group1] + x_offsets["DMSO"]),
         xmax = as.numeric(x_centres[group2] + x_offsets["DMSO"]),
         y.position = y_max_rp107 * c(1.12, 1.20, 1.28))

# samples size
sample_sizes_rp107 <- df_rp107 %>%
  count(Genotype,
        Condition,
        name = "n") %>%
  mutate(x = as.numeric(x_centres[as.character(Genotype)]) +
           x_offsets[as.character(Condition)],
         y = 0.08 * y_max_rp107,
         label = paste0("n = ", n))

# Figure 3D: RP107 application
genotype_colours <- c("+/+" = "blue",
                      "s64w/+" = "green",
                      "s64w/s64w" = "red")

condition_shapes_rp107 <- c(DMSO = 21,
                            RP107 = 24)

plot_ymax_rp107 <- max(df_rp107$RpH_ratio_plot,
                       stats_condition_rp107$y.position,
                       stats_dmso_genotype_rp107$y.position,
                       na.rm = TRUE)

set.seed(1)

figure_3d <- ggplot(df_rp107,
                    aes(x = Genotype,
                        y = RpH_ratio_plot,
                        fill = Genotype)) +
  geom_boxplot(aes(group = interaction(Genotype, Condition)),
               width = 0.45,
               outlier.shape = NA,
               colour = "black",
               position = position_dodge(width = dodge_width)) +
  geom_point(aes(shape = Condition,
                 group = Condition),
             position = position_jitterdodge(jitter.width = 0.20,
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 3,
             fill = "black",
             colour = "black",
             stroke = 0) +
  geom_text(data = sample_sizes_rp107,
            aes(x = x,
                y = y,
                label = label),
            inherit.aes = FALSE,
            size = 4) +
  stat_pvalue_manual(stats_condition_rp107,
                     label = "label",
                     xmin = "xmin",
                     xmax = "xmax",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  stat_pvalue_manual(stats_dmso_genotype_rp107,
                     label = "label",
                     xmin = "xmin",
                     xmax = "xmax",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  scale_fill_manual(values = genotype_colours,
                    guide = "none") +
  scale_shape_manual(values = condition_shapes_rp107,
                     name = "Condition") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(0, 1.02 * plot_ymax_rp107),
                  clip = "off") +
  labs(x = "Genotype",
       y = expression("RpH ratio (" * CTF[pHluorin] / CTF[mCherry] * ") " %*% 10^4),
       title = "D: Pharmacological application with RP107") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_3d

# Supplementary Figure S8E: Imaging quality control for the RP107 application
cbr_threshold <- 5
snr_threshold <- 5

df_qc_rp107 <- df_rp107 %>%
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
         qc_pass = CBR_mCherry >= cbr_threshold & SNR_mCherry >= snr_threshold) %>%
  filter(is.finite(CBR_mCherry),
         is.finite(SNR_mCherry),
         is.finite(delta_background))

qc_summary_rp107 <- df_qc_rp107 %>%
  summarise(n = n(),
            pass_rate = mean(qc_pass),
            delta_background_median = median(delta_background),
            delta_background_q1 = quantile(delta_background, 0.25),
            delta_background_q3 = quantile(delta_background, 0.75))

qc_subtitle_rp107 <- with(qc_summary_rp107,
                          sprintf(paste0("CBR \u2265 %g and SNR \u2265 %g: %d%% passed; ",
                                         "\u0394B median %.2f%% (IQR %.2f\u2013%.2f%%)"),
                                  cbr_threshold,
                                  snr_threshold,
                                  round(100 * pass_rate),
                                  delta_background_median,
                                  delta_background_q1,
                                  delta_background_q3))

x_range_rp107 <- range(df_qc_rp107$CBR_mCherry,
                       na.rm = TRUE)

y_range_rp107 <- range(df_qc_rp107$SNR_mCherry,
                       na.rm = TRUE)

x_limits_rp107 <- c(min(floor(x_range_rp107[1]) - 1, cbr_threshold - 1),
                    ceiling(x_range_rp107[2]) + 1)

y_limits_rp107 <- c(max(0,
                        min(floor(y_range_rp107[1]) - 1, snr_threshold - 2)),
                    ceiling(y_range_rp107[2]) + 1)

# Supplementary Figure S8E plot 
condition_shapes_rp107 <- c(DMSO = 21,
                            RP107 = 24)

figure_s7e <- ggplot(df_qc_rp107,
                     aes(x = CBR_mCherry,
                         y = SNR_mCherry,
                         colour = Genotype,
                         shape = Condition,
                         fill = after_scale(colour))) +
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
             alpha = 0.85,
             stroke = 0.2) +
  annotate("text",
           x = cbr_threshold + 0.4,
           y = snr_threshold - 0.8,
           label = "Reference thresholds",
           hjust = 0,
           size = 4) +
  scale_colour_manual(values = genotype_colours,
                      breaks = names(genotype_colours),
                      name = "Genotype") +
  scale_shape_manual(values = condition_shapes_rp107,
                     breaks = names(condition_shapes_rp107),
                     name = "Condition") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = x_limits_rp107,
                  ylim = y_limits_rp107,
                  clip = "off") +
  labs(x = "Contrast-to-background ratio, mCherry",
       y = "Signal-to-noise ratio, mCherry",
       title = "SE: Signal and background adequacy (RP107 application)",
       subtitle = qc_subtitle_rp107) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         fill = "none") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.position = "right",
        legend.title = element_text(face = "bold")) + theme(legend.position="none")

figure_s7e

# Supplementary Figure S8F: Sensitivity checks for the RP107 application
background_threshold <- 0.5
mask_threshold <- 7

df_sensitivity_rp107 <- df_rp107 %>%
  mutate(CTF_mCherry_bg5 = IntDen_mC - Area * mean_back_g5_mC,
         CTF_pHluorin_bg5 = IntDen_pH - Area * mean_back_g5_pH,
         CTF_mCherry_mask1 = IntDen_mC_tight - Area_tight * mean_back_mC,
         CTF_pHluorin_mask1 = IntDen_pH_tight - Area_tight * mean_back_pH) %>%
  mutate(RpH_ratio_bg5 = if_else(CTF_mCherry_bg5 > 0 &  CTF_pHluorin_bg5 > 0,
                                 CTF_pHluorin_bg5 / CTF_mCherry_bg5,
                                 NA_real_),
         RpH_ratio_mask1 = if_else(CTF_mCherry_mask1 > 0 &  CTF_pHluorin_mask1 > 0,
                                   CTF_pHluorin_mask1 / CTF_mCherry_mask1,
                                   NA_real_)) %>%
  mutate(background_change = 100 * abs(RpH_ratio_bg5 - RpH_ratio) / RpH_ratio,
         mask_change = 100 * abs(RpH_ratio_mask1 - RpH_ratio) / RpH_ratio,
         sensitivity_pass = background_change <= background_threshold & mask_change <= mask_threshold) %>%
  filter(RpH_ratio > 0,
         is.finite(background_change),
         is.finite(mask_change))

sensitivity_summary_rp107 <- df_sensitivity_rp107 %>%
  summarise(n = n(),
            pass_rate = mean(sensitivity_pass),
            background_median = median(background_change),
            background_q1 = quantile(background_change, 0.25),
            background_q3 = quantile(background_change, 0.75),
            mask_median = median(mask_change),
            mask_q1 = quantile(mask_change, 0.25),
            mask_q3 = quantile(mask_change, 0.75))

sensitivity_subtitle_rp107 <- with(sensitivity_summary_rp107,
                                   sprintf("Background +5 \u2264 %.1f%% and Mask \u22121 px \u2264 %.0f%%: %d%% passed",
                                           background_threshold,
                                           mask_threshold,
                                           round(100 * pass_rate)))

x_max_sensitivity_rp107 <- 1.10 * max(mask_threshold,
                                      df_sensitivity_rp107$mask_change,
                                      na.rm = TRUE)

y_max_sensitivity_rp107 <- 1.10 * max(background_threshold,
                                      df_sensitivity_rp107$background_change,
                                      na.rm = TRUE)

# Supplementary Figure S8F plot
condition_shapes_rp107 <- c(DMSO = 21,
                            RP107 = 24)

figure_s7f <- ggplot(df_sensitivity_rp107,
                     aes(x = mask_change,
                         y = background_change,
                         colour = Genotype,
                         shape = Condition,
                         fill = after_scale(colour))) +
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
             alpha = 0.90,
             stroke = 0.2) +
  annotate("text",
           x = 0.55 * mask_threshold,
           y = 1.04 * background_threshold,
           label = "Reference thresholds",
           hjust = 0.5,
           vjust = 0,
           size = 4) +
  scale_colour_manual(values = genotype_colours,
                      breaks = names(genotype_colours),
                      name = "Genotype") +
  scale_shape_manual(values = condition_shapes_rp107,
                     breaks = names(condition_shapes_rp107),
                     name = "Condition") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = c(0, x_max_sensitivity_rp107),
                  ylim = c(0, y_max_sensitivity_rp107),
                  clip = "off") +
  labs(x = "Mask \u22121 px, absolute % change",
       y = "Background +5, absolute % change",
       title = "SF: Sensitivity checks (RP107 application)",
       subtitle = sensitivity_subtitle_rp107) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         fill = "none") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_s7f