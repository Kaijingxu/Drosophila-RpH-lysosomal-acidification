library(dplyr)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(readxl)
library(sandwich)

# Import and prepare data
df_raw <- readxl::read_excel("data/cg_RpH/BafA1.xlsx")

scale_factor <- 1e4
dodge_width <- 0.60

df_baf <- df_raw %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("+/+", "s64w/+", "s64w/s64w")),
         Condition = factor(Condition,
                            levels = c("DMSO", "BafA1")),
         CTF_mCherry = IntDen_mC - Area * mean_back_mC,
         CTF_pHluorin = IntDen_pH - Area * mean_back_pH) %>%
  filter(Area > 0,
         CTF_mCherry > 0,
         CTF_pHluorin > 0) %>%
  mutate(RpH_ratio = CTF_pHluorin / CTF_mCherry,
         RpH_ratio_plot = scale_factor * RpH_ratio,
         log_RpH_ratio = log(RpH_ratio))

# Statistical analysis
#linear model
fit_baf <- lm(log_RpH_ratio ~ Genotype * Condition,
              data = df_baf)
vcov_hc3 <- sandwich::vcovHC(fit_baf, type = "HC3")

# BafA1 versus DMSO within each genotype
emm_condition <- emmeans(fit_baf,
                         ~ Condition | Genotype,
                         vcov. = vcov_hc3)

condition_contrasts <- contrast(emm_condition,
                                method = "revpairwise",
                                adjust = "holm") %>%
  summary(infer = TRUE) %>%
  as.data.frame()

# Genotype comparisons under DMSO
emm_genotype <- emmeans(fit_baf,
                        ~ Genotype | Condition,
                        vcov. = vcov_hc3)

dmso_genotype_contrasts <- contrast(emm_genotype,
                                    method = "pairwise",
                                    adjust = "holm") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  filter(Condition == "DMSO")

# Format statistical annotations
format_p <- function(p) {case_when(is.na(p)   ~ NA_character_,
                                   p < 1e-4   ~ paste0("p = ", format(signif(p, 3), scientific = TRUE)),
                                   p < 0.01   ~ paste0("p = ", formatC(p, format = "f", digits = 4)),
                                   p < 0.10   ~ paste0("p = ", formatC(p, format = "f", digits = 3)),
                                   TRUE       ~ paste0("p = ", formatC(p, format = "f", digits = 2)))}

genotype_levels <- levels(df_baf$Genotype)

x_centres <- setNames(seq_along(genotype_levels),
                      genotype_levels)

x_offsets <- c(DMSO = -dodge_width / 4,
               BafA1 = dodge_width / 4)

y_max <- max(df_baf$RpH_ratio_plot, na.rm = TRUE)

genotype_ymax <- df_baf %>%
  group_by(Genotype) %>%
  summarise(ymax = max(RpH_ratio_plot, na.rm = TRUE),
            .groups = "drop")

# BafA1 versus DMSO annotations
stats_condition <- condition_contrasts %>%
  transmute(Genotype = as.character(Genotype),
            group1 = "DMSO",
            group2 = "BafA1",
            label = format_p(p.value),
            xmin = unname(x_centres[Genotype]) + x_offsets["DMSO"],
            xmax = unname(x_centres[Genotype]) + x_offsets["BafA1"]) %>%
  left_join(genotype_ymax, by = "Genotype") %>%
  mutate(y.position = 1.12 * ymax)

# DMSO genotype annotations
dmsо_pairs <- tibble(group1 = c("+/+", "+/+", "s64w/+"),
                     group2 = c("s64w/+", "s64w/s64w", "s64w/s64w"))

contrast_pairs <- dmso_genotype_contrasts %>%
  tidyr::separate(contrast,
                  into = c("group1", "group2"),
                  sep = " - ",
                  remove = FALSE) %>%
  mutate(group1 = trimws(gsub("[`()]", "", group1)),
         group2 = trimws(gsub("[`()]", "", group2)))

stats_dmsо_genotype <- dmsо_pairs %>%
  left_join(contrast_pairs %>%
              select(group1, group2, p.value),
            by = c("group1", "group2")) %>%
  mutate(label = format_p(p.value),
         xmin = x_centres[group1] + x_offsets["DMSO"],
         xmax = x_centres[group2] + x_offsets["DMSO"],
         y.position = c(1.12, 1.20, 1.28) * y_max)

# Sample-size annotations
sample_sizes <- df_baf %>%
  count(Genotype, Condition, name = "n") %>%
  mutate(x = x_centres[as.character(Genotype)] + x_offsets[as.character(Condition)],
         y = 0.05 * y_max,
         label = paste0("n = ", n))

# Figure 3B
genotype_colours <- c("+/+" = "blue",
                      "s64w/+" = "green",
                      "s64w/s64w" = "red")

set.seed(1)

figure_3b <- ggplot(df_baf,
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
  geom_text(data = sample_sizes,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 4) +
  stat_pvalue_manual(stats_condition,
                     label = "label",
                     xmin = "xmin",
                     xmax = "xmax",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  stat_pvalue_manual(stats_dmsо_genotype,
                     label = "label",
                     xmin = "xmin",
                     xmax = "xmax",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  scale_fill_manual(values = genotype_colours,
                    guide = "none") +
  scale_shape_manual(values = c(DMSO = 21,
                                BafA1 = 25),
                     name = "Condition") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(0, 1.30 * y_max),
                  clip = "off") +
  labs(x = "Genotype",
       y = expression("RpH ratio (" * CTF[pHluorin] / CTF[mCherry] * ") " %*%  10^4),
       title = "B: Validation with BafA1") +
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

figure_3b

# Supplementary Figure S8A: cg-RpH fat-body imaging quality control
cbr_threshold <- 5
snr_threshold <- 5
background_threshold <- 1

df_qc_fatbody <- df_baf %>%
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

qc_summary_fatbody <- df_qc_fatbody %>%
  summarise(n = n(),
            pass_rate = mean(qc_pass),
            delta_background_median = median(delta_background),
            delta_background_q1 = quantile(delta_background, 0.25),
            delta_background_q3 = quantile(delta_background, 0.75),
            background_pass_rate = mean(delta_background <= background_threshold))

qc_subtitle <- with(qc_summary_fatbody,
                    sprintf(paste0("CBR \u2265 %g and SNR \u2265 %g: %d%% passed; ",
                                   "\u0394B median %.2f%% (IQR %.2f\u2013%.2f%%)"),
                            cbr_threshold,
                            snr_threshold,
                            round(100 * pass_rate),
                            delta_background_median,
                            delta_background_q1,
                            delta_background_q3))

# Supplementary Figure S8A plot
condition_shapes <- c(DMSO = 21,
                      BafA1 = 25)

x_range <- range(df_qc_fatbody$CBR_mCherry, na.rm = TRUE)
y_range <- range(df_qc_fatbody$SNR_mCherry, na.rm = TRUE)

x_limits <- c(min(floor(x_range[1]) - 1, cbr_threshold - 1),
              ceiling(x_range[2]) + 1)

y_limits <- c(max(0, min(floor(y_range[1]) - 1, snr_threshold - 2)),
              ceiling(y_range[2]) + 1)

figure_s7a <- ggplot(df_qc_fatbody,
                       aes(x = CBR_mCherry,
                           y = SNR_mCherry,
                           colour = Genotype,
                           fill = after_scale(colour),
                           shape = Condition)) +
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
  scale_fill_manual(values = genotype_colours,
                    breaks = names(genotype_colours),
                    guide = "none") +
  scale_shape_manual(values = condition_shapes,
                     breaks = names(condition_shapes),
                     name = "Condition") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = x_limits,
                  ylim = y_limits,
                  clip = "off") +
  labs(x = "Contrast-to-background ratio, mCherry",
       y = "Signal-to-noise ratio, mCherry",
       title = "SA: Signal and background adequacy (BafA1 validation)",
       subtitle = qc_subtitle) +
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

figure_s7a

# Supplementary Figure S8B: cg-RpH fat-body imaging Sensitivity to background and mask definitions
background_threshold <- 0.5
mask_threshold <- 7

df_sensitivity_baf <- df_baf %>%
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
         mask_change = 100 * abs(RpH_ratio_mask1 - RpH_ratio) / RpH_ratio,
         sensitivity_pass = (background_change <= background_threshold &  mask_change <= mask_threshold)) %>%
  filter(RpH_ratio > 0,
         is.finite(background_change),
         is.finite(mask_change))

sensitivity_summary_baf <- df_sensitivity_baf %>%
  summarise(n = n(),
            pass_rate = mean(sensitivity_pass),
            background_median = median(background_change),
            background_q1 = quantile(background_change, 0.25),
            background_q3 = quantile(background_change, 0.75),
            mask_median = median(mask_change),
            mask_q1 = quantile(mask_change, 0.25),
            mask_q3 = quantile(mask_change, 0.75))

sensitivity_subtitle_baf <- with(sensitivity_summary_baf,
                                 sprintf("Background +5 \u2264 %.1f%% and Mask \u22121 px \u2264 %.0f%%: %d%% passed",
                                         background_threshold,
                                         mask_threshold,
                                         round(100 * pass_rate)))

# Supplementary Figure S8B plot
x_range_sensitivity <- range(df_sensitivity_baf$mask_change,
                             na.rm = TRUE)

y_range_sensitivity <- range(df_sensitivity_baf$background_change,
                             na.rm = TRUE)

x_limits_sensitivity <- c(0,
                          1.10 * max(mask_threshold, x_range_sensitivity[2]))

y_limits_sensitivity <- c(0,
                          1.10 * max(background_threshold, y_range_sensitivity[2]))

figure_s7b <- ggplot(df_sensitivity_baf,
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
  scale_shape_manual(values = condition_shapes,
                     breaks = names(condition_shapes),
                     name = "Condition") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = x_limits_sensitivity,
                  ylim = y_limits_sensitivity,
                  clip = "off") +
  labs(x = "Mask \u22121 px, absolute % change",
       y = "Background +5, absolute % change",
       title = "SB: Sensitivity checks (BafA1 validation)",
       subtitle = sensitivity_subtitle_baf) +
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

figure_s7b