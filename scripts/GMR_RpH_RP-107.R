library(readxl)
library(dplyr)
library(ggplot2)
library(emmeans)
library(sandwich)
library(ggpubr)

# Import data
df <- readxl::read_excel("data/GMR_RpH/RP107.xlsx")

# Settings
scale_factor <- 100
eps <- 1e-12
dodge_width <- 0.60

pd <- position_dodge(width = dodge_width)

pos_jit <- position_jitterdodge(jitter.width = 0.20,
                                jitter.height = 0,
                                dodge.width = dodge_width)

genotype_levels <- c("+/+",
                     "s64w/+",
                     "s64w/s64w")

condition_levels <- c("DMSO",
                      "RP107")

cols <- c("+/+" = "blue",
          "s64w/+" = "green",
          "s64w/s64w" = "red")

x_off <- c("DMSO" = -dodge_width / 4,
           "RP107" =  dodge_width / 4)

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
df1 <- df %>%
  mutate(Genotype = factor(Genotype,
                           levels = genotype_levels),
         Condition = factor(Condition,
                            levels = condition_levels)) %>%
  filter(puncta_area > 0,
         IntDen_corr_mC > 0,
         IntDen_corr_pH > 0) %>%
  mutate(RpH_ratio     = IntDen_corr_pH / IntDen_corr_mC,
         ratio_plot    = scale_factor * RpH_ratio,
         log_RpH_ratio = log(RpH_ratio + eps)) %>%
  filter(is.finite(RpH_ratio),
         is.finite(log_RpH_ratio))

# Linear model with HC3 robust covariance
fit <- lm(log_RpH_ratio ~ Genotype * Condition,
                    data = df1)

Vhc3 <- sandwich::vcovHC(fit,
                         type = "HC3")

# RP107 versus DMSO within each genotype
emm_treat <- emmeans(fit,
                     ~ Condition | Genotype,
                     vcov. = Vhc3)

p_within <- contrast(emm_treat,
                     method = "revpairwise",
                     adjust = "holm") %>%
  summary() %>%
  as.data.frame()

# Genotype comparisons within DMSO
emm_geno <- emmeans(fit,
                    ~ Genotype | Condition,
                    vcov. = Vhc3)

p_geno <- contrast(emm_geno,
                   method = "pairwise",
                   adjust = "holm") %>%
  summary() %>%
  as.data.frame() %>%
  filter(Condition == "DMSO")

# Plot positions
geno_lv <- levels(df1$Genotype)

x_center <- setNames(seq_along(geno_lv),
                     geno_lv)

y_all <- max(df1$ratio_plot,
             na.rm = TRUE)

y_byG <- df1 %>%
  group_by(Genotype) %>%
  summarise(ymax = max(ratio_plot,
                       na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Genotype = as.character(Genotype))

# DMSO versus RP107 within each genotype
stats_within <- p_within %>%
  mutate(Genotype = as.character(Genotype),
         group1 = "DMSO",
         group2 = "RP107",
         xmin = unname(x_center[Genotype] + x_off["DMSO"]),
         xmax = unname(x_center[Genotype] + x_off["RP107"]),
         label = fmt_p(p.value)) %>%
  left_join(y_byG,
            by = "Genotype") %>%
  mutate(y.position = ymax + 0.055 * y_all) %>%
  filter(Genotype == "s64w/s64w")

# DMSO genotype comparisons
split_con <- strsplit(as.character(p_geno$contrast),
                      " - ",
                      fixed = TRUE)

g1 <- vapply(split_con,
             `[`,
             character(1),
             1)

g2 <- vapply(split_con,
             `[`,
             character(1),
             2)

g1 <- trimws(gsub("[`()]", "",
                  g1))

g2 <- trimws(gsub("[`()]", "",
                  g2))

target <- tibble::tibble(group1 = c("+/+",
                                    "+/+",
                                    "s64w/+"),
                         group2 = c("s64w/+",
                                    "s64w/s64w",
                                    "s64w/s64w"))

match_idx <- match(paste(target$group1,
                         target$group2),
                   paste(g1,
                         g2))

stats_dmsogeno <- target %>%
  mutate(label = fmt_p(p_geno$p.value[match_idx]),
         xmin = unname(x_center[group1] + x_off["DMSO"]),
         xmax = unname(x_center[group2] + x_off["DMSO"]),
         y.position = c(0.95,
                        1.04,
                        1.13) * y_all)

# Sample-size labels
nlab <- df1 %>%
  count(Genotype,
        Condition,
        name = "n") %>%
  mutate(x = unname(x_center[as.character(Genotype)] + x_off[as.character(Condition)]),
         y = 0.13 * y_all,
         label = paste0("n=",
                        n))

# Y-axis limit
annotation_max <- max(stats_within$y.position,
                      stats_dmsogeno$y.position,
                      na.rm = TRUE)

y_limit <- max(100,
               annotation_max * 1.05)

# Figure 4C
set.seed(1)

figure_4c <- ggplot(df1,
                    aes(x = Genotype,
                        y = ratio_plot,
                        fill = Genotype)) +
  geom_boxplot(aes(group = interaction(Genotype,
                                       Condition)),
               width = 0.45,
               outlier.shape = NA,
               colour = "black",
               coef = Inf,
               position = pd) +
  geom_point(aes(shape = Condition,
                 group = Condition),
             position = pos_jit,
             size = 3,
             fill = "black",
             colour = "black",
             stroke = 0) +
  geom_text(data = nlab,
            aes(x = x,
                y = y,
                label = label),
            inherit.aes = FALSE,
            size = 4) +
  ggpubr::stat_pvalue_manual(stats_within,
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             tip.length = 0.01,
                             size = 4.5,
                             inherit.aes = FALSE) +
  ggpubr::stat_pvalue_manual(stats_dmsogeno,
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             tip.length = 0.01,
                             size = 4.5,
                             inherit.aes = FALSE) +
  scale_fill_manual(values = cols,
                    guide = "none") +  
  scale_shape_manual(values = c("DMSO" = 21,
                                "RP107" = 24),
                     name = "Condition") +
  scale_y_continuous(limits = c(0,
                                y_limit),
                     expand = expansion(mult = c(0.02,
                                                 0.04))) +
  labs(x = "Genotype",
       y = expression("RpH ratio (" * CTF[pHluorin] / CTF[mCherry] * ") × 100"),
       title = "D: Pharmacological application with RP107") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(fill = "black",
                                                  colour = "black",
                                                  size = 3))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_4c

# Supplementary Figure S9B: Background CV and mask coverage: RP107 application
cv_threshold       <- 0.30
coverage_threshold <- 25
eps                 <- 1e-12

condition_shapes <- c("DMSO"  = 16,
                      "RP107" = 17)

eye_rp107_qc <- df1 %>%
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
cv_pass_rate <- mean(eye_rp107_qc$CV_mCherry <= cv_threshold & eye_rp107_qc$CV_pHluorin <= cv_threshold)

coverage_pass_rate <- mean(eye_rp107_qc$mask_coverage >= coverage_threshold)

coverage_median <- median(eye_rp107_qc$mask_coverage,
                          na.rm = TRUE)

coverage_iqr <- quantile(eye_rp107_qc$mask_coverage,
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
x_upper <- max(eye_rp107_qc$CV_mCherry,
               cv_threshold,
               na.rm = TRUE) * 1.15

y_upper <- max(eye_rp107_qc$CV_pHluorin,
               cv_threshold,
               na.rm = TRUE) * 1.15

# Supplementary Figure S9B
figure_s9b <- ggplot(eye_rp107_qc,
                     aes(x = CV_mCherry,
                         y = CV_pHluorin,
                         colour = Genotype,
                         shape = Condition)) +
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
  geom_point(size = 2.8,
             alpha = 0.90,
             stroke = 0.3) +
  annotate("text",
           x = 0.52 * cv_threshold,
           y = 1.06 * cv_threshold,
           label = "Reference thresholds",
           hjust = 0,
           size = 4.5) +
  scale_colour_manual(values = cols,
                      breaks = names(cols),
                      name = "Genotype") +
  scale_shape_manual(values = condition_shapes,
                     breaks = names(condition_shapes),
                     name = "Condition") +
  scale_x_continuous(limits = c(0, x_upper),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(limits = c(0, y_upper),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(title = "SB: Background CV and mask coverage (RP107 application)",
       subtitle = qc_subtitle,
       x = "Background CV, mCherry",
       y = "Background CV, pHluorin") +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

figure_s9b