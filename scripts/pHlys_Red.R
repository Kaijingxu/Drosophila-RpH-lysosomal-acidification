library(tidyverse)
library(glmmTMB)
library(emmeans)
library(ggpubr)
library(readr)

# Import and prepare replicate 1
df_r1 <- readxl::read_excel("data/pHlys_Red/pHlys red 1.xlsx") %>%
  transmute(replicate = "R1",
            fly_id = `Fly number`,
            genotype = factor(Genotype,
                              levels = c("+/+", "s64w/+", "s64w/s64w")),
            gfp = GFP,
            overlap = AND,
            phlys_red = pHlysred,
            # Panel B: percentage of Lamp1-GFP area positive for pHlys Red
            pct_acidified = 100 * overlap / gfp,
            # Diagnostic variables
            relative_red_load = 100 * phlys_red / gfp,
            pct_red_in_gfp = 100 * overlap / phlys_red) %>%
  filter(gfp > 0,
         phlys_red > 0,
         is.finite(pct_acidified),
         is.finite(relative_red_load),
         is.finite(pct_red_in_gfp))

# Import and prepare replicate 2
df_r2 <- readxl::read_excel("data/pHlys_Red/pHlys red 2.xlsx") %>%
  transmute(replicate = "R2",
            fly_id = `Fly number`,
            genotype = factor(Genotype,
                              levels = c("+/+", "s64w/+", "s64w/s64w")),
            gfp = GFP,
            overlap = AND,
            phlys_red = pHlysred,
            # Panel B: percentage of Lamp1-GFP area positive for pHlys Red
            pct_acidified = 100 * overlap / gfp,
            # Diagnostic variables
            relative_red_load = 100 * phlys_red / gfp,
            pct_red_in_gfp = 100 * overlap / phlys_red) %>%
  filter(gfp > 0,
         phlys_red > 0,
         is.finite(pct_acidified),
         is.finite(relative_red_load),
         is.finite(pct_red_in_gfp))

# Beta regression and pairwise comparisons
# Pairwise contrast definitions shared by R1 and R2
contrast_definitions <- list("+/+ vs s64w/+"       = c(1, -1, 0),
                             "+/+ vs s64w/s64w"    = c(1, 0, -1),
                             "s64w/+ vs s64w/s64w" = c(0, 1, -1))

# Format adjusted p values for plotting
format_p <- function(p) {case_when(is.na(p)   ~ NA_character_,
                                   p < 0.0001 ~ paste0("p = ",
                                                       format(p, scientific = TRUE, digits = 3)),
                                   p < 0.01   ~ sprintf("p = %.4f", p),
                                   p < 0.10   ~ sprintf("p = %.3f", p),
                                   TRUE       ~ sprintf("p = %.2f", p))}

# Replicate 1
# Convert proportions from [0, 1] to (0, 1), because beta
# regression does not permit exact zero or one values
n_r1 <- nrow(df_r1)

model_data_r1 <- df_r1 %>%
  mutate(proportion = pct_acidified / 100,
         proportion_adjusted = (proportion * (n_r1 - 1) + 0.5) / n_r1)

# Beta regression with logit link
fit_beta_r1 <- glmmTMB(proportion_adjusted ~ genotype,
                       family = beta_family(link = "logit"),
                       data = model_data_r1)

# Estimated marginal means
emm_r1 <- emmeans(fit_beta_r1,
                  ~ genotype,
                  type = "response")

# Pairwise comparisons with Holm correction within R1
contrasts_r1 <- contrast(emm_r1,
                         method = contrast_definitions,
                         adjust = "holm") %>%
  summary() %>%
  as.data.frame()

plot_stat_r1 <- contrasts_r1 %>%
  mutate(replicate = "R1",
         group1 = case_when(contrast == "+/+ vs s64w/+" ~ "+/+",
                            contrast == "+/+ vs s64w/s64w" ~ "+/+",
                            contrast == "s64w/+ vs s64w/s64w" ~ "s64w/+"),
         group2 = case_when(contrast == "+/+ vs s64w/+" ~ "s64w/+",
                            contrast == "+/+ vs s64w/s64w" ~ "s64w/s64w",
                            contrast == "s64w/+ vs s64w/s64w" ~ "s64w/s64w"),
         label = format_p(p.value),
         y.position = case_when(contrast == "+/+ vs s64w/+" ~ 94,
                                contrast == "s64w/+ vs s64w/s64w" ~ 114,
                                contrast == "+/+ vs s64w/s64w" ~ 104)) %>%
  select(replicate,
         group1,
         group2,
         p.adj = p.value,
         label,
         y.position)

plot_stat_r1

# Sample sizes
sample_sizes_r1 <- df_r1 %>%
  count(genotype, name = "n") %>%
  mutate(replicate = "R1",
         y.position = 2,
         label = paste0("n = ", n))

sample_sizes_r1

# Replicate 2
n_r2 <- nrow(df_r2)

model_data_r2 <- df_r2 %>%
  mutate(proportion = pct_acidified / 100,
         proportion_adjusted = (proportion * (n_r2 - 1) + 0.5) / n_r2)

# Beta regression with logit link
fit_beta_r2 <- glmmTMB(proportion_adjusted ~ genotype,
                       family = beta_family(link = "logit"),
                       data = model_data_r2)

# Estimated marginal means
emm_r2 <- emmeans(fit_beta_r2,
                  ~ genotype,
                  type = "response")

# Pairwise comparisons with Holm correction within R2
contrasts_r2 <- contrast(emm_r2,
                         method = contrast_definitions,
                         adjust = "holm") %>%
  summary() %>%
  as.data.frame()

plot_stat_r2 <- contrasts_r2 %>%
  mutate(replicate = "R2",
         group1 = case_when(contrast == "+/+ vs s64w/+" ~ "+/+",
                            contrast == "+/+ vs s64w/s64w" ~ "+/+",
                            contrast == "s64w/+ vs s64w/s64w" ~ "s64w/+"),
         group2 = case_when(contrast == "+/+ vs s64w/+" ~ "s64w/+",
                            contrast == "+/+ vs s64w/s64w" ~ "s64w/s64w",
                            contrast == "s64w/+ vs s64w/s64w" ~ "s64w/s64w"),
         label = format_p(p.value),
         y.position = case_when(contrast == "+/+ vs s64w/+" ~ 94,
                                contrast == "s64w/+ vs s64w/s64w" ~ 114,
                                contrast == "+/+ vs s64w/s64w" ~ 104)) %>%
  select(replicate,
         group1,
         group2,
         p.adj = p.value,
         label,
         y.position)

plot_stat_r2

# Sample sizes
sample_sizes_r2 <- df_r2 %>%
  count(genotype, name = "n") %>%
  mutate(replicate = "R2",
         y.position = 2,
         label = paste0("n = ", n))

sample_sizes_r2

# Combine R1 and R2 for plotting
df_plot <- bind_rows(df_r1,
                     df_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

plot_stats <- bind_rows(plot_stat_r1,
                        plot_stat_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

sample_sizes <- bind_rows(sample_sizes_r1,
                          sample_sizes_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

# Generate Figure 2B
genotype_colours <- c("+/+" = "blue",
                      "s64w/+" = "green",
                      "s64w/s64w" = "red")

set.seed(1)

figure_2b <- ggplot(df_plot,
                    aes(x = genotype, y = pct_acidified, fill = genotype)) +
  geom_boxplot(width = 0.45,
               outlier.shape = NA,
               colour = "black") +
  geom_jitter(width = 0.20,
              size = 2,
              alpha = 1,
              show.legend = FALSE) +
  geom_text(data = sample_sizes,
            aes(x = genotype,
                y = y.position,
                label = label),
            inherit.aes = FALSE,
            size = 5) +
  stat_pvalue_manual(plot_stats,
                     label = "label",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 5,
                     inherit.aes = FALSE) +
  facet_wrap(~ replicate,
             nrow = 1) +
  scale_fill_manual(values = genotype_colours,
                    name = "Genotype") +
  scale_y_continuous(breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0.02, 0.03))) +
  coord_cartesian(ylim = c(0, 118),
                  clip = "off") +
  labs(x = "Genotype",
       y = "Lamp1-GFP area positive for pHlys Red (%)",
       title = "B: Quantification of pHlys Red localisation") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor    = element_blank(),
        plot.title.position = "plot",
        plot.title          = element_text(hjust = 0, face = "bold"),
        legend.position     = "right",
        legend.box          = "vertical",
        legend.title        = element_text(face = "bold"))

figure_2b

# Supplementary Figure S2A
# Replicate 1
gfp_r1 <- df_r1 %>%
  filter(is.finite(gfp),
         gfp > 0)

kw_gfp_r1 <- kruskal.test(gfp ~ genotype,
                          data = gfp_r1)

pairwise_gfp_r1 <- pairwise.wilcox.test(x = gfp_r1$gfp,
                                        g = gfp_r1$genotype,
                                        p.adjust.method = "holm",
                                        exact = FALSE)

# Replicate 2
gfp_r2 <- df_r2 %>%
  filter(is.finite(gfp),
         gfp > 0)

kw_gfp_r2 <- kruskal.test(gfp ~ genotype,
                          data = gfp_r2)

pairwise_gfp_r2 <- pairwise.wilcox.test(x = gfp_r2$gfp,
                                        g = gfp_r2$genotype,
                                        p.adjust.method = "holm",
                                        exact = FALSE)

gfp_comparisons <- list(c("+/+", "s64w/+"),
                        c("+/+", "s64w/s64w"),
                        c("s64w/+", "s64w/s64w"))

run_pairwise_wilcox <- function(data, replicate_label) {
  purrr::map_dfr(gfp_comparisons,
                 function(groups) {
                   test_data <- data %>%
                     filter(genotype %in% groups) %>%
                     droplevels()
                   test_result <- wilcox.test(gfp ~ genotype,
                                              data = test_data,
                                              exact = FALSE)
                   tibble(replicate = replicate_label,
                          group1 = groups[1],
                          group2 = groups[2],
                          p.raw = test_result$p.value)}) %>%
    mutate(p.adj = p.adjust(p.raw, method = "holm"),
           label = format_p(p.adj))}

gfp_stats_r1 <- run_pairwise_wilcox(gfp_r1,
                                    "R1")

gfp_stats_r2 <- run_pairwise_wilcox(gfp_r2,
                                    "R2")

# Set bracket positions on the original area scale
gfp_max_r1 <- max(gfp_r1$gfp, na.rm = TRUE)
gfp_max_r2 <- max(gfp_r2$gfp, na.rm = TRUE)

gfp_range_r1 <- diff(range(gfp_r1$gfp, na.rm = TRUE))
gfp_range_r2 <- diff(range(gfp_r2$gfp, na.rm = TRUE))

gfp_stats_r1 <- gfp_stats_r1 %>%
  mutate(y.position = case_when(group1 == "+/+" & group2 == "s64w/+" ~  gfp_max_r1 + 0.10 * gfp_range_r1,
                                group1 == "s64w/+" & group2 == "s64w/s64w" ~ gfp_max_r1 + 0.34 * gfp_range_r1,
                                group1 == "+/+" & group2 == "s64w/s64w" ~ gfp_max_r1 + 0.22 * gfp_range_r1))

gfp_stats_r2 <- gfp_stats_r2 %>%
  mutate(y.position = case_when(group1 == "+/+" & group2 == "s64w/+" ~ gfp_max_r2 + 0.10 * gfp_range_r2,
                                group1 == "s64w/+" & group2 == "s64w/s64w" ~  gfp_max_r2 + 0.34 * gfp_range_r2,
                                group1 == "+/+" & group2 == "s64w/s64w" ~ gfp_max_r2 + 0.22 * gfp_range_r2))

# Combine R1 and R2 for plotting
gfp_plot_data <- bind_rows(gfp_r1,
                           gfp_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

# Keep and display all three pairwise comparisons
gfp_plot_stats <- bind_rows(gfp_stats_r1,
                            gfp_stats_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

# Supplementary Figure S2A
genotype_colours <- c("+/+" = "blue",
                      "s64w/+" = "green",
                      "s64w/s64w" = "red")

set.seed(2)

figure_s6a <- ggplot(gfp_plot_data,
                     aes(x = genotype,
                         y = gfp,
                         fill = genotype)) +
  geom_boxplot(width = 0.45,
               outlier.shape = NA,
               colour = "black") +
  geom_jitter(width = 0.20,
              size = 2,
              alpha = 0.9,
              show.legend = FALSE) +
  ggpubr::stat_pvalue_manual(gfp_plot_stats,
                             label = "label",
                             xmin = "group1",
                             xmax = "group2",
                             y.position = "y.position",
                             tip.length = 0.01,
                             size = 4.5,
                             inherit.aes = FALSE) +
  facet_wrap(~ replicate,
             nrow = 1) +
  scale_fill_manual(values = genotype_colours,
                    name = "Genotype") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "Genotype",
       y = expression("Total Lamp1-GFP area (" * mu * "m"^2 * ")"),
       title = "SA: Total Lamp1-GFP area") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0,
                                  face = "bold"),
        strip.background = element_rect(fill = "grey85",
                                        colour = "black"),
        strip.text = element_text(face = "plain"),
        legend.position = "right",
        legend.title = element_text(face = "bold")) + theme(legend.position="none")

figure_s6a

# Supplementary Figure S2B
# Replicate 1
overlap_r1 <- df_r1 %>%
  filter(is.finite(overlap),
         overlap >= 0)

kw_overlap_r1 <- kruskal.test(overlap ~ genotype,
                              data = overlap_r1)

# Replicate 2
overlap_r2 <- df_r2 %>%
  filter(is.finite(overlap),
         overlap >= 0)

kw_overlap_r2 <- kruskal.test(overlap ~ genotype,
                              data = overlap_r2)

overlap_comparisons <- list(c("+/+", "s64w/+"),
                            c("+/+", "s64w/s64w"),
                            c("s64w/+", "s64w/s64w"))

run_pairwise_overlap <- function(data, replicate_label) {
  purrr::map_dfr(overlap_comparisons,
                 function(groups) {test_data <- data %>%
                   filter(genotype %in% groups) %>%
                   droplevels()
                 test_result <- wilcox.test(overlap ~ genotype,
                                            data = test_data,
                                            exact = FALSE)
                 tibble(replicate = replicate_label,
                        group1 = groups[1],
                        group2 = groups[2],
                        p.raw = test_result$p.value)}) %>%
    mutate(p.adj = p.adjust(p.raw, method = "holm"),
           label = format_p(p.adj))}

overlap_stats_r1 <- run_pairwise_overlap(overlap_r1,
                                         "R1")

overlap_stats_r2 <- run_pairwise_overlap(overlap_r2,
                                         "R2")

# Set bracket positions on the original area scale
overlap_max_r1 <- max(overlap_r1$overlap, na.rm = TRUE)
overlap_max_r2 <- max(overlap_r2$overlap, na.rm = TRUE)

overlap_range_r1 <- diff(range(overlap_r1$overlap, na.rm = TRUE))
overlap_range_r2 <- diff(range(overlap_r2$overlap, na.rm = TRUE))

overlap_stats_r1 <- overlap_stats_r1 %>%
  mutate(y.position = case_when(group1 == "+/+" & group2 == "s64w/+" ~ overlap_max_r1 + 0.10 * overlap_range_r1,
                                group1 == "+/+" & group2 == "s64w/s64w" ~ overlap_max_r1 + 0.22 * overlap_range_r1,
                                group1 == "s64w/+" & group2 == "s64w/s64w" ~ overlap_max_r1 + 0.34 * overlap_range_r1))

overlap_stats_r2 <- overlap_stats_r2 %>%
  mutate(y.position = case_when(group1 == "+/+" & group2 == "s64w/+" ~ overlap_max_r2 + 0.10 * overlap_range_r2,
                                group1 == "+/+" & group2 == "s64w/s64w" ~ overlap_max_r2 + 0.22 * overlap_range_r2,
                                group1 == "s64w/+" & group2 == "s64w/s64w" ~ overlap_max_r2 + 0.34 * overlap_range_r2))

# Combine R1 and R2 for plotting
overlap_plot_data <- bind_rows(overlap_r1,
                               overlap_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

# Keep and display all three pairwise comparisons
overlap_plot_stats <- bind_rows(overlap_stats_r1,
                                overlap_stats_r2) %>%
  mutate(replicate = factor(replicate,
                            levels = c("R1", "R2")))

# Supplementary Figure S2B
set.seed(3)

figure_s6b <- ggplot(overlap_plot_data,
                     aes(x = genotype,
                         y = overlap,
                         fill = genotype)) +
  geom_boxplot(width = 0.45,
               outlier.shape = NA,
               colour = "black") +
  geom_jitter(width = 0.20,
              size = 2,
              alpha = 0.9,
              show.legend = FALSE) +
  ggpubr::stat_pvalue_manual(overlap_plot_stats,
                             label = "label",
                             xmin = "group1",
                             xmax = "group2",
                             y.position = "y.position",
                             tip.length = 0.01,
                             size = 4.5,
                             inherit.aes = FALSE) +
  facet_wrap(~ replicate,
             nrow = 1) +
  scale_fill_manual(values = genotype_colours,
                    name = "Genotype") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "Genotype",
       y = expression("Total overlap area (" * mu * "m"^2 * ")"),
       title = "SB: Total Lamp1-GFP–pHlys Red overlap area") +
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

figure_s6b