##############################################################################
# 07_shap_dependence_plots.R
#
# SHAP-based variable interpretation:
#   Fig 6:  Along-wind gradient × coast distance (directional foehn signal)
#   Fig 7:  Upwind-downwind DEM diff × eastness (directional terrain effect)
#   Fig 8:  Terrain complexity × elevation (turbulence threshold)
#   Fig 9:  Broadleaf fraction × season (phenological attenuation)
#   Fig 10: Distance to coast × eastness (foehn amplification zone)
#   Fig 11: SHAP summary plot (variable importance ranking)
#
# Input:  01_intermediate/shap_long.rds
#         01_intermediate/X_train.rds
# Output: 02_output/figures/Fig6~11
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  07. SHAP Dependence & Summary Plots\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Load SHAP data
# =============================================================================
cat("[1/8] Loading SHAP data...\n")

shap_long <- readRDS(inter_path("shap_long.rds"))
X_train   <- readRDS(inter_path("X_train.rds"))

all_vars <- unique(shap_long$variable)
cat(sprintf("  Variables: %d\n", length(all_vars)))

# Identify key column names (exact_extract naming)
find_col <- function(pattern) {
  matched <- grep(pattern, all_vars, value = TRUE)[1]
  if (is.na(matched)) warning(sprintf("Column not found: %s", pattern))
  matched
}

dist_col      <- find_col("mean.*Dist.Coast|mean.*Dist_Coast")
max_dist_col  <- find_col("max.*Dist.Coast|max.*Dist_Coast")
east_col      <- find_col("mean.*Eastness|mean.*eastness")
stdev_dem_col <- find_col("stdev.*DEM|stdev\\.DEM")
mean_dem_col  <- find_col("mean.*DEM|mean\\.DEM")
frtp2_col     <- find_col("FRTP.*2.*frac|2.*FRTP.*frac")
season_sum    <- find_col("season.*Summer|Summer")

# New directional variables
gradient_col  <- find_col("along_wind_gradient")
updown_col    <- find_col("dem_updown_diff")
fetch_col     <- find_col("upwind_fetch")
barrier_col   <- find_col("barrier_index")

cat(sprintf("  Coast dist:  %s\n", dist_col))
cat(sprintf("  Eastness:    %s\n", east_col))
cat(sprintf("  DEM stdev:   %s\n", stdev_dem_col))
cat(sprintf("  DEM mean:    %s\n", mean_dem_col))
cat(sprintf("  Broadleaf:   %s\n", frtp2_col))
cat(sprintf("  Summer:      %s\n", season_sum))

# =============================================================================
# 2. Helper: extract SHAP pair data
# =============================================================================
extract_shap_pair <- function(shap_long, x_var, color_var) {
  x_data <- shap_long %>%
    filter(variable == x_var) %>%
    dplyr::select(ID, shap_x = value, raw_x = rfvalue)

  c_data <- shap_long %>%
    filter(variable == color_var) %>%
    dplyr::select(ID, raw_color = rfvalue)

  inner_join(x_data, c_data, by = "ID")
}

# =============================================================================
# 3. Fig 6: Along-wind gradient × coast distance (SHAP #1, directional)
# =============================================================================
cat("\n[2/8] Fig 6: Along-wind gradient × coast distance...\n")

if (!is.na(gradient_col) && !is.na(dist_col)) {
  df6 <- extract_shap_pair(shap_long, gradient_col, dist_col) %>%
    mutate(
      dist_km = raw_color / 1000,
      coast_group = factor(
        ifelse(dist_km <= 40, "\u2264 40 km (coastal)", "> 40 km (inland)"),
        levels = c("\u2264 40 km (coastal)", "> 40 km (inland)")
      )
    )
  
  fig6 <- ggplot(df6, aes(x = raw_x, y = shap_x)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.5) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "#FF6B6B", alpha = 0.05) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#6B9FFF", alpha = 0.05) +
    geom_point(aes(color = coast_group, shape = coast_group), alpha = 0.5, size = 1.8) +
    geom_smooth(aes(color = coast_group, fill = coast_group),
                method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    scale_color_manual(name = "Coast\nproximity",
                       values = c("\u2264 40 km (coastal)" = "#D62728", "> 40 km (inland)" = "#1F77B4")) +
    scale_fill_manual(name = "Coast\nproximity",
                      values = c("\u2264 40 km (coastal)" = "#D62728", "> 40 km (inland)" = "#1F77B4")) +
    scale_shape_manual(name = "Coast\nproximity",
                       values = c("\u2264 40 km (coastal)" = 17, "> 40 km (inland)" = 16)) +
    annotate("text", x = min(df6$raw_x, na.rm = TRUE) * 0.7,
             y = max(df6$shap_x, na.rm = TRUE) * 0.85,
             label = "Negative gradient\n(upwind higher\n= descending flow)",
             size = 2.8, color = "#D62728", fontface = "italic") +
    annotate("text", x = max(df6$raw_x, na.rm = TRUE) * 0.7,
             y = min(df6$shap_x, na.rm = TRUE) * 0.85,
             label = "Positive gradient\n(upwind lower\n= ascending flow)",
             size = 2.8, color = "#1F77B4", fontface = "italic") +
    labs(
      x = "Along-wind elevation gradient (m/m)",
      y = "SHAP value (contribution to extreme wind speed)"
    )
  
  fig6_margin <- ggExtra::ggMarginal(fig6, type = "histogram", margins = "x",
                                      size = 8, fill = "gray80", color = "gray50")
  
  save_figure(fig6_margin, fig_path("Fig6_SHAP_gradient.tiff"), width = 200, height = 150)
} else {
  cat("  [SKIP] along_wind_gradient not found in SHAP data\n")
}

# =============================================================================
# 4. Fig 7: Upwind-downwind DEM difference × eastness (SHAP #3, directional)
# =============================================================================
cat("[3/8] Fig 7: Upwind-downwind DEM diff × eastness...\n")

if (!is.na(updown_col) && !is.na(east_col)) {
  df7 <- extract_shap_pair(shap_long, updown_col, east_col) %>%
    mutate(
      aspect_group = factor(
        ifelse(raw_color > 0.3, "Leeward (East-facing)", "Other aspects"),
        levels = c("Other aspects", "Leeward (East-facing)")
      )
    )
  
  fig7 <- ggplot(df7, aes(x = raw_x, y = shap_x)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.5) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "#6B9FFF", alpha = 0.05) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#FF6B6B", alpha = 0.05) +
    geom_point(aes(color = aspect_group, shape = aspect_group), alpha = 0.5, size = 1.8) +
    geom_smooth(aes(color = aspect_group, fill = aspect_group),
                method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
    scale_color_manual(name = "Aspect group",
                       values = SHAP_GROUP_COLORS$aspect) +
    scale_fill_manual(name = "Aspect group",
                      values = SHAP_GROUP_COLORS$aspect) +
    scale_shape_manual(name = "Aspect group",
                       values = c("Other aspects" = 16, "Leeward (East-facing)" = 17)) +
    annotate("text", x = min(df7$raw_x, na.rm = TRUE) * 0.6,
             y = max(df7$shap_x, na.rm = TRUE) * 0.85,
             label = "Downwind higher\n(sheltered valley)",
             size = 2.8, color = "#1F77B4", fontface = "italic") +
    annotate("text", x = max(df7$raw_x, na.rm = TRUE) * 0.6,
             y = max(df7$shap_x, na.rm = TRUE) * 0.85,
             label = "Upwind higher\n(mountain barrier\n= foehn descent)",
             size = 2.8, color = "#D62728", fontface = "italic") +
    labs(
      x = "Upwind \u2212 downwind mean elevation (m)",
      y = "SHAP value (contribution to extreme wind speed)"
    )
  
  fig7_margin <- ggExtra::ggMarginal(fig7, type = "histogram", margins = "x",
                                      size = 8, fill = "gray80", color = "gray50")
  
  save_figure(fig7_margin, fig_path("Fig7_SHAP_updown.tiff"), width = 200, height = 150)
} else {
  cat("  [SKIP] dem_updown_diff not found in SHAP data\n")
}

# =============================================================================
# 5. Fig 8: Terrain complexity × mean elevation (Hypothesis 1)
# =============================================================================
cat("[4/8] Fig 8: Terrain complexity × elevation...\n")

df8 <- extract_shap_pair(shap_long, stdev_dem_col, mean_dem_col) %>%
  mutate(
    elev_group = factor(
      ifelse(raw_color >= 600, "\u2265 600 m a.s.l.", "< 600 m a.s.l."),
      levels = c("< 600 m a.s.l.", "\u2265 600 m a.s.l.")
    )
  )

fig8 <- ggplot(df8, aes(x = raw_x, y = shap_x)) +
  annotate("rect", xmin = 200, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#8B4513", alpha = 0.06) +
  geom_vline(xintercept = 200, linetype = "dotted", color = "brown", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_point(aes(color = elev_group, shape = elev_group), alpha = 0.55, size = 1.8) +
  geom_smooth(aes(color = elev_group, fill = elev_group),
              method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
  scale_color_manual(name = "Elevation\ngroup",
                     values = SHAP_GROUP_COLORS$elevation) +
  scale_fill_manual(name = "Elevation\ngroup",
                    values = SHAP_GROUP_COLORS$elevation) +
  scale_shape_manual(name = "Elevation\ngroup",
                     values = c("< 600 m a.s.l." = 16, "\u2265 600 m a.s.l." = 17)) +
  annotate("text", x = 200, y = max(df8$shap_x, na.rm = TRUE) * 0.9,
           label = "Threshold\n(200 m)", size = 2.8, color = "brown",
           hjust = -0.1, fontface = "italic") +
  labs(
    x = "Sub-grid topographic complexity (SD of elevation, m)",
    y = "SHAP value (contribution to extreme wind speed)"
  )

fig8_margin <- ggExtra::ggMarginal(fig8, type = "histogram", margins = "x",
                                    size = 8, fill = "gray80", color = "gray50")

save_figure(fig8_margin, fig_path("Fig8_SHAP_terrain.tiff"), width = 200, height = 150)

# =============================================================================
# 4. Fig 7: Broadleaf fraction × season (Hypothesis 2)
# =============================================================================
cat("[5/8] Fig 9: Broadleaf fraction × season...\n")

df9 <- extract_shap_pair(shap_long, frtp2_col, season_sum) %>%
  mutate(
    season_label = factor(
      ifelse(raw_color == 1, "Leaf-on (Summer)", "Leaf-off (Winter/Spring)"),
      levels = c("Leaf-off (Winter/Spring)", "Leaf-on (Summer)")
    )
  )

fig9 <- ggplot(df9, aes(x = raw_x, y = shap_x)) +
  geom_vline(xintercept = 0.3, linetype = "dotted", color = "darkgreen", linewidth = 0.6) +
  annotate("rect", xmin = 0.3, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#FFA500", alpha = 0.06) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_point(aes(color = season_label, shape = season_label), alpha = 0.55, size = 1.8) +
  geom_smooth(aes(color = season_label, fill = season_label),
              method = "loess", se = TRUE, alpha = 0.12, linewidth = 0.9) +
  scale_color_manual(name = "Phenological\nperiod",
                     values = SHAP_GROUP_COLORS$phenology) +
  scale_fill_manual(name = "Phenological\nperiod",
                    values = SHAP_GROUP_COLORS$phenology) +
  scale_shape_manual(name = "Phenological\nperiod",
                     values = c("Leaf-off (Winter/Spring)" = 17, "Leaf-on (Summer)" = 16)) +
  annotate("text", x = 0.3, y = max(df9$shap_x, na.rm = TRUE) * 0.9,
           label = "Threshold\n(30%)", size = 2.8, color = "darkgreen",
           hjust = -0.1, fontface = "italic") +
  labs(
    x = "Broadleaf forest fraction (area ratio within footprint)",
    y = "SHAP value (contribution to extreme wind speed)"
  )

fig9_margin <- ggExtra::ggMarginal(fig9, type = "histogram", margins = "x",
                                    size = 8, fill = "gray80", color = "gray50")

save_figure(fig9_margin, fig_path("Fig9_SHAP_phenology.tiff"), width = 200, height = 150)

# =============================================================================
# 5. Fig 8: Distance to coast × eastness (Hypothesis 3)
# =============================================================================
cat("[6/8] Fig 10: Distance to coast × eastness...\n")

df10 <- extract_shap_pair(shap_long, dist_col, east_col) %>%
  mutate(
    dist_km = raw_x / 1000,
    aspect_group = factor(
      ifelse(raw_color > 0.3, "Leeward (East-facing)", "Other aspects"),
      levels = c("Other aspects", "Leeward (East-facing)")
    )
  )

fig10 <- ggplot(df10, aes(x = dist_km, y = shap_x)) +
  annotate("rect", xmin = 0, xmax = 25, ymin = -Inf, ymax = Inf,
           fill = "#FF6B6B", alpha = 0.08) +
  geom_vline(xintercept = 25, linetype = "dotted", color = "red3", linewidth = 0.6) +
  geom_vline(xintercept = 50, linetype = "dotted", color = "gray40", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_point(aes(color = aspect_group), alpha = 0.5, size = 1.5) +
  geom_smooth(aes(color = aspect_group, fill = aspect_group),
              method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.9) +
  scale_color_manual(name = "Aspect group",
                     values = SHAP_GROUP_COLORS$aspect) +
  scale_fill_manual(name = "Aspect group",
                    values = SHAP_GROUP_COLORS$aspect) +
  annotate("text", x = 12.5, y = max(df10$shap_x, na.rm = TRUE) * 0.85,
           label = "Foehn\namplification\nzone", size = 3, color = "red3",
           fontface = "italic", hjust = 0.5) +
  annotate("text", x = 25, y = min(df10$shap_x, na.rm = TRUE) * 0.85,
           label = "25 km", size = 2.8, color = "red3", hjust = -0.1) +
  annotate("text", x = 50, y = min(df10$shap_x, na.rm = TRUE) * 0.85,
           label = "50 km", size = 2.8, color = "gray40", hjust = -0.1) +
  labs(
    x = "Distance to coast (km)",
    y = "SHAP value (contribution to extreme wind speed)"
  )

fig10_margin <- ggExtra::ggMarginal(fig10, type = "histogram", margins = "x",
                                    size = 8, fill = "gray80", color = "gray50")

save_figure(fig10_margin, fig_path("Fig10_SHAP_foehn.tiff"), width = 200, height = 150)

# =============================================================================
# 6. Fig 9: SHAP Summary Plot (variable importance)
# =============================================================================
cat("[7/8] Fig 11: SHAP summary plot...\n")

imp_top <- shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
  arrange(desc(mean_abs_shap)) %>%
  slice_head(n = 15)

fig11_data <- shap_long %>%
  filter(variable %in% imp_top$variable) %>%
  mutate(variable = factor(variable, levels = rev(imp_top$variable)))

fig11 <- ggplot(fig11_data, aes(x = value, y = variable, color = rfvalue)) +
  geom_jitter(alpha = 0.4, size = 0.8, height = 0.2) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Feature\nvalue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "SHAP value (impact on prediction)", y = "") +
  theme(axis.text.y = element_text(size = 10))

save_figure(fig11, fig_path("Fig11_SHAP_summary.tiff"), width = 200, height = 180)

# =============================================================================
# 7. Summary
# =============================================================================
cat("\n[8/8] Variable importance ranking (top 10):\n")

imp <- shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = round(mean(abs(value)), 4), .groups = "drop") %>%
  arrange(desc(mean_abs_shap))

print(head(imp, 10))

cat("\n================================================================\n")
cat("  07 Complete\n")
cat("================================================================\n\n")
cat(sprintf("  Fig 6:  %s (along-wind gradient)\n",   fig_path("Fig6_SHAP_gradient.tiff")))
cat(sprintf("  Fig 7:  %s (upwind-downwind diff)\n",  fig_path("Fig7_SHAP_updown.tiff")))
cat(sprintf("  Fig 8:  %s (terrain complexity)\n",     fig_path("Fig8_SHAP_terrain.tiff")))
cat(sprintf("  Fig 9:  %s (phenology)\n",              fig_path("Fig9_SHAP_phenology.tiff")))
cat(sprintf("  Fig 10: %s (foehn)\n",                  fig_path("Fig10_SHAP_foehn.tiff")))
cat(sprintf("  Fig 11: %s (summary)\n",                fig_path("Fig11_SHAP_summary.tiff")))
