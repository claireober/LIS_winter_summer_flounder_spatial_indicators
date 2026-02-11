

# LIS Change Point Analysis (CT DEEP trawl survey)
# Supporting code for:
# "Evaluating spatial and temporal dynamics of cold- and warm-adapted fish in a changing Long Island Sound"
#
# This script:
#   1) Reads spatial indicator time series (res_table.csv) for one species
#   2) Fits GAM-smoothed trends by season
#   3) Detects change points using:
#        - Bayesian Change Point (BCP)
#        - Nonparametric PELT (changepoint.np)
#   4) Makes figures + a change point summary table
#
# NOTE ON DATA ACCESS:
# Raw CT DEEP trawl data are not publicly available. This script assumes
# the user has obtained permission and has local access to derived inputs
# the res.table you read in is from running and saving out the results from the spatial indicator code
# ============================================================


library(bcp)
library(changepoint.np) 
library(ggplot2)
library(gridExtra)
library(mgcv)  
library(dplyr)

# Load data
res_table <- read.csv("C:/Users/oberc/OneDrive/Desktop/Flounders/Winter/res_table.csv")

# Rename columns 
names(res_table)[names(res_table) == "Year"] <- "year"

# Separate data by season
fall_data <- res_table[res_table$Season == "FA", ]
spring_data <- res_table[res_table$Season == "SP", ]

# Fit GAMs for each indicator for Fall and Spring separately
fall_data$AI_gam <- predict(gam(AI ~ s(year), data = fall_data, method = "REML"))
fall_data$CG_Lat_gam <- predict(gam(CG_Lat ~ s(year), data = fall_data, method = "REML"))
fall_data$CG_Lon_gam <- predict(gam(CG_Lon ~ s(year), data = fall_data, method = "REML"))
fall_data$Evenness_gam <- predict(gam(Evenness ~ s(year), data = fall_data, method = "REML"))
fall_data$Inertia_gam <- predict(gam(Inertia ~ s(year), data = fall_data, method = "REML"))
fall_data$Pos_Area_gam <- predict(gam(Pos_Area ~ s(year), data = fall_data, method = "REML"))

spring_data$AI_gam <- predict(gam(AI ~ s(year), data = spring_data, method = "REML"))
spring_data$CG_Lat_gam <- predict(gam(CG_Lat ~ s(year), data = spring_data, method = "REML"))
spring_data$CG_Lon_gam <- predict(gam(CG_Lon ~ s(year), data = spring_data, method = "REML"))
spring_data$Evenness_gam <- predict(gam(Evenness ~ s(year), data = spring_data, method = "REML"))
spring_data$Inertia_gam <- predict(gam(Inertia ~ s(year), data = spring_data, method = "REML"))
spring_data$Pos_Area_gam <- predict(gam(Pos_Area ~ s(year), data = spring_data, method = "REML"))

#  empty data frame to store change point summary
change_point_summary <- data.frame(
  Indicator = character(),
  Season = character(),
  Method = character(),
  Change_Point_Year = integer(),
  Metric = numeric(),
  stringsAsFactors = FALSE
)

# add BCP and PELT change points to summary
add_change_points <- function(gam_data, years, season_label, indicator_name) {
  # BCP Analysis
  bcp_result <- bcp(gam_data, return.mcmc = TRUE)
  bcp_years <- years[which(bcp_result$posterior.prob > 0.55)]
  bcp_probs <- bcp_result$posterior.prob[which(bcp_result$posterior.prob > 0.55)]
  
  for (i in seq_along(bcp_years)) {
    change_point_summary <<- rbind(change_point_summary, data.frame(
      Indicator = indicator_name,
      Season = season_label,
      Method = "BCP",
      Change_Point_Year = bcp_years[i],
      Metric = bcp_probs[i]
    ))
  }
  
  # Nonparametric PELT Analysis
  cpt_result <- cpt.np(gam_data, method = "PELT", minseglen = 5)
  pelt_years <- years[cpts(cpt_result)]
  
  for (i in seq_along(pelt_years)) {
    change_point_summary <<- rbind(change_point_summary, data.frame(
      Indicator = indicator_name,
      Season = season_label,
      Method = "PELT",
      Change_Point_Year = pelt_years[i],
      Metric = NA
    ))
  }
}

# Process each indicator for Fall season
add_change_points(fall_data$AI_gam, fall_data$year, "Fall", "Abundance Index")
add_change_points(fall_data$CG_Lat_gam, fall_data$year, "Fall", "Center of Gravity Latitude")
add_change_points(fall_data$CG_Lon_gam, fall_data$year, "Fall", "Center of Gravity Longitude")
add_change_points(fall_data$Evenness_gam, fall_data$year, "Fall", "Evenness")
add_change_points(fall_data$Inertia_gam, fall_data$year, "Fall", "Inertia")
add_change_points(fall_data$Pos_Area_gam, fall_data$year, "Fall", "Positive Area")

# Process each indicator for Spring season
add_change_points(spring_data$AI_gam, spring_data$year, "Spring", "Abundance Index")
add_change_points(spring_data$CG_Lat_gam, spring_data$year, "Spring", "Center of Gravity Latitude")
add_change_points(spring_data$CG_Lon_gam, spring_data$year, "Spring", "Center of Gravity Longitude")
add_change_points(spring_data$Evenness_gam, spring_data$year, "Spring", "Evenness")
add_change_points(spring_data$Inertia_gam, spring_data$year, "Spring", "Inertia")
add_change_points(spring_data$Pos_Area_gam, spring_data$year, "Spring", "Positive Area")


# Function to perform BCP and PELT analysis on GAM-smoothed predictions and plot
plot_change_points <- function(raw_data, gam_data, years, season_label, indicator_name, y_label) {
  # BCP Analysis
  bcp_result <- bcp(gam_data, return.mcmc = TRUE)
  bcp_years <- years[which(bcp_result$posterior.prob > 0.55)]
  
  # Nonparametric PELT Analysis
  cpt_result <- cpt.np(gam_data, method = "PELT", minseglen = 5)
  pelt_years <- years[cpts(cpt_result)]
  
  # Plot data with change points
  plot_data <- data.frame(Year = years, Raw = raw_data, GAM = gam_data)
  p <- ggplot(plot_data, aes(x = Year)) +
    geom_line(aes(y = Raw), color = "black", size = 0.8) +
    geom_line(aes(y = GAM), color = "green", size = 0.8) +
    geom_vline(xintercept = bcp_years, color = "red", linetype = "dashed", size = 1.0) +
    geom_vline(xintercept = pelt_years, color = "blue", linetype = "dotted", size = 1.0) +
    labs(title = paste("Change Points -", indicator_name, "(", season_label, ")"),
         x = "Year", y = y_label) +
    theme_minimal()
  return(p)
}


# Plotting
fall_plot_ai <- plot_change_points(fall_data$AI, fall_data$AI_gam, fall_data$year, "Fall", "Abundance Index", "Abundance Index")
fall_plot_cg_lat <- plot_change_points(fall_data$CG_Lat, fall_data$CG_Lat_gam, fall_data$year, "Fall", "Center of Gravity Latitude", "Center of Gravity Latitude")
fall_plot_cg_lon <- plot_change_points(fall_data$CG_Lon, fall_data$CG_Lon_gam, fall_data$year, "Fall", "Center of Gravity Longitude", "Center of Gravity Longitude")
fall_plot_evenness <- plot_change_points(fall_data$Evenness, fall_data$Evenness_gam, fall_data$year, "Fall", "Evenness", "Evenness")
fall_plot_inertia <- plot_change_points(fall_data$Inertia, fall_data$Inertia_gam, fall_data$year, "Fall", "Inertia", "Inertia")
fall_plot_pos_area <- plot_change_points(fall_data$Pos_Area, fall_data$Pos_Area_gam, fall_data$year, "Fall", "Positive Area", "Positive Area")

spring_plot_ai <- plot_change_points(spring_data$AI, spring_data$AI_gam, spring_data$year, "Spring", "Abundance Index", "Abundance Index")
spring_plot_cg_lat <- plot_change_points(spring_data$CG_Lat, spring_data$CG_Lat_gam, spring_data$year, "Spring", "Center of Gravity Latitude", "Center of Gravity Latitude")
spring_plot_cg_lon <- plot_change_points(spring_data$CG_Lon, spring_data$CG_Lon_gam, spring_data$year, "Spring", "Center of Gravity Longitude", "Center of Gravity Longitude")
spring_plot_evenness <- plot_change_points(spring_data$Evenness, spring_data$Evenness_gam, spring_data$year, "Spring", "Evenness", "Evenness")
spring_plot_inertia <- plot_change_points(spring_data$Inertia, spring_data$Inertia_gam, spring_data$year, "Spring", "Inertia", "Inertia")
spring_plot_pos_area <- plot_change_points(spring_data$Pos_Area, spring_data$Pos_Area_gam, spring_data$year, "Spring", "Positive Area", "Positive Area")

# separate legend
legend_data <- data.frame(
  x = c(1, 2),
  y = c(1, 1),
  line_type = factor(c("BCP Change Point", "PELT Change Point"), 
                     levels = c("BCP Change Point", "PELT Change Point"))
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = line_type, linetype = line_type)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c("BCP Change Point" = "red", "PELT Change Point" = "blue")) +
  scale_linetype_manual(values = c("BCP Change Point" = "dashed", "PELT Change Point" = "dotted")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.box.just = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))


legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = line_type, linetype = line_type)) +
  geom_line(size = 3) +  # Make legend lines thicker
  scale_color_manual(values = c("BCP Change Point" = "red", "PELT Change Point" = "blue")) +
  scale_linetype_manual(values = c("BCP Change Point" = "dashed", "PELT Change Point" = "dotted")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),  
    legend.key.width = unit(2, "cm"),       
    legend.key.height = unit(0.7, "cm"),    
    legend.box.just = "center",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Combine legend with fall and spring plots
fall_final <- grid.arrange(
  fall_plot_ai, fall_plot_cg_lat, fall_plot_cg_lon,
  fall_plot_evenness, fall_plot_inertia, fall_plot_pos_area,
  legend_plot,  
  ncol = 3, 
  nrow = 3,
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6),
                        c(7, 7, 7)),  
  heights = c(1, 1, 0.1)  # Adjust
)

spring_final <- grid.arrange(
  spring_plot_ai, spring_plot_cg_lat, spring_plot_cg_lon,
  spring_plot_evenness, spring_plot_inertia, spring_plot_pos_area,
  legend_plot,  
  ncol = 3, 
  nrow = 3,
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6),
                        c(7, 7, 7)),  
  heights = c(1, 1, 0.1)  
)



# Save 
ggsave("Change_Point_Analysis_summer_fall3.png", fall_final, width = 15, height = 12)
ggsave("Change_Point_Analysis_summer_spring3.png", spring_final, width = 15, height = 12)

# Save the change point summary
write.csv(change_point_summary, "Change_Point_Summary_summer2.csv", row.names = FALSE)


# Update the change_point_summary
change_point_summary <- change_point_summary %>%
  group_by(Indicator, Season) %>%
  mutate(
    Overlap_BCP_PELT = Change_Point_Year %in% Change_Point_Year[Method == "BCP"] & 
      Change_Point_Year %in% Change_Point_Year[Method == "PELT"]
  ) %>%
  ungroup()


# Print 
print(change_point_summary)



# Combine all plots into a single figure: spring on top, fall on bottom
all_plots_combined <- grid.arrange(
  spring_plot_ai, spring_plot_cg_lat, spring_plot_cg_lon,
  spring_plot_evenness, spring_plot_inertia, spring_plot_pos_area,
  fall_plot_ai, fall_plot_cg_lat, fall_plot_cg_lon,
  fall_plot_evenness, fall_plot_inertia, fall_plot_pos_area,
  legend_plot,
  ncol = 3,
  nrow = 5,
  layout_matrix = rbind(
    c(1, 2, 3),   
    c(4, 5, 6),   
    c(7, 8, 9),  
    c(10, 11, 12),
    c(13, 13, 13) 
  ),
  heights = c(1, 1, 1, 1, 0.1)
)

# Save 
#gsave("Change_Point_Analysis_summer_flounder_combined.png", all_plots_combined, width = 15, height = 21)





