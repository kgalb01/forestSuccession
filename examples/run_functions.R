# Delete all variables in the workspace
rm(list = ls())

# Set working directory, change this to your specific directory
setwd("C:/Users/kgalb/Documents/Workspace/R/LiDAR")

# set and load library
librarys <- c("lidR", "TreeLS", "RCSF", "sf", "rgl", "raster", "sp", "dplyr", "ggplot2")

# function that loads all libraries
load_librarys <- function(librarys) {
  for (library in librarys) {
    # Check, if lib is installed
    if (!require(library, character.only = TRUE)) {
      # install, if not
      install.packages(library)
      library(library, character.only = TRUE)
    }
  }
}

# call function
load_librarys(librarys)

# read data
loc1_t1_cropped <- lidR::readLAS("O1/t1_cropped.laz")
loc1_t2_cropped <- lidR::readLAS("O1/t2_cropped.laz")
loc1_t3_cropped <- lidR::readLAS("O1/t3_cropped.laz")

# run functions for all three locations
########################################################################################################
# location 1
# calculate DOM for each time set
loc1_t1_dsm <- createDSM(loc1_t1_cropped)
loc1_t2_dsm <- createDSM(loc1_t2_cropped)
loc1_t3_dsm <- createDSM(loc1_t3_cropped)

# calculate DGM for each time set
loc1_t1_dtm <- createDTM(loc1_t1_cropped)
loc1_t2_dtm <- createDTM(loc1_t2_cropped)
loc1_t3_dtm <- createDTM(loc1_t3_cropped)

# normalize the data of each time set
loc1_t1_norm <- normalize(loc1_t1_cropped, loc1_t1_dtm)
loc1_t2_norm <- normalize(loc1_t2_cropped, loc1_t2_dtm)
loc1_t3_norm <- normalize(loc1_t3_cropped, loc1_t3_dtm)

loc1_t1_ttops <- locate_trees(loc1_t1_norm, lmf(ws = 5))
loc1_t2_ttops <- locate_trees(loc1_t2_norm, lmf(ws = 5))
loc1_t3_ttops <- locate_trees(loc1_t3_norm, lmf(ws = 5))

# plot the normalized data
x <- plot(loc1_t1_norm, size = 4, color = "Z", bg = "white", axis = TRUE, legend = FALSE)
add_treetops3d(x, loc1_t1_ttops)
rgl::title3d("Location 1 - 2011", line = -0.1)

y <- plot(loc1_t2_norm, size = 4, color = "Z", bg = "white", axis = TRUE, legend = FALSE)
add_treetops3d(y, loc1_t2_ttops)
rgl::title3d("Location 1 - 2015", line = -0.1)

z <- plot(loc1_t3_norm, size = 4, color = "Z", bg = "white", axis = TRUE, legend = FALSE)
add_treetops3d(z, loc1_t3_ttops)
rgl::title3d("Location 1 - 2021", line = -0.1)

# calculate CHM for each time set
loc1_t1_chm <- createCHM(loc1_t1_norm)
loc1_t2_chm <- createCHM(loc1_t2_norm)
loc1_t3_chm <- createCHM(loc1_t3_norm)

# remove trees for each time set
loc1_t1_no_tree <- remove_trees(loc1_t1_norm)
loc1_t2_no_tree <- remove_trees(loc1_t2_norm)
loc1_t3_no_tree <- remove_trees(loc1_t3_norm)

# plot the data without trees
plot(loc1_t1_no_tree, size = 4, color = "Z", bg = "white", axis = TRUE, legend = FALSE)
rgl::title3d("Location 1 - 2011, No Trees", line = -0.1)
plot(loc1_t2_no_tree, size = 4, color = "Z", bg = "white", axis = TRUE, legend = FALSE)
rgl::title3d("Location 1 - 2015, No Trees", line = -0.1)
plot(loc1_t3_no_tree, size = 4, color = "Z", bg = "white", axis = TRUE, legend = FALSE)
rgl::title3d("Location 1 - 2021, No Trees", line = -0.1)

# Plot the low vegetation differently ...
loc1_t1_no_tree_no_ground <- remove_ground(loc1_t1_no_tree)
loc1_t1_no_tree_no_ground_sf <- convertToSF(loc1_t1_no_tree_no_ground)
coords <- st_coordinates(loc1_t1_no_tree_no_ground_sf)

plot1 <- ggplot(loc1_t1_no_tree_no_ground_sf) +
  geom_sf(aes(color = Z), size = 2) +
  scale_color_viridis_c(name = "Height") +
  theme_minimal() +
  labs(title = "Spatial Distribution of low Vegetation of Location 1 between in 2011",
       subtitle = "Color represents the height (Z-value)",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
        plot.subtitle = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 12, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  coord_sf(xlim = c(min(coords[, 1]) - 10, max(coords[, 1]) + 10))

# 2015 ...
loc1_t2_no_tree_no_ground <- remove_ground(loc1_t2_no_tree)
loc1_t2_no_tree_no_ground_sf <- convertToSF(loc1_t2_no_tree_no_ground)
coords <- st_coordinates(loc1_t2_no_tree_no_ground_sf)

ggplot(loc1_t2_no_tree_no_ground_sf) +
  geom_sf(aes(color = Z), size = 2) +
  scale_color_viridis_c(name = "Height") +
  theme_minimal() +
  labs(title = "Spatial Distribution of low Vegetation of Location 1 in 2015",
       subtitle = "Color represents the height (Z-value)",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
        plot.subtitle = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 12, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  coord_sf(xlim = c(min(coords[, 1]) - 10, max(coords[, 1]) + 10))

# 2021 ..
loc1_t3_no_tree_no_ground <- remove_ground(loc1_t3_no_tree)
loc1_t3_no_tree_no_ground_sf <- convertToSF(loc1_t3_no_tree_no_ground)
coords <- st_coordinates(new_points_loc1_t2 )

ggplot(new_points_loc1_t2 ) +
  geom_sf(aes(color = Z), size = 2) +
  scale_color_viridis_c(name = "Height") +
  theme_minimal() +
  labs(title = "Spatial Distribution of low Vegetation of Location 1 in 2021",
       subtitle = "Color represents the height (Z-value)",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
        plot.subtitle = element_text(size = 12, family = "Arial"),
        axis.title = element_text(size = 12, family = "Arial"),
        axis.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial")) +
  coord_sf(xlim = c(min(coords[, 1]) - 10, max(coords[, 1]) + 10))

# Convert LAS Objects to RasterLayer Objects
loc1_t1_norm_sf <- convertToSF(loc1_t1_norm)
loc1_t2_norm_sf <- convertToSF(loc1_t2_norm)
loc1_t3_norm_sf <- convertToSF(loc1_t3_norm)

loc1_t1_norm_rast <- createRaster(loc1_t1_norm_sf)
loc1_t2_norm_rast <- createRaster(loc1_t2_norm_sf)
loc1_t3_norm_rast <- createRaster(loc1_t3_norm_sf)

loc1_t1_norm_rast <- fillRaster(loc1_t1_norm_sf, loc1_t1_norm_rast)
loc1_t2_norm_rast <- fillRaster(loc1_t2_norm_sf, loc1_t2_norm_rast)
loc1_t3_norm_rast <- fillRaster(loc1_t3_norm_sf, loc1_t3_norm_rast)

# Convert low Vegetation LAS Object to RasterLayer Object
loc1_t1_no_tree_sf <- convertToSF(loc1_t1_no_tree)
loc1_t2_no_tree_sf <- convertToSF(loc1_t2_no_tree)
loc1_t3_no_tree_sf <- convertToSF(loc1_t3_no_tree)

loc1_t1_no_tree_rast <- createRaster(loc1_t1_no_tree_sf)
loc1_t2_no_tree_rast <- createRaster(loc1_t2_no_tree_sf)
loc1_t3_no_tree_rast <- createRaster(loc1_t3_no_tree_sf)

loc1_t1_no_tree_rast <- fillRaster(loc1_t1_no_tree_sf, loc1_t1_no_tree_rast)
loc1_t2_no_tree_rast <- fillRaster(loc1_t2_no_tree_sf, loc1_t2_no_tree_rast)
loc1_t3_no_tree_rast <- fillRaster(loc1_t3_no_tree_sf, loc1_t3_no_tree_rast)

# append further info about the trees on the data set
loc1_t1_norm_tree_features <- tree_info(loc1_t1_norm)
loc1_t2_norm_tree_features <- tree_info(loc1_t2_norm)
loc1_t3_norm_tree_features <- tree_info(loc1_t3_norm)

# calculate the biomass and intensity for each time set
ggplot(veg_cluster(loc1_t1_norm), aes(x = Z, y = Intensity, color = factor(ClusterLabel))) +
  geom_point(size = 3) +
  labs(x = "Höhe (Meter)", y = "Intensität",
       title = "Vegetations-Clustering under 5 Meters at Location 1, 2011") +
  scale_color_discrete(name = "Cluster")

ggplot(veg_cluster(loc1_t2_norm), aes(x = Z, y = Intensity, color = factor(ClusterLabel))) +
  geom_point(size = 3) +
  labs(x = "Höhe (Meter)", y = "Intensität",
       title = "Vegetations-Clustering under 5 Meters at Location 1, 2015") +
  scale_color_discrete(name = "Cluster")

ggplot(veg_cluster(loc1_t3_norm), aes(x = Z, y = Intensity, color = factor(ClusterLabel))) +
  geom_point(size = 3) +
  labs(x = "Höhe (Meter)", y = "Intensität",
       title = "Vegetations-Clustering under 5 Meters at Location 1, 2021") +
  scale_color_discrete(name = "Cluster")

# rasterize the vegetation for each time set
loc1_t1_df <- convertToSF(loc1_t1_norm)
loc1_t1_raster <- createRaster(loc1_t1_df)
loc1_t1_veg_raster <- fillRaster(loc1_t1_df, loc1_t1_raster)

loc1_t2_df <- convertToSF(loc1_t2_norm)
loc1_t2_raster <- createRaster(loc1_t2_df)
loc1_t2_veg_raster <- fillRaster(loc1_t2_df, loc1_t2_raster)

loc1_t3_df <- convertToSF(loc1_t3_norm)
loc1_t3_raster <- createRaster(loc1_t3_df)
loc1_t3_veg_raster <- fillRaster(loc1_t3_df, loc1_t3_raster)

# rasterize the low vegetation for each time set
loc1_t1_no_trees_df <- convertToSF(loc1_t1_no_tree)
loc1_t1_no_trees_raster <- createRaster(loc1_t1_no_trees_df)
loc1_t1_no_trees_veg_raster <- fillRaster(loc1_t1_no_trees_df, loc1_t1_no_trees_raster)

loc1_t2_no_trees_df <- convertToSF(loc1_t2_no_tree)
loc1_t2_no_trees_raster <- createRaster(loc1_t2_no_trees_df)
loc1_t2_no_trees_veg_raster <- fillRaster(loc1_t2_no_trees_df, loc1_t2_no_trees_raster)

loc1_t3_no_trees_df <- convertToSF(loc1_t3_no_tree)
loc1_t3_no_trees_raster <- createRaster(loc1_t3_no_trees_df)
loc1_t3_no_trees_veg_raster <- fillRaster(loc1_t3_no_trees_df, loc1_t3_no_trees_raster)

# calculate vegetation volume for each time set
loc1_t1_veg_volume <- volume(loc1_t1_veg_raster)
loc1_t1_veg_volume_total <- loc1_t1_veg_volume$total_volume

loc1_t2_veg_volume <- volume(loc1_t2_veg_raster)
loc1_t2_veg_volume_total <- loc1_t2_veg_volume$total_volume

loc1_t3_veg_volume <- volume(loc1_t3_veg_raster)
loc1_t3_veg_volume_total <- loc1_t3_veg_volume$total_volume

# print and plot volume values
print(loc1_t1_veg_volume)
print(loc1_t2_veg_volume)
print(loc1_t3_veg_volume)
par(mfrow=c(1,3))
plot(loc1_t1_veg_volume$volume_raster, main = "Vegetation volume at Location 1, 2011", col = terrain.colors(100))
plot(loc1_t2_veg_volume$volume_raster, main = "Vegetation volume at Location 1, 2015", col = terrain.colors(100))
plot(loc1_t3_veg_volume$volume_raster, main = "Vegetation volume at Location 1, 2021", col = terrain.colors(100))

# calculate vegetation volume for each time set without trees
loc1_t1_no_trees_veg_volume <- volume(loc1_t1_no_trees_veg_raster)
loc1_t1_no_trees_veg_volume_total <- loc1_t1_no_trees_veg_volume$total_volume

loc1_t2_no_trees_veg_volume <- volume(loc1_t2_no_trees_veg_raster)
loc1_t2_no_trees_veg_volume_total <- loc1_t2_no_trees_veg_volume$total_volume

loc1_t3_no_trees_veg_volume <- volume(loc1_t3_no_trees_veg_raster)
loc1_t3_no_trees_veg_volume_total <- loc1_t3_no_trees_veg_volume$total_volume

# print and plot volume values
print(loc1_t1_no_trees_veg_volume)
print(loc1_t2_no_trees_veg_volume)
print(loc1_t3_no_trees_veg_volume)
par(mfrow=c(1,3))
plot(loc1_t1_no_trees_veg_volume$volume_raster, main = "Vegetation volume of low vegetation at Location 1, 2011", col = terrain.colors(100))
plot(loc1_t2_no_trees_veg_volume$volume_raster, main = "Vegetation volume of low vegetation at Location 1, 2015", col = terrain.colors(100))
plot(loc1_t3_no_trees_veg_volume$volume_raster, main = "Vegetation volume of low vegetation at Location 1, 2021", col = terrain.colors(100))

# calculate average vegetation values of low vegetation
loc1_t1_avg_veg <- cellStats(loc1_t1_veg_raster, stat = mean, na.rm = TRUE)
loc1_t2_avg_veg <- cellStats(loc1_t2_veg_raster, stat = mean, na.rm = TRUE)
loc1_t3_avg_veg <- cellStats(loc1_t3_veg_raster, stat = mean, na.rm = TRUE)

data <- data.frame(
  Timepoint = c("2011", "2015", "2021"),
  AverageVegetation = c(loc1_t1_avg_veg, loc1_t2_avg_veg, loc1_t3_avg_veg)
)

# Creating the bar chart
ggplot(data, aes(x = Timepoint, y = AverageVegetation, fill = Timepoint)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Timepoint",
    y = "Average Vegetation Height",
    title = "Average Vegetation Height Over Timepoints in Location 1"
  ) +
  theme_minimal()

# calculate and print number of points for low vegetation
loc1_t1_low_veg <- sum(loc1_t1_veg_raster[] <= 3, na.rm = TRUE)
loc1_t2_low_veg <- sum(loc1_t2_veg_raster[] <= 3, na.rm = TRUE)
loc1_t3_low_veg <- sum(loc1_t3_veg_raster[] <= 3, na.rm = TRUE)

data <- data.frame(
  Timepoint = c("2011", "2015", "2021"),
  low_veg_points = c(loc1_t1_low_veg, loc1_t2_low_veg, loc1_t3_low_veg)
)

# Creating the bar chart
ggplot(data, aes(x = Timepoint, y = low_veg_points, fill = Timepoint)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Timepoint",
    y = "Number of low Vegetation Points",
    title = "Number of Points for Low Vegetation in Location 1"
  ) +
  theme_minimal()

# calculate sd for low vegetation as a degree for diversity
loc1_t1_low_veg_sd <- sd(loc1_t1_veg_raster[], na.rm = TRUE)
loc1_t2_low_veg_sd <- sd(loc1_t2_veg_raster[], na.rm = TRUE)
loc1_t3_low_veg_sd <- sd(loc1_t3_veg_raster[], na.rm = TRUE)

data <- data.frame(
  Timepoint = c("2011", "2015", "2021"),
  low_veg_points = c(loc1_t1_low_veg_sd, loc1_t2_low_veg_sd, loc1_t3_low_veg_sd)
)

# Creating the bar chart
ggplot(data, aes(x = Timepoint, y = low_veg_points, fill = Timepoint)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Timepoint",
    y = "Standard Deviation of low Vegetation",
    title = "Standard deviation of low Vegetation in Location 1"
  ) +
  theme_minimal()

# Calculate and plot differences:
# Height difference as RasterLayer Object
height_diff_rast_loc1_t2_t1 <- loc1_t2_norm_rast - loc1_t1_norm_rast
height_diff_rast_loc1_t3_t1 <- loc1_t3_norm_rast - loc1_t1_norm_rast
height_diff_rast_loc1_t3_t2 <- loc1_t3_norm_rast - loc1_t2_norm_rast

# remove negativ values, only show where vegetation grew
rcl_pos <- matrix(c(-Inf, 6, NA), ncol=3, byrow=TRUE)
new_veg_rast_loc1_t2_t1 <- reclassify(height_diff_rast_loc1_t2_t1, rcl_pos)
new_veg_rast_loc1_t3_t1 <- reclassify(height_diff_rast_loc1_t3_t1, rcl_pos)
new_veg_rast_loc1_t3_t2 <- reclassify(height_diff_rast_loc1_t3_t2, rcl_pos)

par(mfrow=c(1,2))
plot(new_veg_rast_loc1_t2_t1, main="New Vegetation at Location 1 between 2011 and 2015", col = terrain.colors(50), legend = TRUE)
plot(new_veg_rast_loc1_t3_t2, main="New Vegetation at Location 1 between 2015 and 2021", col = terrain.colors(50), legend = TRUE)

# remove positive values, only show where vegetation is lost
rcl_neg <- matrix(c(-5, Inf, NA), ncol=3, byrow=TRUE)
lost_veg_rast_loc1_t2_t1 <- reclassify(height_diff_rast_loc1_t2_t1, rcl_neg)
lost_veg_rast_loc1_t3_t2 <- reclassify(height_diff_rast_loc1_t3_t2, rcl_neg)

par(mfrow=c(1,2))
plot(lost_veg_rast_loc1_t2_t1, main="Lost Vegetation at Location 1 between 2011 and 2015", col = rev(height.colors(50)), legend = TRUE)
plot(lost_veg_rast_loc1_t3_t2, main="Lost Vegetation at Location 1 between 2015 and 2021", col = rev(height.colors(50)), legend = TRUE)

# Height difference in low vegetation as RasterLayer Object
height_diff_low_rast_loc1_t2_t1 <- loc1_t2_no_tree_rast  - loc1_t1_no_tree_rast
height_diff_low_rast_loc1_t3_t1 <- loc1_t3_no_tree_rast  - loc1_t1_no_tree_rast
height_diff_low_rast_loc1_t3_t2 <- loc1_t3_no_tree_rast  - loc1_t2_no_tree_rast

# remove negativ values, only show where vegetation grew
rcl_pos <- matrix(c(-Inf, 0.1, NA), ncol=3, byrow=TRUE)
new_low_veg_rast_loc1_t2_t1 <- reclassify(height_diff_low_rast_loc1_t2_t1, rcl_pos)
new_low_veg_rast_loc1_t3_t1 <- reclassify(height_diff_low_rast_loc1_t3_t1, rcl_pos)
new_low_veg_rast_loc1_t3_t2 <- reclassify(height_diff_low_rast_loc1_t3_t2, rcl_pos)

par(mfrow=c(1,2))
plot(new_low_veg_rast_loc1_t2_t1, main="New low Vegetation at Location 1 between 2011 and 2015", col = terrain.colors(10), legend = TRUE)
plot(new_low_veg_rast_loc1_t3_t2, main="New low Vegetation at Location 1 between 2015 and 2021", col = terrain.colors(10), legend = TRUE)

# remove positive values, only show where vegetation is lost
rcl_neg <- matrix(c(-0.1, Inf, NA), ncol=3, byrow=TRUE)
lost_low_veg_rast_loc1_t2_t1 <- reclassify(height_diff_low_rast_loc1_t2_t1, rcl_neg)
lost_low_veg_rast_loc1_t3_t2 <- reclassify(height_diff_low_rast_loc1_t3_t2, rcl_neg)

par(mfrow=c(1,2))
plot(lost_low_veg_rast_loc1_t2_t1, main="Lost low Vegetation at Location 1 between 2011 and 2015", col = topo.colors(10), legend = TRUE)
plot(lost_low_veg_rast_loc1_t3_t2, main="Lost low Vegetation at Location 1 between 2015 and 2021", col = topo.colors(10), legend = TRUE)

# average height between the three time frames
timepoints <- c("2011", "2015", "2021")
min_length <- min(length(loc1_t2_norm$Z), length(loc1_t1_norm$Z))
average_height_t1 <- mean(loc1_t1_norm$Z, na.rm = TRUE)
average_height_t2 <- mean(loc1_t2_norm$Z, na.rm = TRUE)
average_height_t3 <- mean(loc1_t3_norm$Z, na.rm = TRUE)
average_heights <- c(average_height_t1, average_height_t2, average_height_t3)
height_data <- data.frame(Timepoint = timepoints, AverageHeight = average_heights)
barplot(height_data$AverageHeight, names.arg = height_data$Timepoint, col = "skyblue",
        main = "Average height at different timepoints",
        xlab = "Timepoint", ylab = "Average Height")

# height differences between the three time frames
min_length <- min(length(loc1_t2_norm$Z), length(loc1_t1_norm$Z))
height_difference_loc1_t2_t1 <- loc1_t2_norm$Z[1:min_length] - loc1_t1_norm$Z[1:min_length]
colors <- ifelse(height_difference_loc1_t2_t1 >= 0, "blue", "red")
plot(loc1_t2_norm$X, loc1_t2_norm$Y, col = colors, pch = 16,
     main = "Scatterplot of Height Differences of Location 1, 2011 & 2015", xlab = "X", ylab = "Y")

min_length <- min(length(loc1_t3_norm$Z), length(loc1_t1_norm$Z))
height_difference_loc1_t3_t1 <- loc1_t3_norm$Z[1:min_length] - loc1_t1_norm$Z[1:min_length]
colors <- ifelse(height_difference_loc1_t3_t1 >= 0, "blue", "red")
plot(loc1_t3_norm$X, loc1_t3_norm$Y, col = colors, pch = 16,
     main = "Scatterplot of Height Differences of Location 1, 2021 & 2011", xlab = "X", ylab = "Y")

min_length <- min(length(loc1_t3_norm$Z), length(loc1_t2_norm$Z))
height_difference_loc1_t3_t2 <- loc1_t3_norm$Z[1:min_length] - loc1_t2_norm$Z[1:min_length]
colors <- ifelse(height_difference_loc1_t3_t2 >= 0, "blue", "red")
plot(loc1_t3_norm$X, loc1_t3_norm$Y, col = colors, pch = 16,
     main = "Scatterplot of Height Differences of Location 1, 2021 & 2015", xlab = "X", ylab = "Y")

# Volume differences
vol_diff_loc1_t3_t1 <- loc1_t3_veg_volume$volume_raster - loc1_t1_veg_volume$volume_raster
vol_diff_loc1_t3_t2 <- loc1_t3_veg_volume$volume_raster - loc1_t2_veg_volume$volume_raster
vol_diff_loc1_t2_t1 <- loc1_t2_veg_volume$volume_raster - loc1_t1_veg_volume$volume_raster
par(mfrow = c(1, 3))
plot(vol_diff_loc1_t2_t1, main = "Volume Difference of Location 1, 2011 & 2015", xlab = "X", ylab = "Y")
plot(vol_diff_loc1_t3_t2, main = "Volume Difference of Location 1, 2015 & 2021", xlab = "X", ylab = "Y")
plot(vol_diff_loc1_t3_t1, main = "Volume Difference of Location 1, 2011 & 2021", xlab = "X", ylab = "Y")

# DGM differences
dtm_diff_loc1_t2_t1 <- loc1_t2_dtm - loc1_t1_dtm
dtm_diff_loc1_t3_t2 <- loc1_t3_dtm - loc1_t2_dtm
par(mfrow = c(1, 2))
plot(dtm_diff_loc1_t2_t1, main = "DGM Difference of Location 1, 2011 & 2015", col = topo.colors(255), legend = TRUE)
plot(dtm_diff_loc1_t3_t2, main = "DGM Difference of Location 1, 2015 & 2021", col = topo.colors(255), legend = TRUE)

# DOM differences
dsm_diff_loc1_t2_t1 <- loc1_t2_dsm - loc1_t1_dsm
dsm_diff_loc1_t3_t2 <- loc1_t3_dsm - loc1_t2_dsm
par(mfrow = c(1, 2))
plot(dsm_diff_loc1_t2_t1, main = "DOM Difference of Location 1, 2011 & 2015", col = topo.colors(255), legend = TRUE)
plot(dsm_diff_loc1_t3_t2, main = "DOM Difference of Location 1, 2015 & 2021", col = topo.colors(255), legend = TRUE)

# CHM differences
chm_diff_loc1_t2_t1 <- loc1_t2_chm - loc1_t1_chm
chm_diff_loc1_t3_t2 <- loc1_t3_chm - loc1_t2_chm
par(mfrow = c(1, 2))
plot(chm_diff_loc1_t2_t1, main = "CHM Difference of Location 1, 2011 & 2015", col = height.colors(50), axes=T, legend = TRUE)
plot(chm_diff_loc1_t3_t2, main = "CHM Difference of Location 1, 2015 & 2021", col = height.colors(50), axes=T, legend = TRUE)

# Vegetation differences
veg_diff_loc1_t2_t1 <- loc1_t2_veg_raster - loc1_t1_veg_raster
veg_diff_loc1_t3_t2 <- loc1_t3_veg_raster - loc1_t2_veg_raster
par(mfrow = c(1, 2))
plot(veg_diff_loc1_t2_t1, main = "Vegetation Difference of Location 1, 2011 & 2015", col = topo.colors(255), legend = TRUE)
plot(veg_diff_loc1_t3_t2, main = "Vegetation Difference of Location 1, 2015 & 2021", col = topo.colors(255), legend = TRUE)
########################################################################################################


########################################################################################################
# location 2
# repeat all this for the second dataset
