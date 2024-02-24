#' createDSM:
#' This function calculates a digital surface model out of lidar data using the grid_canopy function
#' from the "lidR" - Package and the pitfree() algorithm
#' @param data LiDAR point cloud data
#'
#' @return A RasterLayer Object containing the digital surface model
#' @export
#'
#' @examples
#' dsm <- createDSM(myLiDARData)
createDSM <- function(data){
  raster_resolution = 1
  dsm <- grid_canopy(data, res = raster_resolution, algorithm = pitfree())
  return(dsm)
}

#' createDTM:
#' This function calculates a digital terrain model out of lidar data using the grid_terrain function
#' from the "lidR" package and the tin() algorithm
#' @param data LiDAR point cloud data
#'
#' @return A RasterLayer object containing the digital terrain model
#' @export
#'
#' @examples
#' dtm <- createDTM(myLiDARData)
createDTM <- function(data) {
  raster_resolution = 1
  dtm <- grid_terrain(data, res = raster_resolution, algorithm = tin())
  return(dtm)
}

#' createCHM:
#' This function creates a canopy height model out of lidar data using the rasterize_canopy function
#' from the "lidR" package and the p2r() algorithm
#' @param data LiDAR point cloud data
#'
#' @return A RasterLayer object containing the canopy height model
#' @export
#'
#' @examples
#' chm <- createCHM(myLiDARData)
createCHM <- function(data) {
  chm = rasterize_canopy(data, res = 0.5, algorithm = p2r(), pkg = 'stars')
  return(chm)
}

#' normalize:
#' This function normalizes the height of lidar data using the normalize_height function
#' from the "lidR" package and a digital terrain model and classifies ground and noise points
#' @param data LiDAR point cloud data
#' @param dtm Digital terrain model
#'
#' @return A normalized LiDAR point cloud data object
#' @export
#'
#' @examples
#' normalized_data <- normalize(myLiDARData, dtm)
normalize <- function(data, dtm){
  raster_resolution <- 1
  norm <- normalize_height(data, dtm)
  crs(norm) <- sp::CRS("EPSG:25832")
  tree_height_threshold <- 5
  norm@data$Classification[norm@data$Z > tree_height_threshold] <- 5L
  norm@data$Classification[norm@data$Z <= tree_height_threshold] <- 3L
  norm <- classify_ground(norm, csf())
  noise_class = classify_noise(norm, algorithm = ivf(res = 3, n = 50))
  norm = filter_poi(noise_class, Classification != 7)
  return(norm)
}

#' convertToSF:
#' This function converts LiDAR point cloud data to a simple features (SF) data frame
#' using the lidR::as.spatial and st_as_sf functions
#' @param data LiDAR point cloud data
#'
#' @return A simple features (SF) data frame
#' @export
#'
#' @examples
#' sf_data <- convertToSF(myLiDARData)
convertToSF <- function(data) {
  spdf <- lidR::as.spatial(data)
  sf_obj <- st_as_sf(spdf)
  return(sf_obj)
}

#' createRaster:
#' This function creates a rasterized version of LiDAR data using the raster function
#' from the "raster" package
#' @param data LiDAR point cloud data
#' @param resolution Resolution of the raster
#'
#' @return A RasterLayer object
#' @export
#'
#' @examples
#' raster_data <- createRaster(myLiDARData, resolution = 1)
createRaster <- function(data, resolution = 1) {
  ext <- extent(data)
  r <- raster(ext, resolution = resolution)
  projection(r) <- crs(data)
  return(r)
}

#' fillRaster:
#' This function fills a rasterized version of LiDAR data with actual data using the rasterize function
#' from the "raster" package
#' @param data LiDAR point cloud data
#' @param raster RasterLayer object
#'
#' @return A RasterLayer object with filled data
#' @export
#'
#' @examples
#' filled_raster <- fillRaster(myLiDARData, myRaster)
fillRaster <- function(data, raster) {
  rast <- rasterize(data, raster, field = "Z", fun = mean, na.rm = TRUE)
  return(rast)
}

#' volume:
#' This function calculates the volume of the LiDAR point cloud data based on a height threshold
#' @param data LiDAR point cloud data
#' @param height_threshold Height threshold for calculating volume
#'
#' @return A list containing a RasterLayer object and the total volume
#' @export
#'
#' @examples
#' volume_data <- volume(myLiDARData, height_threshold = 5)
volume <- function(data, height_threshold = 5) {
  data[] <- data[] * 1
  total_volume <- sum(data[], na.rm = TRUE)
  return(list("volume_raster" = data, "total_volume" = total_volume))
}

#' identify_trees:
#' This function identifies trees in the LiDAR point cloud data using canopy height model and tree segmentation
#' @param data LiDAR point cloud data
#'
#' @return A list containing tree points and segmentation information
#' @export
#'
#' @examples
#' tree_info <- identify_trees(myLiDARData)
identify_trees <- function(data){
  chm <- createCHM(data)
  tree_points <- data[data@data$Classification == 5L, ]
  treetops <- locate_trees(chm, algorithm = lmf(ws = 3))
  seg <- segment_trees(tree_points, algorithm = dalponte2016(chm, treetops), attribute = "treeID", uniqueness = "incremental")
  return(list(tree_points = tree_points, seg = seg))
}

#' remove_trees:
#' This function removes identified trees from the LiDAR point cloud data and removes them
#' and performs noise removal for optimal results
#' @param data LiDAR point cloud data
#'
#' @return LiDAR point cloud data without trees
#' @export
#'
#' @examples
#' no_trees_data <- remove_trees(myLiDARData)
remove_trees <- function(data) {
  result <- identify_trees(data)
  tree_points <- result$tree_points
  seg <- result$seg
  seg@data$treeID <- seg$treeID
  data@data$treeID <- rep(NA, nrow(data))
  data@data$treeID[data@data$Classification == 5L] <- seg@data$treeID
  no_trees <- data[is.na(data$treeID), ]
  err <- runif(20, -50, 50)
  id <- round(runif(20, 0, npoints(no_trees)))
  no_trees$Z[id] <- no_trees$Z[id] + err
  no_trees <- classify_noise(no_trees, ivf(5, 2))
  no_trees <- filter_poi(no_trees, Classification != LASNOISE)
  no_trees <- filter_poi(no_trees, no_trees$Z <= 3, no_trees$Z >= 0)
  return(no_trees)
}

#' remove_ground:
#' This function removes ground points from the LiDAR point cloud data
#' @param data LiDAR point cloud data
#'
#' @return LiDAR point cloud data without ground points
#' @export
#'
#' @examples
#' no_ground_data <- remove_ground(myLiDARData)
remove_ground <- function(data){
  data <- filter_poi(data, Classification != 2L)
  return(data)
}

#' tree_info:
#' This function provides information about trees in the LiDAR point cloud data, including tree diameters,
#' tree shapes, tree intensities, and tree heights
#' @param data LiDAR point cloud data
#'
#' @return A data frame containing information about trees
#' @export
#'
#' @examples
#' tree_information <- tree_info(myLiDARData)
tree_info <- function(data){
  result <- identify_trees(data)
  tree_points <- result$tree_points
  seg <- result$seg
  seg@data$treeID <- seg$treeID
  chm <- createCHM(data)
  treetops <- locate_trees(chm, algorithm = lmf(ws = 3))
  tree_points@data$treeID <- seg$treeID
  seg_df <- as.data.frame(seg@data)
  tree_points_df <- as.data.frame(tree_points@data)
  tree_points@data$treeID <- rep(seg$treeID, length.out = nrow(tree_points@data))
  tree_diameters <- seg_df %>%
    group_by(treeID) %>%
    summarise(
      MinX = min(X),
      MaxX = max(X),
      MinY = min(Y),
      MaxY = max(Y)
    ) %>%
    mutate(Diameter = sqrt((MaxX - MinX)^2 + (MaxY - MinY)^2))
  tree_shapes <- seg_df %>%
    group_by(treeID) %>%
    summarise(VarianceZ = var(Z))
  tree_intensities <- tree_points_df %>%
    group_by(treeID) %>%
    summarise(
      MeanIntensity = mean(Intensity, na.rm = TRUE),
      MaxIntensity = max(Intensity, na.rm = TRUE)
    )
  tree_height <- seg_df %>%
    group_by(treeID) %>%
    summarise(MaxHeight = max(Z, na.rm = TRUE))
  tree_features <- left_join(tree_diameters, tree_shapes, by = "treeID") %>%
    left_join(tree_intensities, by = "treeID") %>%
    left_join(tree_height, by = "treeID")
  return(tree_features)
}

#' veg_cluster:
#' This function performs k-means clustering on low vegetation points in the LiDAR point cloud data
#' based on height and intensity attributes
#' @param data LiDAR point cloud data
#'
#' @return A data frame containing low vegetation points with cluster labels
#' @export
#'
#' @examples
#' clustered_vegetation <- veg_cluster(myLiDARData)
veg_cluster <- function(data){
  raster_resolution <- 1
  tree_height_threshold <- 5
  low_veg <- data[data@data$Z <= tree_height_threshold,]
  dom_low_veg <- grid_canopy(low_veg, res = raster_resolution, algorithm = pitfree())
  low_veg <- as.data.frame(low_veg@data)
  low_veg <- low_veg %>%
    dplyr::filter(Z > 0 & Z <= 5)
  low_veg$Volume <- sum(low_veg$Z, na.rm = TRUE) * raster_resolution^2
  low_veg_data <- low_veg %>%
    dplyr::filter(Z <= 5) %>%
    dplyr::select(Z, Intensity)
  low_veg_cluster <- kmeans(low_veg_data, centers = 3)
  low_veg$ClusterLabel <- low_veg_cluster$cluster
  return(low_veg)
}
