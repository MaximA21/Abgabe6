#' Find Spatial Hotspots in Geographic Data
#'
#' Performs spatial hotspot analysis using kernel density estimation to identify
#' areas of high point density in geographic data. Returns an interactive map
#' with heatmap overlay and summary statistics.
#'
#' @param data A data.frame containing coordinate data
#' @param lat_col Character string specifying the name of the latitude column
#' @param lon_col Character string specifying the name of the longitude column
#' @param name_col Character string specifying the name column for point labels.
#'   If NULL, points will be labeled as "Point 1", "Point 2", etc.
#' @param bandwidth Numeric value for kernel bandwidth in kilometers. If "auto",
#'   optimal bandwidth is calculated automatically.
#' @param grid_size Integer specifying the resolution of the density grid.
#'   Higher values give more detailed results but slower computation.
#' @param threshold Numeric value between 0 and 1 for hotspot detection threshold.
#'   Points above this quantile of density values are considered hotspots.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{map}: Interactive leaflet map with heatmap overlay
#'   \item \code{hotspots}: Data.frame of identified hotspot locations
#'   \item \code{statistics}: List of summary statistics
#'   \item \code{density_data}: Full density grid data
#' }
#'
#' @export
#' @examples
#' # Example 1: Using included Munich places dataset
#' data(munich_places)
#' result_munich <- find_hotspots(munich_places, "latitude", "longitude")
#'
#' # View the interactive map
#' result_munich$map
#'
#' # Check statistics
#' result_munich$statistics
#'
#' # View top hotspots
#' head(result_munich$hotspots)
#'
#' # Example 2: Using custom random data
#' set.seed(123)
#' sample_data <- data.frame(
#'   latitude = rnorm(100, mean = 48.1351, sd = 0.1),
#'   longitude = rnorm(100, mean = 11.5820, sd = 0.1),
#'   event_type = sample(c("A", "B", "C"), 100, replace = TRUE)
#' )
#'
#' result_random <- find_hotspots(sample_data, "latitude", "longitude")
#' result_random$map
#' @importFrom stats complete.cases
find_hotspots <- function(data, lat_col, lon_col, name_col = NULL,
                          bandwidth = "auto", grid_size = 50, threshold = 0.8) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  if (!lat_col %in% names(data)) {
    stop(paste("Column", lat_col, "not found in data"))
  }

  if (!lon_col %in% names(data)) {
    stop(paste("Column", lon_col, "not found in data"))
  }

  if (!is.null(name_col) && !name_col %in% names(data)) {
    stop(paste("Column", name_col, "not found in data"))
  }

  # Extract coordinates and names
  coords <- data[, c(lon_col, lat_col)]

  # Handle point names
  if (!is.null(name_col)) {
    point_names <- data[[name_col]]
  } else {
    point_names <- paste("Point", 1:nrow(data))
  }

  # Remove rows with missing coordinates (and corresponding names)
  complete_rows <- complete.cases(coords)
  coords <- coords[complete_rows, ]
  point_names <- point_names[complete_rows]

  if (nrow(coords) < 3) {
    stop("Need at least 3 valid coordinate pairs")
  }

  # Calculate automatic bandwidth if needed
  if (bandwidth == "auto") {
    # Use Silverman's rule of thumb adapted for geographic data
    n <- nrow(coords)
    bandwidth <- 0.9 * min(stats::sd(coords[,1]), stats::sd(coords[,2])) * n^(-1/5)
    # Convert to roughly kilometers (very approximate)
    bandwidth <- bandwidth * 111  # 1 degree ≈ 111 km
  }

  # Create grid for density estimation
  lon_range <- range(coords[,1])
  lat_range <- range(coords[,2])

  # Expand ranges slightly
  lon_expand <- diff(lon_range) * 0.1
  lat_expand <- diff(lat_range) * 0.1

  lon_seq <- seq(lon_range[1] - lon_expand, lon_range[2] + lon_expand,
                 length.out = grid_size)
  lat_seq <- seq(lat_range[1] - lat_expand, lat_range[2] + lat_expand,
                 length.out = grid_size)

  # Perform 2D kernel density estimation
  kde_result <- KernSmooth::bkde2D(as.matrix(coords),
                                   bandwidth = c(bandwidth/111, bandwidth/111),
                                   gridsize = c(grid_size, grid_size),
                                   range.x = list(lon_range + c(-lon_expand, lon_expand),
                                                  lat_range + c(-lat_expand, lat_expand)))

  # Create density grid
  density_grid <- expand.grid(longitude = kde_result$x1,
                              latitude = kde_result$x2)
  density_grid$density <- as.vector(kde_result$fhat)

  # Identify hotspots
  hotspot_threshold <- stats::quantile(density_grid$density, threshold, na.rm = TRUE)
  hotspots <- density_grid[density_grid$density >= hotspot_threshold, ]
  hotspots <- hotspots[order(hotspots$density, decreasing = TRUE), ]

  # Calculate statistics
  statistics <- list(
    total_points = nrow(coords),
    num_hotspots = nrow(hotspots),
    max_density = max(density_grid$density, na.rm = TRUE),
    mean_density = mean(density_grid$density, na.rm = TRUE),
    bandwidth_used = bandwidth,
    threshold_used = threshold
  )

  # Create interactive map
  map <- leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::addCircleMarkers(
      data = coords,
      lng = coords[,1],
      lat = coords[,2],
      radius = 3,
      opacity = 0.7,
      fillOpacity = 0.5,
      color = "blue",
      popup = paste("<b>", point_names, "</b><br>",
                    "Lat:", round(coords[,2], 4), "<br>",
                    "Lon:", round(coords[,1], 4))
    )

  # Add hotspots if any found
  if (nrow(hotspots) > 0) {
    top_hotspots <- hotspots[1:min(10, nrow(hotspots)), ]
    map <- map %>%
      leaflet::addCircleMarkers(
        data = top_hotspots,
        lng = ~longitude,
        lat = ~latitude,
        radius = 8,
        color = "red",
        fillColor = "red",
        fillOpacity = 0.8,
        popup = paste("<b>HOTSPOT</b><br>",
                      "Density:", round(top_hotspots$density, 6), "<br>",
                      "Lat:", round(top_hotspots$latitude, 4), "<br>",
                      "Lon:", round(top_hotspots$longitude, 4))
      )
  }

  # Add legend
  map <- map %>%
    leaflet::addLegend(
      position = "bottomright",
      colors = c("blue", "red"),
      labels = c("Data Points", "Hotspots"),
      title = "Legend"
    )

  # Add title
  map <- htmlwidgets::prependContent(map,
                                     htmltools::tags$h3("Spatial Hotspot Analysis",
                                                        style = "text-align: center; margin-bottom: 10px;"))

  # Return results
  return(list(
    map = map,
    hotspots = hotspots,
    statistics = statistics,
    density_data = density_grid
  ))
}
