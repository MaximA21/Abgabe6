#' Create Intelligent Choropleth Maps
#'
#' Creates interactive choropleth maps with automatic data classification and
#' intelligent color palette selection. Supports multiple classification methods
#' and customizable styling options.
#'
#' @param shapefile An sf object containing geographic boundaries (polygons)
#' @param data_col Character string specifying the name of the numeric column
#'   to visualize. This column should contain the values that will determine
#'   the colors of the regions.
#' @param name_col Character string specifying the name column for region labels
#'   in popups. If NULL, uses row numbers.
#' @param classification Character string specifying the classification method.
#'   Options: "jenks" (natural breaks), "quantile", "equal", "sd" (standard deviation)
#' @param n_classes Integer specifying the number of color classes (2-9 recommended)
#' @param palette Character string specifying the color palette.
#'   Options: "viridis", "plasma", "inferno", "magma", "Blues", "Reds", "Greens"
#' @param reverse_palette Logical. Should the color palette be reversed?
#' @param popup_fields Character vector of additional fields to show in popups.
#'   If NULL, only shows name and data value.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{map}: Interactive leaflet choropleth map
#'   \item \code{classification}: Data.frame with classification breaks and statistics
#'   \item \code{data_summary}: Summary statistics of the input data
#'   \item \code{settings}: List of settings used for the map
#' }
#'
#' @importFrom stats quantile median sd
#' @importFrom htmltools tags
#' @export
#' @examples
#' # Example 1: Using included German states dataset
#' data(german_states)
#' result1 <- create_choropleth(german_states, "population_million", "name")
#' result1$map
#'
#' # Example 2: GDP per capita with custom settings
#' result2 <- create_choropleth(
#'   german_states,
#'   "gdp_per_capita",
#'   "name",
#'   classification = "quantile",
#'   n_classes = 4,
#'   palette = "Blues"
#' )
#' result2$map
create_choropleth <- function(shapefile, data_col, name_col = NULL,
                              classification = "jenks", n_classes = 5,
                              palette = "viridis", reverse_palette = FALSE,
                              popup_fields = NULL) {

  # Input validation
  if (!inherits(shapefile, "sf")) {
    stop("shapefile must be an sf object. Use sf::st_read() to load your shapefile.")
  }

  if (!data_col %in% names(shapefile)) {
    stop(paste("Column", data_col, "not found in shapefile"))
  }

  if (!is.null(name_col) && !name_col %in% names(shapefile)) {
    stop(paste("Column", name_col, "not found in shapefile"))
  }

  if (!is.numeric(shapefile[[data_col]])) {
    stop(paste("Column", data_col, "must be numeric"))
  }

  if (!(classification %in% c("jenks", "quantile", "equal", "sd"))) {
    stop("classification must be one of: 'jenks', 'quantile', 'equal', 'sd'")
  }

  if (n_classes < 2 || n_classes > 9) {
    stop("n_classes must be between 2 and 9")
  }

  # Extract data values and remove missing values
  data_values <- shapefile[[data_col]]
  valid_rows <- !is.na(data_values)

  if (sum(valid_rows) < 2) {
    stop("Need at least 2 valid data values")
  }

  # Filter shapefile to valid rows
  shapefile_clean <- shapefile[valid_rows, ]
  data_values_clean <- data_values[valid_rows]

  # Handle region names
  if (!is.null(name_col)) {
    region_names <- shapefile_clean[[name_col]]
  } else {
    region_names <- paste("Region", 1:nrow(shapefile_clean))
  }

  # Data classification using classInt package with error handling
  n_unique_values <- length(unique(data_values_clean))

  original_n_classes <- n_classes
  if (n_classes >= n_unique_values) {
    n_classes <- max(2, n_unique_values - 1)
    message(paste("Reducing number of classes to", n_classes, "due to limited unique values"))
  }

  breaks <- tryCatch({
    if (classification == "jenks") {
      classInt::classIntervals(data_values_clean, n = n_classes, style = "jenks")$brks
    } else if (classification == "quantile") {
      classInt::classIntervals(data_values_clean, n = n_classes, style = "quantile")$brks
    } else if (classification == "equal") {
      classInt::classIntervals(data_values_clean, n = n_classes, style = "equal")$brks
    } else if (classification == "sd") {
      classInt::classIntervals(data_values_clean, n = n_classes, style = "sd")$brks
    }
  }, error = function(e) {
    # Fallback
    message(paste("Classification method", classification, "failed, falling back to equal intervals"))
    classInt::classIntervals(data_values_clean, n = n_classes, style = "equal")$brks
  })

  breaks <- unique(breaks)
  if (length(breaks) < 2) {
    # Notfall-Fallback: Erstelle manuelle Breaks
    min_val <- min(data_values_clean, na.rm = TRUE)
    max_val <- max(data_values_clean, na.rm = TRUE)
    breaks <- seq(min_val, max_val, length.out = n_classes + 1)
  }

  # Ensure breaks cover the full data range
  breaks[1] <- min(data_values_clean, na.rm = TRUE)
  breaks[length(breaks)] <- max(data_values_clean, na.rm = TRUE)

  if (length(unique(breaks)) != length(breaks)) {
    # Finale Lösung: Verwende Quantile für robuste Breaks
    breaks <- stats::quantile(data_values_clean,
                              probs = seq(0, 1, length.out = n_classes + 1),
                              na.rm = TRUE)
    breaks <- as.numeric(breaks)
    breaks <- unique(breaks)
  }

  if (length(breaks) < 2) {
    breaks <- c(min(data_values_clean, na.rm = TRUE), max(data_values_clean, na.rm = TRUE))
  }

  actual_n_classes <- length(breaks) - 1

  if (palette %in% c("viridis", "plasma", "inferno", "magma")) {
    colors <- switch(palette,
                     "viridis" = viridisLite::viridis(actual_n_classes),
                     "plasma" = viridisLite::plasma(actual_n_classes),
                     "inferno" = viridisLite::inferno(actual_n_classes),
                     "magma" = viridisLite::magma(actual_n_classes))
  } else {
    # Use RColorBrewer palettes
    colors <- RColorBrewer::brewer.pal(min(actual_n_classes, 9), palette)
    # If we need more colors than available, interpolate
    if (actual_n_classes > length(colors)) {
      colors <- grDevices::colorRampPalette(colors)(actual_n_classes)
    }
  }

  if (reverse_palette) {
    colors <- rev(colors)
  }

  # Assign each value to a class
  data_classes <- cut(data_values_clean, breaks = breaks, include.lowest = TRUE, right = FALSE)

  # Add color information to shapefile
  shapefile_clean$map_color <- colors[as.numeric(data_classes)]
  shapefile_clean$map_class <- as.character(data_classes)

  # Create popup content
  popup_content <- paste0("<b>", region_names, "</b><br>",
                          data_col, ": ", format(data_values_clean, big.mark = ",", digits = 2))

  # Add additional popup fields if specified
  if (!is.null(popup_fields)) {
    for (field in popup_fields) {
      if (field %in% names(shapefile_clean)) {
        popup_content <- paste0(popup_content, "<br>",
                                field, ": ", shapefile_clean[[field]])
      }
    }
  }

  # Create the map
  map <- leaflet::leaflet(shapefile_clean) %>%
    leaflet::addTiles() %>%
    leaflet::addPolygons(
      fillColor = ~map_color,
      weight = 1,
      opacity = 1,
      color = "white",
      dashArray = "2",
      fillOpacity = 0.7,
      popup = popup_content,
      highlight = leaflet::highlightOptions(
        weight = 3,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.9,
        bringToFront = TRUE
      )
    )

  # Create legend labels
  legend_labels <- character(actual_n_classes)
  for (i in 1:actual_n_classes) {
    legend_labels[i] <- paste0(format(breaks[i], digits = 3), " - ",
                               format(breaks[i+1], digits = 3))
  }

  # Add legend
  map <- map %>%
    leaflet::addLegend(
      pal = leaflet::colorFactor(colors, domain = data_classes),
      values = data_classes,
      opacity = 0.7,
      title = data_col,
      position = "bottomright"
    )

  # Add title
  map <- htmlwidgets::prependContent(map,
                                     htmltools::tags$h3(paste("Choropleth Map:", data_col),
                                                        style = "text-align: center; margin-bottom: 10px;"))

  # Create classification summary
  class_counts <- as.numeric(table(data_classes))
  classification_summary <- data.frame(
    class = 1:actual_n_classes,
    range = legend_labels,
    count = class_counts,
    color = colors,
    stringsAsFactors = FALSE
  )

  # Data summary statistics
  data_summary <- list(
    min = min(data_values_clean, na.rm = TRUE),
    max = max(data_values_clean, na.rm = TRUE),
    mean = mean(data_values_clean, na.rm = TRUE),
    median = stats::median(data_values_clean, na.rm = TRUE),
    sd = stats::sd(data_values_clean, na.rm = TRUE),
    missing_values = sum(is.na(data_values)),
    total_regions = length(data_values)
  )

  # Settings used
  settings <- list(
    classification = classification,
    n_classes_requested = original_n_classes,
    n_classes_actual = actual_n_classes,
    palette = palette,
    reverse_palette = reverse_palette,
    data_column = data_col,
    name_column = name_col
  )

  # Return results
  return(list(
    map = map,
    classification = classification_summary,
    data_summary = data_summary,
    settings = settings
  ))
}
