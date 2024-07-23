#' Plot a spatial metric per field of view
#'
#' @param metric_df the metric dataframe as calculated by calcMetricPerFov
#' @param theo logical; if the theoretical line should be plotted
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#' @param image_id the ID of the image/fov
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     r_seq = seq(0, 50, length.out = 50), by = c("patient_stage", "patient_id"),
#'     ncores = 2
#' )
#' p <- plotMetricPerFov(metric_res, correction = "rs", x = "r", image_id = 'image_id', ID = 'ID')
#' print(p)
#' @import dplyr ggplot2
plotMetricPerFov <- function(metric_df, theo = FALSE, correction = NULL, x = NULL, image_id = NULL, ID = NULL) {
  p <- ggplot(metric_df, aes(x = .data[[x]], y = .data[[correction]], group = factor(.data[[image_id]]), colour = factor(.data[[ID]]))) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste0(metric_df$fun, " metric for ", unique(metric_df$selection)))
  if (theo == TRUE) p <- p + geom_line(aes(x=.data[[x]],y=theo),linetype = "dashed", color = "black")
  if (!is.null(ID)) p <- p + facet_wrap(ID)
  return(p)
}

#' Creates a nXn plot of the cross metrics per sample
#'
#' @param sub_fov a subset of the dataframe to the respective fov
#' @param theo logical; if the theoretical line should be plotted
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#' @param image_id the ID of the image/fov
#'
#' @return a ggplot object
#' @export
#'
plotCrossFOV <- function(sub_fov, theo, correction, x, image_id, ID) {
  #  Apply plot metric function for each combination
  lp <- lapply(unique(sub_fov$selection), function(sel)
    plotMetricPerFov(sub_fov[sub_fov$selection == sel, ], theo, correction, x, image_id, ID))
  #  Count number of marks
  nMarks <- length(unique(sub_fov$selection))
  # Wraps the plot in an nXn grid
  p <- wrap_plots(lp, ncol = sqrt(nMarks)) +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom')
  return(p)
}


#' Plot a cross type spatial metric per field of view
#'
#' @param metric_df the metric dataframe as calculated by calcMetricPerFov
#' @param theo logical; if the theoretical line should be plotted
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#' @param image_id the ID of the image/fov
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metric_res <- calcCrossMetricPerFov(spe, c("alpha", "beta", "delta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     r_seq = seq(0, 50, length.out = 50), by = c("patient_stage", "patient_id"),
#'     ncores = 2
#' )
#' metric_res <- subset(metric_res, image_id %in% c(138,139,140))
#' p <- plotCrossMetricPerFov(metric_res, theo = TRUE, correction = "rs", x = "r", image_id = 'image_id', ID = 'ID')
#' print(p)
plotCrossMetricPerFov <- function(metric_df,
                                  theo = NULL,
                                  correction = NULL,
                                  x = NULL,
                                  image_id = NULL,
                                  ID = NULL) {
  # Find all unique samples
  samples <- metric_df[[image_id]] |> unique()

  # Applies the function abouve to all samples
  res_p <- lapply(samples, function(fov) {
    sub_fov <- subset(metric_df, image_id %in% fov)
    return(plotCrossFOV(sub_fov, theo, correction, x, image_id))
  })

  return(res_p)
}
