#' Plot a spatial metric per field of view
#'
#' @param spe a spatial experiment object
#' @param selection the mark(s) you want to compare
#' @param fun the spatstat function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param r_seq the range of r values to compute the function over
#' @param ncores the number of cores to use for parallel processing, default is 1
#' @param correction the border correction to apply
#' @param x the x-axis variable
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'  spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#'  p <- plotMetricPerFov(spe, c('alpha', 'beta'), fun = 'Gcross'
#'  , marks = 'cell_type', r_seq = seq(0,50, length.out = 50), ncores = 2,
#'  correction = 'rs', x = 'r')
#'  print(p)
#' @import dplyr ggplot2
plotMetricPerFov <- function(spe, selection, fun, marks = NULL, r_seq = NULL, ncores = 1, correction = NULL, x = NULL){
  metric_df <- calcMetricPerFov(spe, selection, fun, marks, r_seq, ncores = ncores)
  metric_df$ID <- paste0(metric_df$condition,'|' ,metric_df$patient_id)
  p <- ggplot(metric_df, aes(x = .data[[x]], y = .data[[correction]], color = factor(image_id))) +
    geom_line() +
    facet_wrap(~ID) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste0(fun, ' metric for ', paste(selection, collapse = ' and ')))
  return(p)
}
