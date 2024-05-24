#' Plot a spatial metric per field of view
#'
#' @param metric_df the metric dataframe as calculated by calcMetricPerFov
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'  spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#'  metric_res <- calcMetricPerFov(spe, c('alpha', 'beta'),
#'  subsetby = 'image_number',fun = 'Gcross', marks = 'cell_type',
#'  r_seq = seq(0,50, length.out = 50), by = c('patient_stage', 'patient_id'),
#'  ncores = 2)
#'  p <- plotMetricPerFov(metric_res, correction = 'rs', x = 'r')
#'  print(p)
#' @import dplyr ggplot2
plotMetricPerFov <- function(metric_df, correction = NULL, x = NULL){
  p <- ggplot(metric_df, aes(x = .data[[x]], y = .data[[correction]], group = factor(image_id), colour = factor(ID))) +
    geom_line() +
    #geom_line(aes(x=.data[[x]],y=theo),linetype = "dashed")+
    facet_wrap(~ID) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste0(metric_df$fun, ' metric for ', unique(metric_df$selection)))
  return(p)
}
