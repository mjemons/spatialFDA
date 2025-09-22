#' Heuristic for the choice of rMax
#' @param spe a `SpatialExperiment` object
#' @param subsetby the spe `colData` variable to subset the data by. This
#' variable has to be provided, even if there is only one sample.
#' @param marks the marks to consider e.g. cell types
#'
#' @returns a ggplot histogram of the bounding radius of all the
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' p <- rMaxHeuristic(spe,
#' subsetby = "image_number", marks = "cell_type"
#' )
rMaxHeuristic <- function(spe, subsetby, marks){
  #create a dataframe
  df <- .speToDf(spe)
  #split the dataframe
  dfLs <- base::split(df, df[[subsetby]])
  radiusDf <- lapply(dfLs, function(dfSub){
    pp <- .dfToppp(dfSub, marks = marks)
    return(data.frame("bounding_radius" = spatstat.geom::boundingradius(pp)))
  }) %>% dplyr::bind_rows()
 rMax <- min(radiusDf[["bounding_radius"]]) / 2
 p <- ggplot(radiusDf, aes(.data[["bounding_radius"]])) +
   geom_histogram(binwidth = (10), boundary=0) +
   theme_light() +
   geom_vline(xintercept = rMax, color = "red") +
   xlab("minimum bounding radius") +
   ggtitle(paste0("heuristic rMax is ", rMax))
 return(p)
}
