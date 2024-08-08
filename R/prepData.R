#' Prepare data from calcMetricRes to be in the right format for FDA
#'
#' @param metric_res a dataframe as calculated by calcMetricRes - requires
#' the columns ID (unique identifier of each row)
#' @param x the name of the x-axis of the spatial metric
#' @param y the name of the y-axis of the spatial metric
#'
#' @return returns a list with three entries, the unique ID, the functional
#' response Y and the weights
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     r_seq = seq(0, 50, length.out = 50), by = c("patient_stage", "patient_id"
#'     , "image_number"),
#'     ncores = 1
#' )
#'
#' # create a unique ID for each row
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_number
#' )
#' dat <- prepData(metric_res, "r", "rs")
#' @import tidyr
prepData <- function(metric_res, x, y) {
    # extract the functional response matrix
    mat <- metric_res %>%
        select(ID, x, y) %>%
        spread(ID, y) %>%
        select(!x)
    # extract the number of points as weights - are ordered differently,
    # thus order according to image ID
    weights <- metric_res %>%
        group_by(ID) %>%
        select(ID, npoints) %>%
        unique() %>%
        arrange(ID)
    # create a dataframe as required by pffr
    dat <- data.frame(ID = colnames(mat))
    dat$Y <- t(mat)
    dat$weights <- weights
    return(dat)
}
