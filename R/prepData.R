#' Prepare data from calcMetricRes to be in the right format for FDA
#'
#' @param metricRes a dataframe as calculated by calcMetricRes - requires
#' the columns ID (unique identifier of each row)
#' @param x the name of the x-axis of the spatial metric
#' @param y the name of the y-axis of the spatial metric
#' @param sample_id the spe `colData` variable to mark the sample
#' @param image_id the spe `colData` variable to mark the image
#' @param condition the spe `colData` variable to mark the condition
#' @param ID the unique ID corresponding to each row in the metricRes defaults
#' to `"ID"`
#'
#' @return returns a list with three entries, the unique ID, the functional
#' response Y and the weights
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' metricRes <- calcMetricPerFov(spe, c("alpha", "Tc"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), by = c(
#'         "patient_stage", "patient_id",
#'         "image_number"
#'     ),
#'     ncores = 1
#' )
#'
#' # create a unique ID for each row
#' metricRes$ID <- paste0(
#'     metricRes$patient_stage, "|", metricRes$patient_id,
#'     "|", metricRes$image_number
#' )
#' dat <- prepData(metricRes, "r", "rs", sample_id = "patient_id",
#' image_id = "image_number", condition = "patient_stage")
#' @import tidyr
#' @importFrom methods is
prepData <- function(metricRes, x, y, sample_id = NULL, image_id = NULL,
                     condition = NULL, ID = "ID"){
    # type checking
    stopifnot(
        is.data.frame(metricRes),
        is.character(x),
        is.character(y)
    )
    # extract the functional response matrix
    mat <- metricRes %>%
        dplyr::select(ID, x, y) %>%
        tidyr::spread(ID, y) %>%
        dplyr::select(!x) %>%
        as.matrix() 


    storage.mode(mat) <- "numeric"
    # create a dataframe as required by pffr
    # the colnames of the matrix are the new row IDs
    dat <- data.frame(ID = colnames(mat))
    # transpose of the matrix to have the entire response in one row
    dat$Y <- t(mat)
    # extract the number of points as weights
    weights <- metricRes %>%
      dplyr::select("ID", "npoints", "npointsmin", "npointsmax") %>%
      unique()
    # add the weights to the data.frame
    dat <- dat %>% dplyr::left_join(weights, by = "ID")
    # extract the coordinates and meta data and convert to factor
    meta <- metricRes %>%
        dplyr::select("ID", "centroidx", "centroidy", "minIntensity", sample_id,
                      image_id, condition, "intensityQuery") %>%
        unique() %>%
        mutate(across(c(sample_id,
                        image_id, condition), as.factor))
    # add the coordinates to the data.frame
    dat <- dat %>% dplyr::left_join(meta, by = "ID")

    return(dat)
}
