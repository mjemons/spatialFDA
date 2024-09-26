#' Convert SpatialExperiment object to ppp object
#'
#' @param df A dataframe with the x and y coordinates from the corresponding
#' SpatialExperiment and the ColData
#' @param marks A vector of marks to be associated with the points, has to be
#' either named 'cell_type' if you want to compare discrete celltypes or else
#' continous gene expression measurements are assumed as marks.
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param window An observation window of the point pattern of class `owin`.
#' @return A ppp object for use with `spatstat` functions
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe_sub <- subset(spe, , image_number == "138")
#' df_sub <- .speToDf(spe_sub)
#' pp <- .dfToppp(df_sub, marks = "cell_type")
#'
#' @importFrom SummarizedExperiment colData
.dfToppp <- function(df, marks = NULL, continuous = FALSE, window = NULL) {
    # this definition of the window is quite conservative - can be set explicitly
    pp <- spatstat.geom::as.ppp(data.frame(x = df$x, y = df$y),
        W = spatstat.geom::owin(
            c(
                min(df$x) - 1,
                max(df$x) + 1
            ),
            c(
                min(df$y) - 1,
                max(df$y) + 1
            )
        )
    )
    # set the marks
    if (!continuous) {
        spatstat.geom::marks(pp) <- factor(df[[marks]])
    } else {
        spatstat.geom::marks(pp) <- subset(df, select = names(df) %in% marks)
    }
    # if window exist, set is as new window and potentially exclude some points
    if (!is.null(window)) {
        pp <- spatstat.geom::as.ppp(spatstat.geom::superimpose(pp, W = window))
    }

    return(pp)
}

#' Transform a SpatialExperiment into a dataframe
#'
#' @param spe A SpatialExperiment object subset to a single image
#'
#' @return A dataframe with the x and y coordinates from the corresponding
#' SpatialExperiment and the ColData
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe_sub <- subset(spe, , image_number == "138")
#' df_sub <- .speToDf(spe_sub)
.speToDf <- function(spe) {
    df <- data.frame(
        x = SpatialExperiment::spatialCoords(spe)[, 1],
        y = SpatialExperiment::spatialCoords(spe)[, 2]
    )
    df <- cbind(df, colData(spe))
}
