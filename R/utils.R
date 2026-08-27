#' load Example dataset from Damond et al. (2019)
#'
#' @param full a boolean indicating whether to load the entire
#' Damond et al. (2019) or only a subset
#'
#' @returns A SpatialExperiment object as uploaded to `ExperimentHub()`
#' @export
#'
#' @examples
#' # retrieve the Damond et al. (2019) dataset
#' spe <- .loadExample()
#'
#' @importFrom ExperimentHub ExperimentHub
.loadExample <- function(full = FALSE) {
   # retrieve data from EH directly - code adapted from `imcdatasets`
   # Damond et.al (2023) licensed under GPLv3
   # (https://www.gnu.org/licenses/gpl-3.0.txt)
   eh <- ExperimentHub::ExperimentHub()
   if(full){
     title = "Damond_2019_Pancreas - sce - v1 - full"
   }
   else{
     title = "Damond_2019_Pancreas - sce - v1"
   }
   object_id <- eh[eh$title == title]$ah_id
   sce <- eh[[object_id]]
   # rename coordinates
   colData(sce)$x <- colData(sce)$cell_x
   colData(sce)$y <- colData(sce)$cell_y
   # create SPE object
   spe <- toSpatialExperiment(sce,
                              sample_id = "image_name",
                              spatialCoordsNames = c("x", "y"))
   return(spe)
}

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
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' speSub <- subset(spe, , image_number == "138")
#' dfSub <- .speToDf(speSub)
#' pp <- .dfToppp(dfSub, marks = "cell_type")
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
.dfToppp <- function(df, marks = NULL, continuous = FALSE, window = NULL) {
    #type checking
    stopifnot(is.data.frame(df))
    # this definition of the window is quite conservative
    # - can be set explicitly
    pp <- spatstat.geom::as.ppp(data.frame(x = df$x, y = df$y),
        W = spatstat.geom::owin(
            c(
                base::min(df$x) - 1,
                base::max(df$x) + 1
            ),
            c(
                base::min(df$y) - 1,
                base::max(df$y) + 1
            )
        )
    )
    # set the marks
    if (!continuous) {
        spatstat.geom::marks(pp) <- factor(df[[marks]])
    } else {
        spatstat.geom::marks(pp) <- base::subset(df, select =
                                                   names(df) %in% marks)
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
#' SpatialExperiment and the colData
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' speSub <- subset(spe, , image_number == "138")
#' dfSub <- .speToDf(speSub)
#' @importFrom methods is
.speToDf <- function(spe) {
    stopifnot(is(spe, "SpatialExperiment"))
    df <- data.frame(
        x = SpatialExperiment::spatialCoords(spe)[, 1],
        y = SpatialExperiment::spatialCoords(spe)[, 2]
    )
    df <- cbind(df, colData(spe))
}
