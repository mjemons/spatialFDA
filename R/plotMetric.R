#' Plot a spatial metric per field of view
#'
#' A function that plots the output of the function `calcMetricPerFov`. The plot
#' contains one curve per FOV and makes subplots by samples.
#'
#' @param metricDf the metric `dataframe` as calculated by `calcMetricPerFov`
#' @param theo logical; if the theoretical line should be plotted
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#' @param imageId the ID of the image/fov
#' @param ID the (optional) ID for plotting combinations
#' @param nrow the number of rows for the facet wrap
#' @param ncol the number of columns for the facet wrap
#' @param legend.position the position of the legend of the plot
#' @param ... Other parameters passed to `ggplot2` functions
#'
#' @return a `ggplot` object
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
#' # ceate a unique plotting ID
#' metricRes$ID <- paste0(
#'  metricRes$patient_stage, "|", metricRes$patient_id
#' )
#'
#' p <- plotMetricPerFov(metricRes,
#'     correction = "rs", x = "r",
#'     imageId = "image_number", ID = "ID"
#' )
#' print(p)
#' @import dplyr ggplot2
#' @importFrom methods is
plotMetricPerFov <- function(metricDf, theo = FALSE, correction = NULL,
    x = NULL, imageId = NULL, ID = NULL, nrow = NULL, ncol = NULL,
    legend.position = "none", ...) {
    # type checking
    stopifnot(
        is.data.frame(metricDf),
        is.character(correction),
        is.character(x)
    )
    p <- ggplot(metricDf, aes(
        x = .data[[x]], y = .data[[correction]],
        group = factor(.data[[imageId]])
    ))
    if (!is.null(ID)) {
        p <- p +
            geom_line(aes(colour = factor(.data[[ID]])), ...) +
            facet_wrap(selection ~ ID, nrow, ncol)
    } else {
        p <- p +
            geom_line(aes(colour = factor(.data[[imageId]])), ...)
    }
    p <- p +
        theme_minimal() +
        theme(legend.position = legend.position) +
        labs(title = paste0(
            metricDf$fun, " metric for ",
            unique(metricDf$selection)
        ))
    if (theo == TRUE) {
        p <- p + geom_line(aes(x = .data[[x]], y = theo),
            linetype = "dashed", color = "black", ...
        )
    }
    return(p)
}

#' Creates a nXn plot of the cross metrics per sample
#'
#' Helper function for `plotCrossMetricPerFov`. It applies `plotMetricPerFov`
#' to all `n` marks defined in the variable `selection`. This gives an
#' `nxn` plot of all marks.
#'
#' @param subFov a subset of the `dataframe` to the respective fov
#' @param theo logical; if the theoretical line should be plotted
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#' @param imageId the ID of the image/fov
#' @param ID the (optional) ID for plotting combinations
#' @param nrow the number of rows for the facet wrap
#' @param ncol the number of columns for the facet wrap
#' @param legend.position the position of the legend of the plot
#' @param ... Other parameters passed to `ggplot2` functions
#'
#' @return a ggplot object
#' @export
#'
#' @importFrom methods is
plotCrossFOV <- function(subFov, theo, correction, x, imageId, ID = NULL,
                         ncol = NULL, nrow = NULL,
                         legend.position = "none", ...) {
    # type checking
    stopifnot(is.data.frame(subFov))
    #  Apply plot metric function for each combination
    lp <- lapply(unique(subFov$selection), function(sel) {
        plotMetricPerFov(
            metricDf = subFov[subFov$selection == sel, ],
            theo = theo, correction = correction, x = x,
            imageId = imageId, ID = ID, nrow = nrow, ncol = ncol,
            legend.position = legend.position, ...
        )
    })
    #  Count number of marks
    nMarks <- length(unique(subFov$selection))
    # Wraps the plot in an nXn grid
    p <- patchwork::wrap_plots(lp, ncol = sqrt(nMarks)) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    return(p)
}


#' Plot a cross type spatial metric per field of view
#'
#' This function plots the cross function between two marks output from
#' `calcMetricPerFov`. It wraps around helper function and applies this
#' function to all samples.
#'
#' @param metricDf the metric dataframe as calculated by `calcMetricPerFov`
#' @param theo logical; if the theoretical line should be plotted
#' @param correction the border correction to plot
#' @param x the x-axis variable to plot
#' @param imageId the ID of the image/fov
#' @param ID the (optional) ID for plotting combinations
#' @param nrow the number of rows for the facet wrap
#' @param ncol the number of columns for the facet wrap
#' @param legend.position the position of the legend of the plot
#' @param ... Other parameters passed to `ggplot2` functions
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' metricRes <- calcCrossMetricPerFov(spe, c("alpha", "Tc"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), by = c(
#'         "patient_stage", "patient_id",
#'         "image_number"
#'     ),
#'     ncores = 1
#' )
#'
#' metricRes$ID <- paste0(
#'  metricRes$patient_stage, "|", metricRes$patient_id
#' )
#'
#' metricRes <- subset(metricRes, image_number %in% c(138, 139, 140))
#' p <- plotCrossMetricPerFov(metricRes,
#'     theo = TRUE, correction = "rs",
#'     x = "r", imageId = "image_number", ID = "ID"
#' )
#' print(p)
#' @importFrom methods is
plotCrossMetricPerFov <- function(
        metricDf,
        theo = NULL,
        correction = NULL,
        x = NULL,
        imageId = NULL,
        ID = NULL,
        nrow = NULL,
        ncol = NULL,
        legend.position = "none",
        ...) {
    # type checking
    stopifnot(
        is.data.frame(metricDf),
        is.character(imageId)
    )
    # Find all unique samples
    samples <- metricDf[[imageId]] |> unique()

    # Applies the function abouve to all samples
    resP <- lapply(samples, function(fov) {
        subFov <- metricDf[metricDf[[imageId]] %in% fov, ]
        return(plotCrossFOV(
            subFov = subFov, theo = theo, correction = correction,
            x = x, imageId = imageId, ID = ID, nrow = nrow, ncol = ncol,
            legend.position = legend.position, ...
        ))
    })

    return(resP)
}

#' Functional boxplot of spatstat curves
#'
#' This function creates a functional boxplot of the spatial statistics curves.
#' It creates one functional boxplot per aggregation category, e.g. condition.
#'
#' @param metricDf the metric dataframe as calculated by `calcMetricPerFov`
#' @param x the name of the x-axis of the spatial metric
#' @param y the name of the y-axis of the spatial metric
#' @param aggregateBy the criterion by which to aggregate the curves into a
#' functional boxplot. Can be e.g. the condition of the different samples.
#'
#' @return a list of base R plots
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
#' # create a unique ID for the data preparation
#' metricRes$ID <- paste0(
#'   metricRes$patient_stage, "|", metricRes$patient_id,
#'   "|", metricRes$image_number
#' )
#'
#' plotFbPlot(metricRes, 'r', 'rs', 'patient_stage')
#' @importFrom fda fbplot
#' @importFrom graphics title
plotFbPlot <- function(
    metricDf, x, y, aggregateBy) {
  aggregationLs <- metricDf[[aggregateBy]] %>% unique
  ylim <- c(min(metricDf[[y]]), max(metricDf[[y]]))
  lapply(aggregationLs, function(aggregate){
      filteredData <- metricDf %>% filter(.data[[aggregateBy]] == aggregate)
      res <- prepData(filteredData, x, y) %>% drop_na
      fda::fbplot(t(res$Y), ylim = ylim)
      graphics::title(main = aggregate)
    })
}
