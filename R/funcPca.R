#' Functional Principal Component Analysis
#'
#' A function that takes as input the output of `calcMetricPerFov` which has to
#' be converted into the correct format by `prepData`. The output is a list with
#' the `fpca.face` output from refund.
#'
#' @param data a data object for functional data analysis containing at least
#' the functional response $Y$.
#' @param r the functional domain
#' @param ... Other parameters passed to `fpca.sc` functions
#'
#' @return a list with components of fpca.sc
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' # calculate the Gcross metric for alpha and Tc cells
#' metricRes <- calcMetricPerFov(spe, c("alpha", "Tc"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id", "image_number"), ncores = 1
#' )
#' metricRes$ID <- paste0(
#'   metricRes$patient_stage, "|", metricRes$patient_id,
#'   "|", metricRes$image_number
#' )
#' # prepare data for FDA
#' dat <- prepData(metricRes, "r", "rs")
#'
#' # drop rows with NA
#' dat <- dat |> drop_na()
#' # create meta info of the IDs
#' splitData <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(splitData, function(x) x[1]))
#' dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
#' dat$image_id <- factor(sapply(splitData, function(x) x[3]))
#' # calculate fPCA
#' mdl <- functionalPCA(
#'     data = dat, r = metricRes$r |> unique()
#' )
#' @import dplyr
#' @importFrom methods is
functionalPCA <- function(data, r, ...) {
    stopifnot(is(data, "data.frame"))
    stopifnot(is(r, "vector"))
    # calculate the fPCA - this is a bit a pointless wrapper until now
    res <- refund::fpca.sc(
        Y = data$Y, center = TRUE, argvals = r, ...
    )
    return(res)
}

#' Plot a biplot from an fPCA analysis
#'
#' A function that takes the output from the `functionalPCA` function and
#' returns a `ggplot` object of the first two dimensions of the PCA as biplot.
#'
#' @param data a data object for functional data analysis containing at least
#' the functional response $Y$.
#' @param res the output from the fPCA calculation
#' @param colourby the variable by which to colour the PCA plot by
#' @param labelby the variable by which to label the PCA plot by
#'
#' @return a list with components of fpca.face
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' # calculate the Gcross metric for alpha and beta cells
#' metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id", "image_number"), ncores = 1
#' )
#' metricRes$ID <- paste0(
#'   metricRes$patient_stage, "|", metricRes$patient_id,
#'   "|", metricRes$image_number
#' )
#'
#' # prepare data for FDA
#' dat <- prepData(metricRes, "r", "rs")
#'
#' # drop rows with NA
#' dat <- dat |> drop_na()
#' # create meta info of the IDs
#' splitData <- str_split(dat$ID, "|")
#' dat$condition <- factor(sapply(splitData, function(x) x[1]))
#' dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
#' dat$image_id <- factor(sapply(splitData, function(x) x[3]))
#' # calculate fPCA
#' mdl <- functionalPCA(
#'     data = dat, r = metricRes$r |> unique()
#' )
#' p <- plotFpca(
#'     data = dat, res = mdl, colourby = "condition",
#'     labelby = "patient_id"
#' )
#' print(p)
#' @import dplyr
#' @importFrom methods is
plotFpca <- function(data, res, colourby = NULL, labelby = NULL) {
    stopifnot(is(data, "data.frame"))
    stopifnot(is(res, "fpca"))
    scoresDf <- res$scores %>% as.data.frame()
    # plot fCPA results - assumes same order of fPCA results and input data
    p <- ggplot(scoresDf, aes(scoresDf[, 1], scoresDf[, 2],
        colour = factor(data[[colourby]])
    )) +
        geom_point() +
        coord_equal() +
        theme_light() +
        xlab("functional PC1") +
        ylab("functional PC2")
    if (!is.null(labelby)) {
        p <- p +
            geom_text(hjust = 0, vjust = 0,
                      aes(label = factor(data[[labelby]])))
    }
    return(p)
}

#' print the fPCA results
#'
#' this is a function that prints a summary of the fPCA result of class `fpca`
#'
#' @param x the result of function `functionalPCA`
#' @param ... other parameters passed to base generic function `print`
#'
#' @return a formatted overview of the fPCA result
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' # calculate the Gcross metric for alpha and beta cells
#' metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id", "image_number"), ncores = 1
#' )
#' metricRes$ID <- paste0(
#'   metricRes$patient_stage, "|", metricRes$patient_id,
#'   "|", metricRes$image_number
#' )
#' # prepare data for FDA
#' dat <- prepData(metricRes, "r", "rs")
#'
#' # drop rows with NA
#' dat <- dat |> drop_na()
#'
#' # create meta info of the IDs
#' splitData <- strsplit(dat$ID, "|", fixed = TRUE)
#' dat$condition <- factor(sapply(splitData, function(x) x[1]))
#' dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
#' dat$image_id <- factor(sapply(splitData, function(x) x[3]))
#' # calculate fPCA
#' mdl <- functionalPCA(
#'     data = dat, r = metricRes$r |> unique()
#' )
#' mdl
print.fpca <- function(x, ...) {
    cat(paste0("Functional principal components analysis ",
                 "result from `fpca.sc` of `refund` \n"))
    cat(paste0("class: ", class(x), "\n"))
    cat(paste0("pve: ", round(x$pve, digits = 4), "\n"))
    cat(paste0("number of PCs: ", x$npc, "\n"))
    cat(paste0("eigenvalues:\n"))
    cat(x$evalues)
}
