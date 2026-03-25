#' Plot a pffr model object
#'
#' A function that takes a pffr object as calculated in `functionalGam` and
#' plots the functional coefficients. The functions are centered such that their
#' expected value is zero. Therefore, the scalar intercept has to be added to
#' the output with the argument `shift` in order to plot the coefficients in
#' their original range.
#'
#' @param mdl a `pffr` model object
#' @param predictor predictor to plot
#' @param shift the value by which to shift the centered functional intercept.
#' this will most often be the constant intercept
#'
#' @return ggplot object of the functional estimate
#' @export
#'
#' @examples
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
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
#' # create a unique ID for each row
#' metricRes$ID <- paste0(
#'     metricRes$patient_stage, "x", metricRes$patient_id,
#'     "x", metricRes$image_number
#' )
#'
#' dat <- prepData(metricRes, "r", "rs", sample_id = "patient_id",
#'     image_id = "image_number", condition = "patient_stage")
#'
#'#' # drop rows with NA
#' dat <- dat |> drop_na()
#'
#' # create a designmatrix
#' condition <- dat$patient_stage
#' # relevel the condition - can set explicit contrasts here
#' condition <- relevel(condition, "Non-diabetic")
#' designmat <- model.matrix(~condition)
#' # colnames don't work with the '-' sign
#' colnames(designmat) <- c(
#'     "(Intercept)", "conditionLong_duration",
#'     "conditionOnset"
#' )
#' # fit the model
#' mdl <- functionalGam(
#'     data = dat, x = metricRes$r |> unique(),
#'     designmat = designmat, weights = dat$npoints,
#'     formula = formula(Y ~ conditionLong_duration +
#'         conditionOnset + s(patient_id, bs = "re")),
#'         algorithm = "bam"
#' )
#' summary(mdl, re.test = FALSE)
#' plotLs <- lapply(colnames(designmat), plotMdl,
#'     mdl = mdl,
#'     shift = mdl$coefficients[["(Intercept)"]]
#' )
#' @import dplyr
#' @importFrom methods is
plotMdl <- function(mdl, predictor, shift = NULL) {
    # type checking
    stopifnot(is(mdl, "pffr"))
    stopifnot(is(predictor, "character"))
    # extract the coefficients from the model
    coef <- coef(mdl)
    if (predictor == "(Intercept)" && !is.null(shift)) {
        #rename as pffr output is without brackets
        predictor = "Intercept"
        coef$smterms[["Intercept(x)"]]$coef$value <-
          coef$smterms[["Intercept(x)"]]$coef$value + shift
    }
    # get the actual values into a dataframe
    df <- coef$smterms[[paste0(predictor, "(x)")]]$coef
    # plot
    p <- ggplot(df, aes(.data$yindex.vec, .data$value)) +
        geom_line(linewidth = 1) +
        # here, I implement a Wald CI - could be improved
        geom_ribbon(
            data = df,
            aes(ymin = .data$value - 1.96 * .data$se,
                ymax = .data$value + 1.96 * .data$se),
            alpha = 0.3
        ) +
        geom_hline(
            yintercept = 0,
            linetype = "dashed", color = "red", linewidth = 1
        ) +
        ggtitle(predictor) +
        ylab("parameter value") +
        xlab("r") +
        theme_light()
    return(p)
}
