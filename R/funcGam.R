#' General additive model with functional response
#'
#' A function that takes the output of a metric calculation as done by
#' `calcMetricPerFov`. The data has to be prepared into the correct format for
#' the functional analysis by the `prepData` function. The output is a `pffr`
#' object as implemented by `refund`.
#'
#' @param data a dataframe with the following columns: Y = functional response;
#' sample_id = sample ID; image_id = image ID;
#' @param x the x-axis values of the functional response
#' @param designmat a design matrix as defined by model.matrix()
#' @param weights weights as the number of points per image. These weights are
#' normalised by the mean of the weights in the fitting process
#' @param formula the formula for the model. The colnames of the designmatrix
#' have to correspond to the variables in the formula.
#' @param family the distributional family as implemented in `family.mgcv`. For
#' fast computation the default is set to `gaussian` with a log link.
#' other interesting options can be `betar` and `scat`
#'  - for more information see `family.mgcv`.
#' @param algorithm algorithm to fit the refund::pffr method. defaults to `bam`
#' @param sandwich string indicating how and if to adjust for heterscedasticity of the 
#' residuals with a sandwich correction
#' @param bs.yindex a list specifying the spline bases for the index. See `refund::pffr`
#' for more details
#' @param bs.int a list specfying the spline bases for the global function intercept.
#'  See `refund::pffr` for more details
#' @param ... Other parameters passed to `pffr`
#'
#' @return a fitted pffr object which inherits from gam
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
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
#' summary(mdl)

#' @import dplyr
#' @importFrom methods is
#' @importFrom stats terms
functionalGam <- function(data, 
    x, 
    designmat, 
    weights, 
    formula,
    family = stats::gaussian(link = "log"),
    algorithm = "bam", 
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 20, m = c(2, 1)),
    sandwich = "cluster",
    ...) {
    # type checking
    stopifnot(is(data, "data.frame"))
    stopifnot(is(x, "vector"))
    stopifnot(is(designmat, "matrix"))
    stopifnot(is(weights, "integer") || is(weights, "numeric"))
    stopifnot(is(formula, "formula"))
    stopifnot(is(family, "character") || is(family, "family"))

    data <- cbind(data, designmat)
    # Test that length of number of points and nrow of designmat correspond
    stopifnot(length(weights) == nrow(designmat))
    # normalise the weights
    weights <- weights / mean(weights)
    # Test that the colnames of the designmat correspond to formula
    # with the exception of the intercept as this is inferred - there are
    # functional and constant intercepts
    # deselect "Intercept"
    colNames <- colnames(designmat)[!colnames(designmat) %in% c("(Intercept)",
                                                                "Intercept")]
    # stop if the rest terms of design mat are not in the formula
    stopifnot(colNames %in% attr(terms(formula), "term.labels"))
    mdl <- refund::pffr(formula,
        yind = x,
        data = data,
        weights = weights,
        family = family,
        algorithm = algorithm,
        bs.yindex = bs.yindex,
        bs.int = bs.int,
        sandwich = sandwich,
        ...
    )
    #add functional GAM class for custom printing
    class(mdl) <- c("functionalGam", class(mdl))
    return(mdl)
}

#' Summary for functionalGam object
#'
#' @param object a fitted \code{functionalGam}-object
#' @param ... see \code{\link[mgcv]{summary.gam}()} for options.
#'
#' @return A list with summary information, see \code{\link[mgcv]{summary.gam}()}
#' @export
#' @method summary functionalGam
#' @importFrom mgcv summary.gam
#' @author Martin Emons, adapted from \code{\link[refund]{summary.pffr}()} by Fabian Scheipl
summary.functionalGam <- function(object, ...){
    ret <- NextMethod("summary.pffr")
    if(any(ret$s.table[,"p-value"] < 0) | any(ret$s.table[,"p-value"] > 1)){
        warning("p-values outside of [0,1]. Clipped to [0,1]")
        ret$s.table[, "p-value"] <- pmin(pmax(ret$s.table[, "p-value"], 0), 1)
    }
    return(ret)
}
