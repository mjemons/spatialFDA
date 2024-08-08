#' General additive model with functional response
#'
#' @param data a dataframe with the following columns: Y = functional response;
#' sample_id = sample ID; image_id = image ID;
#' @param x the x-axis values of the functional response
#' @param designmat a design matrix as defined by model.matrix()
#' @param weights weights as the number of points per image. These weights are
#' normalised by the mean of the weights in the fitting process
#' @param formula the formula for the model. The colnames of the designmatrix
#' have to correspond to the variables in the formula.
#'
#' @return a fitted pffr object which inherits from gam
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' # calculate the Gcross metric for alpha and beta cells
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", r_seq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id", "image_number"), ncores = 1
#' )
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_number
#' )
#' # prepare data for FDA
#' dat <- prepData(metric_res, "r", "rs")
#'
#' # drop rows with NA
#' dat <- dat |> drop_na()
#'
#' # create meta info of the IDs
#' split_data <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(split_data, `[`, 1))
#' dat$patient_id <- factor(sapply(split_data, `[`, 2))
#' dat$image_id <- factor(sapply(split_data, `[`, 3))
#' # create a designmatrix
#' condition <- dat$condition
#' # relevel the condition - can set explicit contrasts here
#' condition <- relevel(condition, "Non-diabetic")
#' designmat <- model.matrix(~condition)
#' # colnames don't work with the '-' sign
#' colnames(designmat) <- c("Intercept", "conditionLong_duration", "conditionOnset")
#' # fit the model
#' mdl <- functionalGam(
#'     dat = dat, x = metric_res$r |> unique(),
#'     designmat = designmat, weights = dat$weights$npoints,
#'     formula = formula(Y ~ conditionLong_duration +
#'         conditionOnset + s(patient_id, bs = "re"))
#' )
#' summary(mdl)
#' plot_ls <- lapply(colnames(designmat), plotMdl, mdl = mdl)

#' @import dplyr

functionalGam <- function(data, x, designmat, weights, formula) {
    # get the colnames
    colnam <- colnames(designmat)
    for (i in 1:length(colnam)) {
        data[[colnam[i]]] <- designmat[, i]
    }
    # TODO how to make weighting optional?
    # normalise the weights
    weights <- weights / mean(weights)
    # TODO write a test that the colnames of the designmat correspond to formula
    startTime <- Sys.time()
    mdl <- refund::pffr(formula,
        yind = x,
        data = data,
        weights = weights
    )
    endTime <- Sys.time()
    print(endTime - startTime)
    return(mdl)
}
