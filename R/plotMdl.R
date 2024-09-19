#' Plot a pffr model object
#'
#' @param mdl a pffr model object
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
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     r_seq = seq(0, 50, length.out = 50), by = c("patient_stage", "patient_id",
#'     "image_number"),
#'     ncores = 1
#' )
#' # create a unique ID for each row
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_number
#' )
#'
#' dat <- prepData(metric_res, "r", "rs")
#'
#' # create meta info of the IDs
#' split_data <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(split_data, function(x) x[1]))
#' dat$patient_id <- factor(sapply(split_data, function(x) x[2]))
#' dat$image_id <- factor(sapply(split_data, function(x) x[3]))
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
#' plot_ls <- lapply(colnames(designmat), plotMdl, mdl = mdl, shift = mdl$coefficients[['(Intercept)']])
#' @import dplyr
plotMdl <- function(mdl, predictor, shift = NULL) {
    # extract the coefficients from the model
    coef <- coef(mdl)
    if(predictor == 'Intercept' && !is.null(shift)){
      coef$sm[["Intercept(x)"]]$coef$value <- coef$sm[["Intercept(x)"]]$coef$value + shift
    }
    # get the actual values into a dataframe
    df <- coef$sm[[paste0(predictor, "(x)")]]$coef
    # plot
    p <- ggplot(df, aes(.data$x.vec, .data$value)) +
        geom_line(size = 1) +
        # here, I implement a Wald CI - could be g improved
        geom_ribbon(data = df, aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se), alpha = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", size=1) +
        ggtitle(predictor) +
        ylab('parameter value') +
        xlab('r') +
        theme_light()
    return(p)
}
