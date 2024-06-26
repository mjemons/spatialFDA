#' Functional Principal Component Analysis
#'
#' @param dat a data object for functional data analysis containing at least the
#' the functional
#' @param r the functional domain
#' @param knots the number of knots
#' @param pve the proportion of variance explained
#'
#' @return a list with components of fpca.face
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
#'     c("patient_stage", "patient_id"), ncores = 2
#' )
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_id
#' )
#' # extract the functional response matrix
#' mat <- metric_res %>%
#'     select(ID, r, rs) %>%
#'     spread(ID, rs) %>%
#'     select(!r)
#' # create a dataframe as required by pffr
#' dat <- data.frame(ID = colnames(mat))
#' dat$Y <- t(mat)
#' # create meta info of the IDs
#' split_data <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(split_data, `[`, 1))
#' dat$patient_id <- factor(sapply(split_data, `[`, 2))
#' dat$image_id <- factor(sapply(split_data, `[`, 3))
#' # calculate fPCA
#' mdl <- functionalPCA(dat = dat, r = metric_res$r |> unique(), knots = 30, pve = 0.99)
#' @import dplyr
functionalPCA <- function(dat, r, knots, pve = 0.95) {
    # calculate the fPCA - this is a bit a pointless wrapper until now
    res <- refund::fpca.face(Y = dat$Y, center = TRUE, argvals = r, knots = knots, pve = pve)
    return(res)
}

#' Plot a biplot from an fPCA analysis
#'
#' @param dat a data object for functional data analysis containing at least the
#' the functional
#' @param res the output from the fPCA calculation
#'
#' @return a list with components of fpca.face
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
#'     c("patient_stage", "patient_id"), ncores = 2
#' )
#' metric_res$ID <- paste0(
#'     metric_res$patient_stage, "x", metric_res$patient_id,
#'     "x", metric_res$image_id
#' )
#' # extract the functional response matrix
#' mat <- metric_res %>%
#'     select(ID, r, rs) %>%
#'     spread(ID, rs) %>%
#'     select(!r)
#' # create a dataframe as required by pffr
#' dat <- data.frame(ID = colnames(mat))
#' dat$Y <- t(mat)
#' # create meta info of the IDs
#' split_data <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(split_data, `[`, 1))
#' dat$patient_id <- factor(sapply(split_data, `[`, 2))
#' dat$image_id <- factor(sapply(split_data, `[`, 3))
#' # calculate fPCA
#' mdl <- functionalPCA(dat = dat, r = metric_res$r |> unique(), knots = 30, pve = 0.99)
#' p <- plotFpca(dat = dat, res = mdl)
#' print(p)
#' @import dplyr
plotFpca <- function(dat, res) {
    scores_df <- res$scores %>% as.data.frame()
    # plot fCPA results - assumes same order of fPCA results and input data
    p <- ggplot(scores_df, aes(scores_df[, 1], scores_df[, 2], colour = factor(dat$condition), label = factor(dat$patient_id))) +
        geom_point() +
        geom_text(hjust = 0, vjust = 0) +
        coord_equal() +
        theme_light()
    return(p)
}
