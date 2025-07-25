#' Statistical Inference on Spatial Statistics Functions
#'
#' A function to perform spatial statistical inference on spatial omics data.
#' This function works so far only on functions of radius "r".
#'
#' @param spe a `SpatialExperiment` object
#' @param selection the mark(s) you want to compare. NOTE: This is directional.
#' c(A,B) is not the same result as c(B,A).
#' @param subsetby the spe `colData` variable to subset the data by. This
#' variable has to be provided, even if there is only one sample.
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param correction the edge correction to be applied
#' @param sample_id the spe `colData` variable to mark the sample, if not NULL
#' this will result in a mixed model estimation
#' @param image_id the spe `colData` variable to mark the image
#' @param condition the spe `colData` variable to mark the condition
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param assay the assay which is used if `continuous = TRUE`
#' @param transformation the transformation to be applied as exponential e.g. 1/2 for sqrt
#' @param eps some distributional families fail if the response is zero,
#' therefore, zeros can be replaced with a very small value eps
#' @param delta the delta value to remove from the beginning of the spatial
#' statistics functions. Can be reasonable if e.g. cells are always spaced
#' by 10 µm.
#' @param family the distributional family for the functional GAM
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param ... Other parameters passed to `spatstat.explore` functions
#'
#' @returns a list with three objects: i) the dataframe with the spatial
#' statistics results, ii) the designmatrix of the inference and iii) the
#' fitted pffr object
#' @export
#'
#' @examples
#' spe <- .loadExample()
#' #make the condition a factor variable
#' colData(spe)[["patient_stage"]] <- factor(colData(spe)[["patient_stage"]])
#' #relevel to have non-diabetic as the reference category
#' colData(spe)[["patient_stage"]] <- relevel(colData(spe)[["patient_stage"]],
#' "Non-diabetic")
#' res <- spatialInference(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), correction = "rs",
#'     sample_id = "patient_id",
#'     image_id = "image_number", condition = "patient_stage",
#'     ncores = 1
#' )
spatialInference <- function(spe,
                             selection,
                             subsetby,
                             fun,
                             marks = NULL,
                             rSeq = NULL,
                             correction,
                             sample_id,
                             image_id,
                             condition,
                             continuous = FALSE,
                             assay = "exprs",
                             transformation = NULL,
                             eps = NULL,
                             delta = 0,
                             family = stats::gaussian(link = "log"),
                             ncores = 1,
                             ...){
  #small assertion that the condition has to be a factor
  stopifnot(is(colData(spe)[[condition]], "factor"))

  #first, run calcMetricPerFov
  metricRes <- calcMetricPerFov(spe = spe,
                                selection = selection,
                                subsetby = subsetby,
                                fun = fun,
                                marks =marks,
                                rSeq = rSeq,
                                by = c(sample_id, image_id, condition),
                                ncores = ncores
  )

  #second, build the dataframes for pffr and designmatrix
  #the model definitions etc should come from calcMetricPerFov in principle and
  #one of those has to be a factor with correct levels

  metricRes$ID <- paste0(
    metricRes[[condition]], "|", metricRes[[sample_id]],
    "|", metricRes[[image_id]]
  )

  # #removing field of views that have as a curve only zeros - these are cases where
  # #there is no cells of one type
  metricRes <- metricRes %>% dplyr::group_by(ID) %>%
    dplyr::filter(sum(.data[[correction]]) >= 1)
  # if a transformation should be applied to the output
  if(!is.null(transformation)){
    stopifnot(is(transformation, "numeric"))
    metricRes[[correction]] <- pmax((metricRes[[correction]])^(transformation),
                                    eps)
  }
  # else just set zeros to eps if lower than eps.
  else if(!is.null(eps)){
    metricRes[[correction]] <- pmax(metricRes[[correction]], eps)
  }

  metricRes <- metricRes %>% filter(r >= delta)

  # prepare data for FDA
  dat <- prepData(metricRes, "r", correction, sample_id,
                  image_id, condition)

  # drop rows with NA
  dat <- dat |> drop_na()
  #create the designmatrix - condition needs to be a factor with the correct
  #level at position one for the reference category
  condition <- dat[[condition]]
  print(paste0("Creating design matrix with ", levels(condition)[[1]],
               " as reference"))
  mm <- stats::model.matrix(~condition)
  #make sure that the colnames don't have "-" instead of "_"
  colnames(mm) <- gsub("-","_", colnames(mm))

  #create a formula without the first intercept column
  formula <- stats::as.formula(paste("Y ~", paste(colnames(mm)[c(-1)],
                                                  collapse="+")))

  if(!is.null(sample_id)){
    formula <- stats::as.formula(paste("Y ~",
                                paste(c(colnames(mm)[c(-1)],
                                        paste0("s(",sample_id,", bs = 're')")),
                                      collapse="+")))
  }

  #due to the removal of delta, rSeq can be less as well
  r <- metricRes$r |> unique()

  #third, run functionalGam
  mdl <- functionalGam(
    data = dat, x = r,
    designmat = mm, weights = dat$npoints,
    formula = formula,
    family = family,
    ...
  )

  #return pffr object and calcMetricPerFov dataframe in a named list
  return(list(metricRes = metricRes, designmat = mm, mdl = mdl))
}
