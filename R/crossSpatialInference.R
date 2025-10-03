#' Function for Cross Spatial Inference
#'
#' This function is a wrappere function around `spatialInference`. It calculates
#' `spatialInference` results either for all cell types in `marks`
#' (if `selection == NULL`) or for a custom subset defined in `selection`.
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
#' by 10 µm. If set to "minNnDist" it will take the mean of the minimum nearest
#' neighbour distance across all images for this cell type pair.
#' @param family the distributional family for the functional GAM
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param ... Other parameters passed to `spatstat.explore` functions
#'
#' @returns a list of objects created by the function `spatialInference`
#' with three objects: i) the dataframe with the spatial
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
#'
#' selection <- c("acinar", "ductal")
#' resLs <- crossSpatialInference(spe, selection,
#'                      subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'                       rSeq = seq(0, 50, length.out = 50), correction = "rs",
#'                       sample_id = "patient_id",
#'                       image_id = "image_number", condition = "patient_stage",
#'                       algorithm = "bam",
#'                       ncores = 1
#'                   )
#'
crossSpatialInference <- function(spe,
                                  selection = NULL,
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
  #first, create a list of all the celltypes if selection = NULL
  if(is.null(selection)){
    selection <- colData(spe)[[marks]] %>% unique()
  }

  ### code adapted from calcCrossMetricPerFov written by Samuel Gunz
  ls <- apply(base::expand.grid(selection, selection), 1, function(x) {
    return(c(x[1], x[2]))
  }, simplify = FALSE)
  # calculate the metric per FOV
  resLs <- lapply(ls, function(x) {
    res <- spatialInference(spe = spe,
                           selection = x,
                           subsetby = subsetby,
                           fun = fun,
                           marks = marks,
                           rSeq = rSeq,
                           correction = correction,
                           sample_id = sample_id,
                           image_id = image_id,
                           condition = condition,
                           continuous = continuous,
                           assay = assay,
                           transformation = transformation,
                           eps = eps,
                           delta = delta,
                           family = family,
                           ncores = ncores,
                           ...)
    return(res)
  })
  mat <- do.call("cbind",ls) %>% t()
  cellTypes <- paste0(mat[,1], "_", mat[,2])
  names(resLs) <- cellTypes

  return(resLs)
}
