#' Compute a spatial metric on a spatial experiment object
#'
#' @param spe a spatial experiment object subset to a single image
#' @param selection the mark(s) you want to compare
#' @param fov the field of view number
#' @param fun the spatstat function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param r_seq the range of r values to compute the function over
#'
#' @return a spatstat metric object with the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#'  spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#'  spe_sub <- subset(spe, ,image_number == '138')
#'  metric_res <- extractMetric(spe_sub, c('alpha', 'beta'), 138, fun = 'Gcross'
#'  , marks = 'cell_type', r_seq = seq(0,1000, length.out = 100))
#' @import spatstat.explore
extractMetric <- function(spe, selection, fov, fun, marks = NULL, r_seq = NULL){
  pp <- .speToppp(spe, marks = marks)
  patient_id <- spe$patient_id %>% unique()
  condition <- spe$patient_stage %>% unique()
  pp_sub <- subset(pp, marks %in% selection, drop = TRUE)

  #small quality control to only consider pp that have more than 100 points per fov and more than one unique mark
  if(spatstat.geom::npoints(pp_sub)>100 && (length(unique(spatstat.geom::marks(pp_sub)))>1 || length(selection) == 1)){ #& npoints(pp1) > 20 & npoints(pp2) > 20){
    #TODO: Here I just fix the r values in the range between 0 and 500 to have the same values to compare against in the library fda - that is not ideal --> need to reconsider!
    metric_res <- do.call(fun, args = list(X=pp_sub, r = r_seq))
    metric_res$image_id = fov
    metric_res$patient_id = patient_id
    metric_res$condition = condition
    metric_res$npoints = spatstat.geom::npoints(pp_sub)
    centroid <- spatstat.geom::centroid.owin(pp_sub$window)
    metric_res$centroidx <- centroid$x
    metric_res$centroidy <- centroid$y
    return(metric_res)
  }
}

#' Calculate a spatial metric on a spatial experiment object per field of view
#'
#' @param spe a spatial experiment object
#' @param selection the mark(s) you want to compare
#' @param fun the spatstat function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param r_seq the range of r values to compute the function over
#' @param ncores the number of cores to use for parallel processing, default is 1
#'
#' @return a dataframe of the spatstat metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#'  spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#'  metric_res <- calcMetricPerFov(spe, c('alpha', 'beta'), fun = 'Gcross'
#'  , marks = 'cell_type', r_seq = seq(0,50, length.out = 50), ncores = 2)
#' @import dplyr parallel
calcMetricPerFov <- function(spe, selection, fun, marks = NULL, r_seq = NULL, ncores = 1){
  fovs <- spe$image_number %>% as.list() %>% unique()
  metric_df <- parallel::mclapply(fovs, function(x){
    spe_sub <- subset(spe, ,image_number == x)
    metric_res <- extractMetric(spe = spe_sub, selection = selection, fov = x, fun = fun, marks = marks, r_seq = r_seq) %>% as.data.frame()
    return(metric_res)
  }, mc.cores = ncores) %>% bind_rows()
  return(metric_df)
}


