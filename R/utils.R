#' Convert SpatialExperiment object to ppp object
#'
#' @param spe A SpatialExperiment object subset to a single image
#' @param marks A vector of marks to be associated with the points
#' @return A ppp object for use with `spatstat` functions
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe_sub <- subset(spe, ,image_number == '138')
#' pp <- .speToppp(spe_sub, marks = 'cell_type')
#'
#' @importFrom SummarizedExperiment colData
.speToppp <- function(spe, marks = NULL){
  df <- data.frame(x = SpatialExperiment::spatialCoords(spe)[,1], y = SpatialExperiment::spatialCoords(spe)[,2])
  #is this definition of the window actually correct? Do I underestimate it?
  pp <- spatstat.geom::as.ppp(df, W = spatstat.geom::owin(c(min(df$x)-1,max(df$x)+1),c(min(df$y)-1,max(df$y)+1)))
  #set the marks
  spatstat.geom::marks(pp) <- factor(colData(spe)[[marks]])
  return(pp)
}

#' Apply .as_ppp to a list of images
#'
#' @param spe a spatial experiment object
#' @param subset a vector of image numbers to subset the spatial experiment object
#' @param marks a vector of marks to be associated with the points
#'
#' @return a list of ppp objects, one for each image in the subset
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' pp_ls <- .apply_ls_ppp(spe, list('138','139'),marks = 'cell_type' )
#'
.apply_ls_ppp <- function(spe, subset, marks = NULL){
  pp_ls <- lapply(subset, function(x){
    spe_sub <- subset(spe, ,image_number == x)
    pp <- .speToppp(spe_sub, marks = marks)
    return(pp)
  })
  return(pp_ls)
}
