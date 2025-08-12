#' Reshaping the Result of a Cross Spatial Inference to a Dataframe
#'
#' @param resLs a list of objects created by the function `spatialInference`
#' with three objects: i) the dataframe with the spatial
#' statistics results, ii) the designmatrix of the inference and iii) the
#' fitted pffr object
#'
#' @returns a dataframe for plotting with ggplot2
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
#'                       ncores = 1
#'                   )
#' df <- extractCrossInferenceData(resLs)
extractCrossInferenceData <- function(resLs){
  df <- lapply(names(resLs), function(x){
    mdl <- resLs[[x]]$mdl
    table <- summary(mdl)$s.table %>% as.data.frame()
    table$condition <- rownames(table)
    table$combination <- x
    table <- table %>% separate(combination, c("cell1", "cell2"), sep = "_")
    return(table)
  }) %>% dplyr::bind_rows()
  return(df)
}

#' Plotting the Result of a Cross Spatial Inference
#'
#' @param resLs a list of objects created by the function `spatialInference`
#' with three objects: i) the dataframe with the spatial
#' statistics results, ii) the designmatrix of the inference and iii) the
#' fitted pffr object
#'
#' @returns a ggplot2 object
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
#'                       ncores = 1
#'                   )
#' p <- plotCrossHeatmap(resLs)
#'
plotCrossHeatmap <- function(resLs){
 df <- extractCrossInferenceData(resLs)
 p <- ggplot(df, aes(.data[["cell1"]], .data[["cell2"]], fill = -log10(.data[["p-value"]]))) +
   geom_tile(colour="white", size=0.2) +
   facet_wrap(~.data[["condition"]]) +
   theme_light()
 return(p)
}
