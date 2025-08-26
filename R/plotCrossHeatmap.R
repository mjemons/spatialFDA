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
    if(!is.null(mdl)){
      table <- summary(mdl)$s.table %>% as.data.frame()
      table$coefficient <- rownames(table)
      table$combination <- x
      table <- table %>% separate(.data[["combination"]],
                                  c("cell1", "cell2"), sep = "_")

      #extract the mean functional coefficient as effect size measure
      coef <- coef(mdl)
      df <- lapply(rownames(table), function(predictor){
        df <- coef$sm[[predictor]]$coef
        df[["coefficient"]] <- predictor
        df <- df %>% mutate(mean_coefficient = mean(.data[["value"]])) %>%
          select(.data[["mean_coefficient"]], .data[["coefficient"]]) %>%
          unique()
        return(df)
      }) %>% dplyr::bind_rows()
      #join by coefficient
      table <- table %>% dplyr::left_join(df, by = "coefficient")
      return(table)
    }
  }) %>% dplyr::bind_rows()
  return(df)
}

#' Plotting the Result of a Cross Spatial Inference
#'
#' @param resLs a list of objects created by the function `spatialInference`
#' with three objects: i) the dataframe with the spatial
#' statistics results, ii) the designmatrix of the inference and iii) the
#' fitted pffr object
#' @param adj.pvalue a pvalue adjustment method as passed to stats::p.adjust
#' defaults to Benjamini-Hochberg correction of the false discovery rate.
#' @param coefficientsToPlot list of which coefficients to plot in the heatmap
#' defaults to NULL in which case all coefficients are plotted
#' @param ... other parameters passed to `ggplot2` functions
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
#' p <- plotCrossHeatmap(resLs, adj.pvalue = "BH")
#'
plotCrossHeatmap <- function(resLs,
                             adj.pvalue = "BH",
                             coefficientsToPlot = NULL,
                             ...){
 df <- extractCrossInferenceData(resLs)
 if(!is.null(coefficientsToPlot)){
   df <- df %>% filter(.data[["coefficient"]] %in% coefficientsToPlot)
 }
 if(is.null(adj.pvalue)){
   p <- ggplot(df, aes(x = .data[["cell1"]], y = .data[["cell2"]])) +
     geom_point(aes(color = .data[["mean_coefficient"]],
                    size = -log10(.data[["p-value"]] + 0.001)))
 }else{
   df[["adj.p-value"]] <- stats::p.adjust(df[["p-value"]], method = adj.pvalue)
   p <- ggplot(df, aes(x = .data[["cell1"]], y = .data[["cell2"]])) +
     geom_point(aes(color = .data[["mean_coefficient"]],
                    size = -log10(.data[["adj.p-value"]] + 0.001)))
 }
 p <- p + scale_x_discrete(guide = guide_axis(angle = 50)) +
   geom_point(aes(size = -log10(.data[["p-value"]] + 0.001)),
              shape = 1,colour = "black")+
   facet_wrap(~.data[["coefficient"]]) +
   theme_light() +
   scale_colour_gradient2(midpoint = 0,
                          high = scales::muted("red"),
                          mid = "white",
                          low = scales::muted("blue"))
 return(p)
}
