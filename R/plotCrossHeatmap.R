#' Reshaping the Result of a Cross Spatial Inference to a Dataframe
#'
#' @param resLs a list with four objects: i) the dataframe with the spatial
#' statistics results transformed and filtered as used for fitting,
#' ii) the raw spatial statistics results, iii) the designmatrix of the
#' inference and iv) the fitted pffr object v) the residual standard error per
#' condition defined as the residual sum of squares divided by the number of
#' datapoints - sum of the estimated degrees of freedom for the model parameters
#'  as well as other QC metrics
#' @param QCMetric the metric to relate the quality of the fit too.
#' @param QCThreshold the threshold on the QC metric. Depends
#' on the function used.
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
#'                       fun = "Gcross", marks = "cell_type",
#'                       rSeq = seq(0, 50, length.out = 50), correction = "rs",
#'                       sample_id = "patient_id",
#'                       image_id = "image_number", condition = "patient_stage",
#'                       algorithm = "bam",
#'                       ncores = 1
#'                   )
#' df <- extractCrossInferenceData(resLs)
extractCrossInferenceData <- function(resLs,
                                      QCMetric = "medianMinIntensity",
                                      QCThreshold = 0.1){
  df <- lapply(names(resLs), function(x){
    mdl <- resLs[[x]]$mdl
    QCDf <- resLs[[x]]$curveFittingQC
    metricRes <- resLs[[x]]$metricRes
    metricResRaw <- resLs[[x]]$metricResRaw
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
      table <- table %>% dplyr::left_join(df, by = "coefficient")  %>%
        dplyr::left_join(QCDf, by = "coefficient")
      #binarise the QC metric
      table <- table %>%
        mutate(QCBin =
                 if_else(.data[[QCMetric]] <= QCThreshold,
                         paste0(QCMetric, " <= ", QCThreshold),
                         paste0(QCMetric, " > ", QCThreshold))
        )
      return(table)
    }
  }) %>% dplyr::bind_rows()
  return(df)
}

#' Plotting the Result of a Cross Spatial Inference
#'
#' @param resLs a list with four objects: i) the dataframe with the spatial
#' statistics results transformed and filtered as used for fitting,
#' ii) the raw spatial statistics results, iii) the designmatrix of the
#' inference and iv) the fitted pffr object v) the residual standard error per
#' condition defined as the residual sum of squares divided by the number of
#' datapoints - sum of the estimated degrees of freedom for the model parameters
#'  as well as other QC metrics
#' @param adj.pvalue a pvalue adjustment method as passed to stats::p.adjust
#' defaults to Benjamini-Hochberg correction of the false discovery rate.
#' @param coefficientsToPlot list of which coefficients to plot in the heatmap
#' defaults to NULL in which case all coefficients are plotted
#' @param QCMetric the metric to relate the quality of the fit too.
#' @param QCThreshold the threshold on the Quality control metric. Depends
#' on the function used.
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
#'                       fun = "Gcross", marks = "cell_type",
#'                       rSeq = seq(0, 50, length.out = 50), correction = "rs",
#'                       sample_id = "patient_id",
#'                       image_id = "image_number", condition = "patient_stage",
#'                       algorithm = "bam",
#'                       ncores = 1
#'                   )
#' p <- plotCrossHeatmap(resLs, adj.pvalue = "BH")
#'
plotCrossHeatmap <- function(resLs,
                             adj.pvalue = "BH",
                             coefficientsToPlot = NULL,
                             QCThreshold = 1e-5,
                             QCMetric = "medianMinIntensity",
                             ...){
 df <- extractCrossInferenceData(resLs = resLs,
                                 QCMetric = QCMetric,
                                 QCThreshold = QCThreshold)
 #make below QC bin NA so that we can colour it as black later on
 df[["mean_coefficient"]] = ifelse(df[["QCBin"]] ==
                                     paste0(QCMetric," > ", QCThreshold),
                                   df[["mean_coefficient"]], NA_real_)
 if(!is.null(coefficientsToPlot)){
   df <- df %>% filter(.data[["coefficient"]] %in% coefficientsToPlot)
 }
 if(is.null(adj.pvalue)){
   p <- ggplot(df, aes(x = .data[["cell1"]], y = .data[["cell2"]])) +
     geom_point(aes(size = -log10(.data[["p-value"]] + 0.001),
                    shape = .data[["QCBin"]],
                    color = .data[["mean_coefficient"]]))+
     geom_point(aes(size = -log10(.data[["p-value"]] + 0.001)),
                    shape = 1, colour = "black")
 }else{
   df[["adj.p-value"]] <- stats::p.adjust(df[["p-value"]], method = adj.pvalue)
   p <- ggplot(df, aes(x = .data[["cell1"]], y = .data[["cell2"]])) +
     geom_point(aes(size = -log10(.data[["adj.p-value"]] + 0.001),
                    shape = .data[["QCBin"]],
                    color = .data[["mean_coefficient"]]))+
     geom_point(aes(size = -log10(.data[["adj.p-value"]] + 0.001)),
                shape = 1,colour = "black")
 }
 p <- p + scale_x_discrete(guide = guide_axis(angle = 50))
 if(any(df[["QCBin"]] == paste0(QCMetric," <= ", QCThreshold), na.rm = TRUE)){
   p <- p + scale_shape_manual(values = c(13,16))
 }else{
   p <- p + scale_shape_manual(values = c(16,13))
 }

 p <- p + facet_wrap(~.data[["coefficient"]]) +
   theme_light() +
   scale_colour_gradient2(midpoint = 0,
                          high = scales::muted("red"),
                          mid = "white",
                          low = scales::muted("blue"),
                          na.value = "black")
 return(p)
}
