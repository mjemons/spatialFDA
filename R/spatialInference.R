#' Statistical Inference on Spatial Statistics Functions
#'
#' A function to perform spatial statistical inference on spatial omics data.
#' This function works so far only on functions of radius "r".
#'
#' @param spe a `SpatialExperiment` object
#' @param selection the mark(s) you want to compare. NOTE: This is directional.
#' c(A,B) is not the same result as c(B,A).
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
#' or Fisher's variance-stabilising transformation if "Fisher"
#' @param weights the weighting to be applied to the functional GAM. Either NULL
#' (equal weights), total (npoints of total pattern), min (npoints of the smaller
#' subpattern) or max (npoints of the larger subpattern) or a user defined value
#' of same length as the number of curves to be estimated
#' @param eps some distributional families fail if the response is zero,
#' therefore, zeros can be replaced with a very small value eps
#' @param delta the delta value to remove from the beginning of the spatial
#' statistics functions. Can be reasonable if e.g. cells are always spaced
#' by 10 µm. If set to "minNnDist" it will take the mean of the minimum nearest
#' neighbour distance across all images for this cell type pair.
#' @param family the distributional family for the functional GAM
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param verbose logical indicating whether to print all information or not
#' @param upperDeltaProb the quantile to filter out the constant 1 part for `Gest`
#' and `Gcross`. If `NULL` no upper filtering is applied.
#' @param weightTransform logical indicating whether the weights (number of points) 
#' should be sqrt transformed
#' @param AR1 logical indicating whether to calculate the autocorrelation of the 
#' residuals along the domain and account for this in a second fitting step
#' @param sandwich string indicating how and if to adjust for heterscedasticity of the 
#' residuals with a sandwich correction
#' @param algorithm algorithm to fit the refund::pffr method. defaults to `bam`
#' @param ... Other parameters passed to `spatstat.explore` functions for
#' parameters concerning the spatial function calculation and to `refund::pffr`
#' for the functional additive mixed model inference
#'
#' @returns a list with four objects: i) the dataframe with the spatial
#' statistics results transformed and filtered as used for fitting,
#' ii) the raw spatial statistics results, iii) the designmatrix of the
#' inference and iv) the fitted pffr object v) the residual standard error per
#' condition defined as the residual sum of squares divided by the number of
#' datapoints - sum of the estimated degrees of freedom for the model parameters
#'  as well as other QC metrics
#'
#' @export
#'
#' @examples
#' spe <- .loadExample()
#' #make the condition a factor variable
#' colData(spe)[["patient_stage"]] <- factor(colData(spe)[["patient_stage"]])
#' #relevel to have non-diabetic as the reference category
#' colData(spe)[["patient_stage"]] <- relevel(colData(spe)[["patient_stage"]],
#' "Non-diabetic")
#' res <- spatialInference(spe, c("alpha", "Tc"),
#'     fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), correction = "rs",
#'     sample_id = "patient_id",
#'     image_id = "image_number", condition = "patient_stage",
#'     ncores = 1,
#'     algorithm = "bam"
#' )
spatialInference <- function(spe,
                             selection,
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
                             weights = "total",
                             eps = 1e-3,
                             delta = "minNnDist",
                             family = stats::gaussian(link = "log"),
                             verbose = TRUE,
                             upperDeltaProb = NULL,
                             weightTransform = TRUE,
                             AR1 = FALSE,
                             sandwich = "cluster",
                             ncores = 1,
                             algorithm = "bam",
                             ...){
  #for computational reasons, remove the assays as we don't need them
  SummarizedExperiment::assays(spe) <- list()
  #for computational reasons, remove the rowData as we don't need them
  SummarizedExperiment::rowData(spe) <- S4Vectors::DataFrame(row.names = rownames(spe))
  #small assertion that the condition has to be a factor
  stopifnot(is(colData(spe)[[condition]], "factor"))

  #first, run calcMetricPerFov
  metricResRaw <- calcMetricPerFov(spe = spe,
                                selection = selection,
                                subsetby = image_id,
                                fun = fun,
                                marks =marks,
                                rSeq = rSeq,
                                by = c(sample_id, image_id, condition),
                                verbose = verbose,
                                ncores = ncores,
                                correction = correction,
                                ...
  )

  #second, build the dataframes for pffr and designmatrix
  #the model definitions etc should come from calcMetricPerFov in principle and
  #one of those has to be a factor with correct levels
  if(!is.null(sample_id)){
    metricResRaw$ID <- paste0(
      metricResRaw[[condition]], "|", metricResRaw[[sample_id]],
      "|", metricResRaw[[image_id]]
    )
  }else{
    metricResRaw$ID <- paste0(
      metricResRaw[[condition]], "|",
      "|", metricResRaw[[image_id]]
    )
  }

  noConditionsPreFiltering <- (length(unique(metricResRaw[[condition]])))
  # #removing field of views that have as a curve only zeros - these are cases where
  # #there is no cells of one type
  metricRes <- metricResRaw %>% dplyr::group_by(.data[["ID"]]) %>%
    dplyr::filter(sum(.data[[correction]]) >= 1)

  #filter the upper part of the curve 
  if(!is.null(upperDeltaProb) && (fun == "Gest" || fun == "Gcross")){
    res <-metricRes |> 
      filter(round(.data[[correction]], 2) == 1) |> 
      group_by(.data[["ID"]]) |> 
      mutate(lowerRQuartile = stats::quantile(r, probs = upperDeltaProb))
    res2 <- metricRes |>
      filter(round(.data[[correction]], 2) != 1 & r == max(rSeq)) |>
      group_by(.data[["ID"]]) |>
      mutate(lowerRQuartile =  max(rSeq))
    res <- rbind(res, res2)
    upperDelta <- stats::median(res$lowerRQuartile)
    metricRes <- metricRes %>% filter(r < upperDelta)
  }

  # if a transformation should be applied to the output
  if(!is.null(transformation)){
    if(transformation == "Fisher"){
      metricRes[[correction]] <- pmax(asin(sqrt(metricRes[[correction]])),
                                      eps)
    }else{
    stopifnot(is(transformation, "numeric"))
    metricRes[[correction]] <- pmax((metricRes[[correction]])^(transformation),
                                    eps)
    }
  }
  # else just set zeros to eps if lower than eps.
  else if(!is.null(eps)){
    metricRes[[correction]] <- pmax(metricRes[[correction]], eps)
  }

  if(delta == "minNnDist"){
    delta <- stats::weighted.mean(x=metricRes[["minDist"]],
                                  w = metricRes[["npoints"]])
  }


  metricRes <- metricRes %>% filter(r >= delta)
  noConditionsPostFiltering <- (length(unique(metricRes[[condition]])))
  if(noConditionsPreFiltering == noConditionsPostFiltering){
    # prepare data for FDA
    dat <- prepData(metricRes, "r", correction, sample_id,
                    image_id, condition)

    # drop rows with NA
    dat <- dat |> drop_na()
    #create the designmatrix - condition needs to be a factor with the correct
    #level at position one for the reference category
    conditionVariable <- condition
    condition <- dat[[condition]]
    if(verbose){
      message(paste0("Creating design matrix with ", levels(condition)[[1]],
                 " as reference"))
    }
    mm <- stats::model.matrix(~condition)
    #make sure that the colnames don't have "-" instead of "_"
    colnames(mm) <- gsub("-","_", colnames(mm))
    #create a formula without the first intercept column
    formula <- stats::as.formula(paste("Y ~", paste(colnames(mm)[c(-1)],
                                                    collapse="+")), env = emptyenv())

    if(!is.null(sample_id)){
      formula <- stats::as.formula(paste("Y ~",
                                  paste(c(colnames(mm)[c(-1)],
                                          paste0("s(",sample_id,", bs = 're')")),
                                        collapse="+")), env = emptyenv())
    }

    #due to the removal of delta, rSeq can be less as well
    r <- metricRes$r |> unique()

    # define the weights
    if(is.null(weights)){
      # give each observation weight one if no weights are passed
      weights = seq.int(from = 1, to = 1, length.out = nrow(dat))
    }else if(weights == "total"){
      weights = dat$npoints
    }else if(weights == "min"){
      weights = dat$npointsmin
    }else if(weights == "max"){
      weights = dat$npointsmax
    }else{
      stopifnot(length(weights) == nrow(dat))
      weights = weights
    }

    if(weightTransform){
      weights = sqrt(weights)
    }
    
    if(AR1){
      #if there is an AR1 correlation parameter given, fit first a model and compute the median ACF 
      #of the residuals
      mdl <- functionalGam(
        data = dat, x = r,
        designmat = mm, weights = weights,
        formula = formula,
        family = family,
        algorithm = algorithm,
        sandwich = sandwich,
        ...
      )
      if(algorithm == "gamm4"){
        mdl = mdl$gam
      }
      #compute median ACF for a lag of 1
      rho_est <- apply(stats::resid(mdl), 1, function(x) stats::acf(x, plot = FALSE)$acf[2]) |> 
        stats::median(na.rm = TRUE)
      #refit with estimated residual autocorrelation
      mdl <- functionalGam(
        data = dat, x = r,
        designmat = mm, weights = weights,
        formula = formula,
        family = family,
        rho = rho_est,
        algorithm = algorithm,
        sandwich = sandwich,
        ...
      )
    }else{
      mdl <- functionalGam(
        data = dat, x = r,
        designmat = mm, weights = weights,
        formula = formula,
        family = family,
        algorithm = algorithm,
        sandwich = sandwich,
        ...
      )
    }
    
    if(algorithm == "gamm4"){
      mdl = mdl$gam
    }
    
    ### Calculation of metrics assessing the quality of the model fit

    # adj R-squared of the entire model
    Rsq.adj <- summary(mdl, re.test = FALSE)$r.sq
    if(verbose){
      message(paste0("The adjusted R-squared of the model is ", Rsq.adj))
    }

    ##rename the conditions to be the same as in the summary output
    dat <- dat %>%
      mutate(coefficient =
               paste0("condition",
                      gsub("-","_", .data[[conditionVariable]]),"(yindex)")) %>%
      #rename the reference category to be Intercept
      mutate(coefficient =
               case_when(coefficient ==
                           paste0("condition",
                                  gsub("-","_",levels(condition)[[1]]), "(yindex)")
                         ~ "Intercept(yindex)", TRUE ~ coefficient))

    # calculate the median intensity per condition
    dfIntensity <- dat %>%
      group_by(.data[["coefficient"]]) %>%
      mutate(medianMinIntensity = stats::median(.data[["minIntensity"]])) %>%
      select(.data[["coefficient"]], .data[["medianMinIntensity"]]) %>%
      unique()

    #another QC metric of the model fit is inspecting the residuals per condition
    #we compare the residual standard error which is the sqrt residual sum of
    #squares divided by the degrees of freedom of the residuals per condition

    #we need condition specific residual degrees of freedom
    #take the rows of dat as this is filtered. Since we need to have not only
    #the number of curves but also the values per curve, we multiply nrow of dat
    #with the length of the functional domain r

    #Furthermore, we need to get condition specific estimated degrees of freedom
    #of the model parameters. These are in the model summary

    df.edf <- summary(mdl, re.test = FALSE)[["s.table"]] %>% as.data.frame()

    #if it is a mixed model, we need to remove the random effect column
    if(!is.null(sample_id)){
      #filter out sample_id
      df.edf <- df.edf %>% filter(!grepl(sample_id, rownames(df.edf)))
      #filter out image_id
      df.edf <- df.edf %>% filter(!grepl(image_id, rownames(df.edf)))
    }
    #this assumes that the order of the levels is the same as the order of the
    #summary output
    df.edf[["coefficient"]] <- rownames(df.edf)

    #select only the edf and the condition variable
    df.edf <- df.edf %>% select(.data[["edf"]], .data[["coefficient"]])

    #subtract from each condition wise number of curves * number of datapoints
    #per curve the condition wise edf from the summary(mdl) output above

    df.residual <- dat %>%
      group_by(.data[["coefficient"]]) %>%
      summarise(no.datapoints = n() * length(r)) %>%
      left_join(df.edf, by = "coefficient") %>%
      mutate(df.residual.condition = .data[["no.datapoints"]] - .data[["edf"]])

    #now we extract the residuals

    residualPffr <- as.data.frame(stats::residuals(mdl))
    residualPffr[[conditionVariable]] <- dat[[conditionVariable]]

    residualPffr <- residualPffr %>%
      mutate(coefficient = paste0("condition",
                                  gsub("-","_",
                                       .data[[conditionVariable]]),"(yindex)")) %>%
      #rename the reference category to be Intercept
      mutate(coefficient = case_when(coefficient == paste0("condition",
                                                           gsub("-","_",
                                                                levels(condition)[[1]]), "(yindex)")
                                     ~ "Intercept(yindex)",
                                     TRUE ~ coefficient))
    # combine the residuals with the degrees of freedom
    residualPffr <- residualPffr %>% left_join(df.residual, by = "coefficient")
    #calculate the grouped RSS and divide by the grouped condition wise residuals
    #and take the sqrt of this
    residualDf <- residualPffr %>%
      group_by(.data[["coefficient"]]) %>%
      reframe(residual_standard_errors = sqrt(sum(across(where(is.numeric) &
                                                           !c(.data[["df.residual.condition"]],
                                                              .data[["edf"]],
                                                              .data[["no.datapoints"]]))**2)
                                               /.data[["df.residual.condition"]]),
              edf = .data[["edf"]]) %>%
      unique()

    #assemble all QC scores
    QCDf <- residualDf %>%
      left_join(dfIntensity, by = "coefficient")

  }else{
    if(verbose){
      message("Can not fit a model if one condition has no images with curves")
    }
    mdl = NULL
    mm = NULL
    QCDf = NULL
  }

  #return pffr object and calcMetricPerFov dataframe in a named list
  return(list(metricRes = metricRes,
              designmat = mm,
              mdl = mdl,
              curveFittingQC = QCDf))
}