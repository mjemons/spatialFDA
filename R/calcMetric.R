#' Compute a spatial metric on a SpatialExperiment object
#'
#' A function that takes a `SpatialExperiment` object and computes a spatial
#' statistics function as implemented in `spatstat`. The output is a `spatstat`
#' object.
#'
#' @param df A `dataframe` with the x and y coordinates from the corresponding
#' `SpatialExperiment` and the `colData`
#' @param selection the mark(s) you want to compare
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param by the spe `colData` variable(s) to add to the meta data
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param window a observation window for the point pattern of class `owin`.
#' @param ... Other parameters passed to `spatstat.explore` functions
#'
#' @return a `spatstat` metric object with the fov number, the number of
#' points and the centroid of the image
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' speSub <- subset(spe, , image_number == "138")
#' dfSub <- .speToDf(speSub)
#' metricRes <- .extractMetric(dfSub, c("alpha", "Tc"),
#'     fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
#'     by = c("patient_stage", "patient_id", "image_number")
#' )
#' @import spatstat.explore
#' @importFrom methods is
.extractMetric <- function(df,
    selection,
    fun,
    marks = NULL,
    rSeq = NULL,
    by = NULL,
    continuous = FALSE,
    window = NULL,
    ...) {
    # type checking
    stopifnot(is(df, "data.frame"))
    pp <- .dfToppp(df, marks = marks, continuous = continuous, window = window)
    stopifnot("Window size must be greater than the maximum radius length considered"
              =spatstat.geom::boundingradius(pp) >= max(rSeq))
    if (!continuous) {
        ppSub <- pp[pp$marks %in% selection, drop = TRUE]
        spatstat.geom::marks(ppSub) <- factor(spatstat.geom::marks(ppSub),
                                              levels = unique(selection))
        metaData <- df[, by] %>% base::unique() %>% as.data.frame()
        base::colnames(metaData) <- by
    } else {
        ppSub <- pp
        metaData <- df[, by] %>% base::unique() %>% as.data.frame()
        base::colnames(metaData) <- by
        metaData$gene <- base::names(df)[base::names(df) %in% marks]
    }
    # small quality control to only consider pp that have more than 2 points per
    # fov and more than one unique mark and that each mark has more than one point
    if (spatstat.geom::npoints(ppSub) > 2 &&
        ((length(unique(
            spatstat.geom::marks(ppSub)
        )) > 1 &&
            #small QC to only consider points with more than 1 points in general
            base::sum(base::table(ppSub$marks) > 1) > 1) ||
            base::length(selection) == 1)) {
        metricRes <- tryCatch(
            {
                metricRes <- do.call(fun,
                    args = list(X = ppSub, r = rSeq, ...)
                )
            },
            warning = function(w) {
                message(w)
                metricRes <- do.call(fun,
                    args = list(X = ppSub, r = rSeq, ...)
                )
            },
            error = function(e) {
                message(e)
                metricRes <- data.frame(
                  r = rSeq,
                  fun = fun,
                  row.names = seq_along(rSeq)
                )
            }
        )
    # This handles the case when we do cross functions for the same type
    } else if (spatstat.geom::npoints(ppSub) > 2 &&
        base::length(base::unique(selection)) == 1 &&
        base::length(selection) > 1) {
        metricRes <- tryCatch(
            {
                # here we use pp, otherwise there are problems with the
                # mark connection function
                metricRes <- do.call(fun, args = list(
                    X = pp,
                    i = selection[1],
                    j = selection[2],
                    r = rSeq,
                    ...
                ))
            },
            warning = function(w) {
                message(w)
                metricRes <- do.call(fun, args = list(
                    X = pp,
                    i = selection[1],
                    j = selection[2],
                    r = rSeq,
                    ...
                ))
            },
            error = function(e) {
                message(e)
                metricRes <- data.frame(
                  r = rSeq,
                  fun = fun,
                  row.names = seq_along(rSeq)
                )
            }
        )
    } else {
      # TODO: the row.names are off and do funny things - need to fix this still
        metricRes <- data.frame(
            r = rSeq,
            fun = fun,
            row.names = seq_along(rSeq)
        )
    }
    metricRes <- cbind(metricRes, metaData)

    #add the total number of points as weights
    metricRes$npoints <- spatstat.geom::npoints(ppSub)
    #has to be more than 1 point in order to find a nn
    if(unique(metricRes$npoints)>1){
      #add the minimum distance from one point to the nearest neighbour
      metricRes$minDist <- min(spatstat.geom::nndist(ppSub, k = 1))
    }
    #add the individual number of points of the marked pattern if multitype
    if(spatstat.geom::is.multitype(ppSub)){
      #first, split the point pattern
      ppSubSplit <- split(ppSub)
      if(length(ppSubSplit)>1){
        npoints1 <- spatstat.geom::npoints(ppSubSplit[[1]])
        npoints2 <- spatstat.geom::npoints(ppSubSplit[[2]])
        metricRes$npointsmin <- base::pmin(npoints1, npoints2)
        metricRes$npointsmax <- base::pmax(npoints1, npoints2)
      }#else if the pattern is multitype but only one point pattern is in it
      #min is the same as total and max
      else{
        metricRes$npointsmin <- metricRes$npoints
        metricRes$npointsmax <- metricRes$npoints
      }
    }# else if the pattern is not multitype, min is the same thing as max
    # and total
    else{
      metricRes$npointsmin <- metricRes$npoints
      metricRes$npointsmax <- metricRes$npoints
    }
    centroid <- spatstat.geom::centroid.owin(ppSub$window)
    metricRes$centroidx <- centroid$x
    metricRes$centroidy <- centroid$y
    metricRes$minIntensity <- base::min(spatstat.geom::intensity(ppSub))
    if(continuous == FALSE){
        metricRes$intensityQuery <- 
        spatstat.geom::intensity(ppSub)[[as.character(selection[length(selection)])]]
    }
    metricRes$pplevels <- paste(levels(spatstat.geom::marks(ppSub)),
                                collapse = " to ")
    # small assertion that the order of the levels in `ppSub`
    # correspond to the order indicated in `selection`
    if(!continuous){
      stopifnot(levels(spatstat.geom::marks(ppSub)) == selection)
    }
    rownames(metricRes) <- NULL
    return(metricRes)
}

#' Calculate a spatial metric on a `SpatialExperiment` object per field of view
#'
#' A function that takes a `SpatialExperiment` object as input and calculates a
#' spatial metric as implemented by `spatstat` per field of view.
#'
#' @param spe a `SpatialExperiment` object
#' @param selection the mark(s) you want to compare. NOTE: This is directional.
#' c(A,B) is not the same result as c(B,A).
#' @param subsetby the spe `colData` variable to subset the data by. This
#' variable has to be provided, even if there is only one sample.
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param by the spe `colData` variable(s) to add to the meta data
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param assay the assay which is used if `continuous = TRUE`
#' @param verbose logical indicating whether to print all information or not
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param ... Other parameters passed to `spatstat.explore` functions
#'
#' @return a `dataframe` of the `spatstat` metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' metricRes <- calcMetricPerFov(spe, c("alpha", "Tc"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), by = c(
#'         "patient_stage", "patient_id",
#'         "image_number"
#'     ),
#'     ncores = 1
#' )
#' @import dplyr parallel SpatialExperiment
#' @importFrom methods is
calcMetricPerFov <- function(spe, selection, subsetby, fun, marks = NULL,
    rSeq = NULL, by = NULL, continuous = FALSE, assay = "exprs", ncores = 1,
    verbose = TRUE, ...) {
    # type checking of input
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(is(fun, "character"))
    stopifnot(is(marks, "character"))
    stopifnot(is(ncores, "numeric"))

    # check if the provide marks are in the column marks of spe colData
    if (!continuous && base::sum(!(selection %in% colData(spe)[[marks]])) > 0) {
      stop(paste0("not all marks of ", selection,
                  " are in the colData ", marks,  " of the spe"))
    }
    if(continuous) {
      expr <- SummarizedExperiment::assay(spe, assay)[marks, , drop=FALSE] %>%
        as.matrix() %>%
        t() %>%
        data.frame()
      colnames(expr) <- marks

      colData(spe) <- colData(spe) %>% cbind(expr)
    }
    df <- .speToDf(spe)
    if(length(selection)>1){
      # printing the combination calculated
        if(verbose){
            message(paste0("Calculating ", fun, " from ",
                            selection[1], " to ",
                            selection[2]))
        }
    }
    else{
      # printing the combination calculated
        if(verbose){
            message(paste0("Calculating ", fun, " of ", selection[1]))
        }
    }
    # we have one case for discrete cell types where we have one column to subset
    if (length(subsetby) == 1) {
        dfLs <- base::split(df, df[[subsetby]])
    } else {
        dfLs <- purrr::map(subsetby, ~ df %>%
            dplyr::select(dplyr::all_of(
              dplyr::setdiff(base::names(df), subsetby)), .x))
    }
    metricDf <- parallel::mclapply(dfLs, function(dfSub) {
        metricRes <- .extractMetric(
            df = dfSub,
            selection = selection,
            fun = fun,
            marks = marks,
            rSeq = rSeq,
            by = by,
            continuous = continuous,
            ...
        ) %>% as.data.frame()
        return(metricRes)
    }, mc.cores = ncores) %>% dplyr::bind_rows()
    # store metadata of the calculation in the dataframe
    metricDf$fun <- fun
    metricDf$selection <- paste(selection, collapse = " to ")
    return(metricDf)
}


#' Calculate cross spatial metrics for all combinations per FOV
#'
#' A function that takes a `SpatialExperiment` object as input and calculates a
#' cross spatial metric as implemented by `spatstat` per field of view for all
#' combinations provided by the user.
#'
#' @param spe a `SpatialExperiment` object
#' @param selection the mark(s) you want to compare
#' @param subsetby the spe `colData` variable to subset the data by
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param by the spe `colData` variable(s) to add to the meta data
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param assay the assay which is used if `continuous = TRUE`
#' @param ... Other parameters passed to spatstat.explore functions
#'
#' @return a dataframe of the `spatstat` metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' metricRes <- calcCrossMetricPerFov(spe, c("alpha", "Tc"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), by = c(
#'         "patient_stage", "patient_id",
#'         "image_number"
#'     ),
#'     ncores = 1
#' )
#'
#' @importFrom methods is
calcCrossMetricPerFov <- function(
        spe, selection, subsetby = NULL, fun,
        marks = NULL, rSeq = NULL, by = NULL,
        ncores = 1, continuous = FALSE, assay = "exprs", ...) {
    # type checking of input
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(is(fun, "character"))
    stopifnot(is(marks, "character"))
    stopifnot(is(ncores, "numeric"))

    # Special case of dot functions
    if (grepl("dot", fun)) {
        # one vs all other
        ls <- unique(selection)
        # calculate the metric per FOV
        resLs <- lapply(ls, function(x) {
            message(x)
            calcMetricPerFov(
                spe = spe, selection = x, subsetby = subsetby, fun = fun,
                marks = marks, rSeq = rSeq, by = by, ncores = ncores,
                continuous = continuous, assay, ...
            )
        })
        # Bind the data and return
        return(dplyr::bind_rows(resLs))
    } else {
        # This creates a grid with all possible 2 way combinations
        ls <- apply(base::expand.grid(selection, selection), 1, function(x) {
            return(c(x[1], x[2]))
        }) %>% t()

        # calculate the metric per FOV
        resLs <- apply(ls, 1, function(x) {
            calcMetricPerFov(
                spe = spe, selection = x, subsetby = subsetby, fun = fun,
                marks = marks, rSeq = rSeq, by = by, ncores = ncores,
                continuous = continuous, ...
            )
        })

        # Bind the data and return
        return(dplyr::bind_rows(resLs))
    }
}
