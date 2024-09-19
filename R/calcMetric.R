#' Compute a spatial metric on a spatial experiment object
#'
#' @param df A dataframe with the x and y coordinates from the corresponding
#' SpatialExperiment and the ColData
#' @param selection the mark(s) you want to compare
#' @param fun the spatstat function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param r_seq the range of r values to compute the function over
#' @param by the spe colData variable(s) to add to the meta data
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param window a observation window for the point pattern of class `owin`.
#' @param ... Other parameters passed to spatstat.explore functions
#'
#' @return a spatstat metric object with the fov number, the number of
#' points and the centroid of the image
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe_sub <- subset(spe, , image_number == "138")
#' df_sub <- .speToDf(spe_sub)
#' metric_res <- extractMetric(df_sub, c("alpha", "beta"),
#'     fun = "Gcross",
#'     marks = "cell_type", r_seq = seq(0, 1000, length.out = 100),
#'     by = c("patient_stage", "patient_id", "image_number")
#' )
#' @import spatstat.explore
extractMetric <- function(df,
    selection,
    fun,
    marks = NULL,
    r_seq = NULL,
    by = NULL,
    continuous = FALSE,
    window = NULL,
    ...) {
    pp <- .dfToppp(df, marks = marks, continuous = continuous, window = window)
    if (!continuous) {
        pp_sub <- subset(pp, marks %in% selection, drop = TRUE)
        meta_data <- df[, by] %>% unique()
    } else{
        pp_sub <- pp
        meta_data <- df[, by] %>% unique()
        meta_data$gene <- names(df)[names(df) %in% marks]
    }
    # small quality control to only consider pp that have more than 2 points per
    # fov and more than one unique mark and that each mark has more than one point
    if (spatstat.geom::npoints(pp_sub) > 2 &&
        ((length(unique(
            spatstat.geom::marks(pp_sub)
        )) > 1 &&
            sum(table(pp_sub$marks) > 0) > 1) ||
            length(selection) == 1)) {
        # TODO: Here I just fix the r values in the range between 0 and 500 to have
        # the same values to compare against in the library fda - that is not ideal
        metric_res <- tryCatch(
            {
                metric_res <- do.call(fun, args = list(X = pp_sub, r = r_seq, ...))
            },
            warning = function(w) {
                print(w)
                metric_res <- do.call(fun, args = list(X = pp_sub, r = r_seq, ...))
            },
            error = function(e) {
                print(e)
                metric_res <- data.frame(
                    r = r_seq,
                    fun = fun,
                    row.names = seq(1:length(r_seq))
                )
            }
        )
    }
    # This handles the case when we do cross functions for the same type
    else if (spatstat.geom::npoints(pp_sub) > 2 &&
        length(unique(selection)) == 1 &&
        length(selection) > 1) {
        metric_res <- tryCatch(
            {
                # here we use pp, otherwise there are problems with the mark connection function
                metric_res <- do.call(fun, args = list(
                    X = pp,
                    i = selection[1],
                    j = selection[2],
                    r = r_seq,
                    ...
                ))
            },
            warning = function(w) {
                print(w)
                metric_res <- do.call(fun, args = list(
                    X = pp,
                    i = selection[1],
                    j = selection[2],
                    r = r_seq,
                    ...
                ))
            },
            error = function(e) {
                print(e)
                metric_res <- data.frame(
                    r = r_seq,
                    fun = fun,
                    row.names = seq(1:length(r_seq))
                )
            }
        )
    } else {
        metric_res <- data.frame(
            r = r_seq,
            fun = fun,
            row.names = seq(1:length(r_seq))
        )
    }
    # is this needed?
    #metric_res$image_id <- df$image_number %>% unique()
    metric_res <- cbind(metric_res, meta_data)
    metric_res$npoints <- spatstat.geom::npoints(pp_sub)
    centroid <- spatstat.geom::centroid.owin(pp_sub$window)
    metric_res$centroidx <- centroid$x
    metric_res$centroidy <- centroid$y
    return(metric_res)
}

#' Calculate a spatial metric on a spatial experiment object per field of view
#'
#' @param spe a spatial experiment object
#' @param selection the mark(s) you want to compare
#' @param subsetby the spe colData variable to subset the data by
#' @param fun the spatstat function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param r_seq the range of r values to compute the function over
#' @param by the spe colData variable(s) to add to the meta data
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param ... Other parameters passed to spatstat.explore functions
#'
#' @return a dataframe of the spatstat metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metric_res <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     r_seq = seq(0, 50, length.out = 50), by = c("patient_stage", "patient_id",
#'     "image_number"),
#'     ncores = 1
#' )
#' @import dplyr parallel
calcMetricPerFov <- function(spe, selection, subsetby = NULL, fun, marks = NULL,
    r_seq = NULL, by = NULL, continuous = FALSE, ncores = 1, ...) {
    df <- .speToDf(spe)
    # we have one case for discrete cell types where we have one column to subset
    if (length(subsetby) == 1) {
        df_ls <- split(df, df[[subsetby]])
    } else {
        df_ls <- purrr::map(subsetby, ~ df %>% select(all_of(setdiff(names(df), subsetby)), .x))
    }
    metric_df <- parallel::mclapply(df_ls, function(df_sub) {
        metric_res <- extractMetric(
            df = df_sub,
            selection = selection,
            fun = fun,
            marks = marks,
            r_seq = r_seq,
            by = by,
            continuous = continuous,
            ...
        ) %>% as.data.frame()
        return(metric_res)
    }, mc.cores = ncores) %>% bind_rows()
    # store metadata of the calculation in the dataframe
    metric_df$ID <- paste0(metric_df[[by[1]]], "|", metric_df[[by[2]]])
    metric_df$fun <- fun
    metric_df$selection <- paste(selection, collapse = " and ")
    return(metric_df)
}


#' Calculate cross spatial metrics for all combinations on a spatial experiment object per field of view
#'
#' @param spe a spatial experiment object
#' @param selection the mark(s) you want to compare
#' @param subsetby the spe colData variable to subset the data by
#' @param fun the spatstat function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param r_seq the range of r values to compute the function over
#' @param by the spe colData variable(s) to add to the meta data
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param ... Other parameters passed to spatstat.explore functions
#'
#' @return a dataframe of the spatstat metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metric_res <- calcCrossMetricPerFov(spe, c("alpha", "beta", "delta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     r_seq = seq(0, 50, length.out = 50), by = c("patient_stage", "patient_id",
#'     "image_number"),
#'     ncores = 1
#' )
calcCrossMetricPerFov <- function(spe, selection, subsetby = NULL, fun,
    marks = NULL, r_seq = NULL, by = NULL,
    ncores = 1, continuous = FALSE, ...) {
    # Special case of dot functions
    if (grepl("dot", fun)) {
        # one vs all other
        ls <- unique(selection)
        # calculate the metric per FOV
        res_ls <- lapply(ls, function(x) {
            print(x)
            calcMetricPerFov(
                spe = spe, selection = x, subsetby = subsetby, fun = fun,
                marks = marks, r_seq = r_seq, by = by, ncores = ncores,
                continuous = continuous, ...
            )
        })
        # Bind the data and return
        return(bind_rows(res_ls))
    } else {
        # This creates a grid with all possible 2 way combinations
        ls <- apply(expand.grid(selection, selection), 1, function(x) {
            return(c(x[1], x[2]))
        }) |> t()

        # calculate the metric per FOV
        res_ls <- apply(ls, 1, function(x) {
            calcMetricPerFov(
                spe = spe, selection = x, subsetby = subsetby, fun = fun,
                marks = marks, r_seq = r_seq, by = by, ncores = ncores,
                continuous = continuous
            )
        })

        # Bind the data and return
        return(bind_rows(res_ls))
    }
}
