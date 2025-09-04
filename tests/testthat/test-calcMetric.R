# retrieve example data from Damond et al. (2019)
spe <- .loadExample()
rSeq <- seq(0, 50, length.out = 50)

## test function calcMetricPerFov
test_that("Output contains correction for discrete single mark", {
  rSeq <- seq(0, 50, length.out = 50)
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "Gest",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "rs",
                                ncores = 1
  )
  expect_contains(colnames(metricRes), "rs")
  expect_contains(colnames(metricRes), "r")
  expect_contains(colnames(metricRes), "theo")
})

test_that("Output contains correction for continuous single mark", {
  # add continuous mark to colData
  protein <- "CD31"

  rSeq <- seq(0, 50, length.out = 50)
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "markcorr",
                                marks = protein,
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "iso",
                                ncores = 1,
                                continuous = TRUE
  )
  expect_contains(colnames(metricRes), "iso")
  expect_contains(colnames(metricRes), "r")
  expect_contains(colnames(metricRes), "theo")
})

test_that("Output contains correction for two marks", {
  metricRes <- calcMetricPerFov(spe, c("alpha", "Tc"),
                                subsetby = "image_number", fun = "Gcross",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "rs",
                                ncores = 1
  )
  expect_contains(colnames(metricRes), "rs")
  expect_contains(colnames(metricRes), "r")
  expect_contains(colnames(metricRes), "theo")
})

test_that("Output has correct dimensions", {
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "Gest",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "rs",
                                ncores = 1
  )
  expect_length(metricRes$rs, length(rSeq) * length(unique(spe$image_name)))
})

test_that("Function fails if marks not in ColData", {
  expect_error(calcMetricPerFov(spe, c("alpha", "epsilon"),
                                subsetby = "image_number", fun = "Gcross",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"
                                ),
                                ncores = 1
  ))
})

# Test function calcCrossMetricPerFov
test_that("Cross function output has correct dimensions", {
  selection <- c("alpha", "beta", "delta")
  metricRes <- calcCrossMetricPerFov(spe, selection,
                                     subsetby = "image_number", fun = "Gcross",
                                     marks = "cell_type",
                                     rSeq = seq(0, 50, length.out = 50), by = c(
                                       "patient_stage", "patient_id",
                                       "image_number"
                                     ),
                                     ncores = 1
  )

  expect_length(metricRes$rs,
                (length(rSeq) * length(unique(spe$image_name))
                * (length(selection)^2)))
})

test_that("Cross function output has correct dimensions for Kdot", {
  selection <- c("alpha", "beta", "delta")
  metricRes <- calcCrossMetricPerFov(spe, selection,
                                     subsetby = "image_number", fun = "Kdot",
                                     marks = "cell_type",
                                     rSeq = seq(0, 50, length.out = 50), by = c(
                                       "patient_stage", "patient_id",
                                       "image_number"
                                     ),
                                     correction = "border",
                                     ncores = 1
  )

  expect_length(metricRes$border,
                length(rSeq) * (length(unique(spe$image_name)))
                * (length(selection)))
})

test_that("Numeric results are correct for Lcross image 148", {
  selection <- c("alpha", "Tc")
  speSub <- subset(spe, , image_number == "148")
   dfSub <- .speToDf(speSub)
   metricRes <- .extractMetric(dfSub, selection,
       fun = "Lcross",
       marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
       by = c("patient_stage", "patient_id", "image_number")
   ) %>% as.data.frame()
  resObserved <- metricRes %>% dplyr::select(dplyr::all_of(c("r", "theo", "border", "trans", "iso"))) %>%
    dplyr::slice(1,4,19,38,45,50)
  resExpected <- read.csv('lcross_single.csv')

  expect_equal(object = resObserved[complete.cases(resObserved), ],
               expected = resExpected[complete.cases(resExpected), ])
})

test_that("Numeric results are correct for Lcross", {
  selection <- c("alpha", "Tc")
  metricRes <- calcMetricPerFov(spe = spe,
                                selection = selection,
                                subsetby = "image_number",
                                fun = "Lcross",
                                marks = "cell_type",
                                rSeq = seq(0, 50, length.out = 50),
                                by = c("patient_stage", "patient_id",
                                     "image_number"),
                                ncores = 1)
  resObserved <- metricRes %>% dplyr::select(dplyr::all_of(c("r", "theo", "border", "trans", "iso"))) %>%
    dplyr::slice(5,29,432,1100,3872,4983)
  resExpected <- read.csv('lcross.csv')

  expect_equal(object = resObserved[complete.cases(resObserved), ],
               expected = resExpected[complete.cases(resExpected), ])
})

test_that("Numeric results are correct for Gcross", {
  selection <- c("alpha", "Tc")
  metricRes <- calcMetricPerFov(spe = spe,
                                selection,
                                subsetby = "image_number",
                                fun = "Gcross",
                                marks = "cell_type",
                                rSeq = seq(0, 50, length.out = 50),
                                by = c("patient_stage", "patient_id",
                                       "image_number"),
                                ncores = 1)
  resObserved <- metricRes %>% dplyr::select(dplyr::all_of(c("r", "theo", "han", "rs", "km"))) %>%
    dplyr::slice(5,29,432,1100,3872,4983)
  resExpected <- read.csv('gcross.csv')

  expect_equal(object = resObserved[complete.cases(resObserved), ],
               expected = resExpected[complete.cases(resExpected), ])
})

test_that("Min Points + Max Points = Total Points", {
  selection <- c("alpha", "Tc")
  metricRes <- calcMetricPerFov(spe = spe,
                                selection,
                                subsetby = "image_number",
                                fun = "Gcross",
                                marks = "cell_type",
                                rSeq = seq(0, 50, length.out = 50),
                                by = c("patient_stage", "patient_id",
                                       "image_number"),
                                ncores = 1)
  minPoints <- metricRes$npointsmin
  maxPoints <- metricRes$npointsmax
  totalPoints <- metricRes$npoints

  expect_equal(object = minPoints + maxPoints,
               expected = totalPoints)
})

test_that("spatstat.geom::nndist has a working check if only 1 point is available",{
  selection <- "beta"
  metricRes <- calcMetricPerFov(spe = spe,
                                selection,
                                subsetby = "image_number",
                                fun = "Gest",
                                marks = "cell_type",
                                rSeq = seq(0, 50, length.out = 50),
                                by = c("patient_stage", "patient_id",
                                       "image_number"),
                                ncores = 1)

  metricResSub <- metricRes %>% subset(patient_id == 6180)

  expect_true(sum(is.na(metricResSub$minDist))==nrow(metricResSub))
})
