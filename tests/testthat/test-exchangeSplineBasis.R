test_that("exchangeSplineBasis contains new bs in first place", {
  # load the pancreas dataset
  library("tidyr")
  library("dplyr")
  # retrieve example data from Damond et al. (2019)
  spe <- .loadExample()
  # calculate the Gcross metric for alpha and Tc cells
  metricRes <- calcMetricPerFov(spe, c("alpha", "Tc"),
      subsetby = "image_number", fun = "Gcross",
      marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
      c("patient_stage", "patient_id", "image_number"), ncores = 1
  )
  metricRes$ID <- paste0(
    metricRes$patient_stage, "|", metricRes$patient_id,
    "|", metricRes$image_number
  )
  dat <- prepData(metricRes, "r", "rs", sample_id = "patient_id",
      image_id = "image_number", condition = "patient_stage")

  #' # drop rows with NA
  dat <- dat |> drop_na()

  # create a designmatrix
  condition <- dat$patient_stage
  # relevel the condition - can set explicit contrasts here
  condition <- relevel(condition, "Non-diabetic")
  designmat <- model.matrix(~condition)
  # colnames don't work with the '-' sign
  colnames(designmat) <- c(
      "(Intercept)", "conditionLong_duration",
      "conditionOnset"
  )
  # fit the model
  mdl <- functionalGam(
      data = dat, x = metricRes$r |> unique(),
      designmat = designmat, weights = dat$npoints,
      formula = formula(Y ~ conditionLong_duration +
          conditionOnset + s(patient_id, bs = "re")),
          fit = FALSE
  )
  newFormula <- exchangeSplineBasis(mdl = mdl, var = "x.vec", bs = "mpi")
  expect_true(grepl('bs\\s*=\\s*"mpi"', formula.tools::rhs.vars(newFormula)[1]))
})
