# load the pancreas dataset
library("tidyr")
library("stringr")
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
  metricRes$patient_stage, "x", metricRes$patient_id,
  "x", metricRes$image_number
)
# prepare data for FDA
dat <- prepData(metricRes, "r", "rs")

# drop rows with NA
dat <- dat |> drop_na()

# create meta info of the IDs
splitData <- str_split(dat$ID, "x")
dat$condition <- factor(sapply(splitData, function(x) x[1]))
dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
dat$image_id <- factor(sapply(splitData, function(x) x[3]))
# create a designmatrix
condition <- dat$condition
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
                      conditionOnset + s(patient_id, bs = "re"))
)

test_that("Output is of correct type", {
  expect_equal(is(mdl), "pffr")
  expect_true(!is.null(mdl$coefficients))
  expect_true(!is.null(mdl$residuals))
  expect_true(!is.null(mdl$weights))
  expect_true(!is.null(mdl$fitted.values))
})

test_that("Fails if designmat and formula arguments don't correspond", {
  expect_error(functionalGam(
    data = dat, x = metricRes$r |> unique(),
    designmat = designmat, weights = dat$npoints,
    formula = formula(Y ~ conditionLong_duration_diabetes +
                        conditionOnset_diabetes + s(patient_id, bs = "re"))
  ))
})

test_that("Can handle missingnis in response - still pffr object", {
  dat[1, 2][9] <- NA
  mdl <- functionalGam(
    data = dat, x = metricRes$r |> unique(),
    designmat = designmat, weights = dat$npoints,
    formula = formula(Y ~ conditionLong_duration +
                        conditionOnset + s(patient_id, bs = "re"))
  )
  expect_equal(is(mdl), "pffr")
  expect_true(!is.null(mdl$coefficients))
  expect_true(!is.null(mdl$residuals))
  expect_true(!is.null(mdl$weights))
  expect_true(!is.null(mdl$fitted.values))
})

test_that("Should fail if weights contain NA", {
  dat$npoints[1] <- NA
  expect_error(functionalGam(
    data = dat, x = metricRes$r |> unique(),
    designmat = designmat, weights = dat$npoints,
    formula = formula(Y ~ conditionLong_duration +
                        conditionOnset + s(patient_id, bs = "re"))
  ))
})
