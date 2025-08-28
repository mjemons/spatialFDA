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
  metricRes$patient_stage, "x", metricRes$patient_id,
  "x", metricRes$image_number
)
#removing field of views that have as a curve only zeros - these are cases where
#there is no cells of one type
metricRes <- metricRes %>% dplyr::group_by(ID) %>% dplyr::filter(sum(rs) >= 1)

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
mdl1 <- functionalGam(
  data = dat, x = metricRes$r |> unique(),
  designmat = designmat, weights = dat$npoints,
  formula = formula(Y ~ conditionLong_duration +
                      conditionOnset + s(patient_id, bs = "re")),
  family = gaussian(link = "log")
)

## do the same with spatialInference

#make the condition a factor variable
colData(spe)[["patient_stage"]] <- factor(colData(spe)[["patient_stage"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["patient_stage"]] <- relevel(colData(spe)[["patient_stage"]],
                                           "Non-diabetic")
res <- spatialInference(spe, c("alpha", "Tc"),
                        subsetby = "image_number", fun = "Gcross", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1
)

mdl2 <- res$mdl

test_that("wrapper function gives same result as manual computation", {
  expect_true(identical(mdl1$coefficients, mdl2$coefficients))
})

resBeta <- spatialInference(spe, "beta",
                        subsetby = "image_number", fun = "Gest", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1
)

test_that("spatialInference handels case when one condition has no images with
          calculated curves ", {
  expect_true(is.null(resBeta$mdl) && is.null(resBeta$designmat))
})

res <- spatialInference(spe, c("alpha", "Tc"),
                        subsetby = "image_number", fun = "Gcross", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        weights = "min",
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1
)

#order
metricRes <- res$metricRes %>% arrange(patient_id)

manual_weights <- metricRes$npointsmin / mean(metricRes$npointsmin)

modelframe <- res$mdl$model %>% arrange(patient_id)

test_that("weights of the model are really the min weights expected", {
  expect_equal(sort(modelframe$`(weights)`), sort(manual_weights))
})

res <- spatialInference(spe, c("alpha", "Tc"),
                        subsetby = "image_number", fun = "Gcross", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        weights = "max",
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1
)

#order
metricRes <- res$metricRes %>% arrange(patient_id)

manual_weights <- metricRes$npointsmax / mean(metricRes$npointsmax)

modelframe <- res$mdl$model %>% arrange(patient_id)

test_that("weights of the model are really the max weights expected", {
  expect_equal(sort(modelframe$`(weights)`), sort(manual_weights))
})

res <- spatialInference(spe, c("alpha", "Tc"),
                        subsetby = "image_number", fun = "Gcross", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        weights = NULL,
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1
)

#order
metricRes <- res$metricRes %>% arrange(patient_id)

manual_weights <- seq.int(from = 1, to = 1, length.out = nrow(metricRes))

modelframe <- res$mdl$model %>% arrange(patient_id)

test_that("weights of the model are really the max weights expected", {
  expect_equal(sort(modelframe$`(weights)`), sort(manual_weights))
})


