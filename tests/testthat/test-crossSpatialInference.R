# load the pancreas dataset
library("tidyr")
library("dplyr")
# retrieve example data from Damond et al. (2019)
spe <- .loadExample()

#make the condition a factor variable
colData(spe)[["patient_stage"]] <- factor(colData(spe)[["patient_stage"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["patient_stage"]] <- relevel(colData(spe)[["patient_stage"]],
                                           "Non-diabetic")
res <- spatialInference(spe, c("alpha"),
                        subsetby = "image_number", fun = "Gest", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1,
                        algorithm = "bam"
)

mdl1 <- res$mdl

resLs <- crossSpatialInference(spe, c("alpha", "acinar"),
                        subsetby = "image_number", fun = "Gcross", marks = "cell_type",
                        rSeq = seq(0, 50, length.out = 50), correction = "rs",
                        sample_id = "patient_id",
                        image_id = "image_number", condition = "patient_stage",
                        ncores = 1,
                        algorithm = "bam"
)

mdl2 <- resLs$alpha_alpha$mdl

test_that("cross function with one element is identical to the single
          function call", {
  expect_true(identical(mdl1$coefficients, mdl2$coefficients))
})

