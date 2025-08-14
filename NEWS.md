# spatialFDA 1.1.7
* Fixing error in `plotCrossHeatmap` when the model is NULL.

# spatialFDA 1.1.6
* Error handling in `spatialInference` when one condition has no images with curves.
* Wrote as well a new test for this case

# spatialFDA 1.1.5
* New overview vignette, showing the usability of `spatialInference`
* Plotting of a heatmap from extracted from `crossSpatialInference`
* Adding of the package sticker to the README

# spatialFDA 1.1.4
* `functionalGAM` accepts flexible naming of the intercept column.

# spatialFDA 1.1.3
* changed the default of the GAM estimation to be a Gaussian with log link for
positivty of the response

# spatialFDA 1.1.2
* wrote a new convenience function `crossSpatialInference` which makes the 
estimation across a range of cell types easier.

# spatialFDA 1.1.1
* wrote a new convenience function `spatialInference` which makes the FDA estimation
of the spatial statistics functions simpler
* in order to do this various changes to `prepData` were introduced
* furthermore, small syntactic changes for the retrieval of the functional
intercept were added to `functionalGam` and `plotMdl`
* The vignette was simplified to reflect the changes in `prepData`

# spatialFDA 0.99.14
* changed the example in the Vignette to scaled t-dist. with log link due to 
lower AIC than Gaussian.

# spatialFDA 0.99.13
* changed the example in the vignette to Gaussian with log link and smooth
random intercept.

# spatialFDA 0.99.12
* Changed indexing of spatial coordinates to positional index. Names can differ
between objects, e.g. "x" or "coord_x" etc. therefore, positional index 
generalises better.

# spatialFDA 0.99.11
* major changes in fixing the levels of the mark factors in the point pattern.
Prior to this fix, it could happen that the marks in `selection` do not correspond
to the marks in the point pattern and therefore give the inverse result.

# spatialFDA 0.99.10
* minor bug fix 

# spatialFDA 0.99.9
* remove `imcdatasets` dependency and instead download example data directly
from `ExperimentHub`
* added small reader function `.loadExample` to download the data mentioned 
above

# spatialFDA 0.99.8
* more packages explicitly named
* vignette shows use of PCA based random errors and exchanged one sample

# spatialFDA 0.99.7
* explicit package naming for functions to make compatible with linux distributions
* possible to pass `...` parameters to fPCA method

# spatialFDA 0.99.6
* small error correction in show fPCA method

# spatialFDA 0.99.5
* change in how failures in `calcMetricPerFov` are handled. Works now as well
if `rSeq = NULL`

# spatialFDA 0.99.4
* variable `ID` is no longer created in `calcMetricPerFov`
* adjusted examples and vignette to the above
* added option to extract gene expression directly in `calcMetricPerFov` if
`continuous = TRUE`
* small bug fix in `calcMetricPerFov` to pass `...` in all cases
* added variables `legend.position` and `ncol` to `plotMetricPerFov`

# spatialFDA 0.99.3
* adjusted `by` variable to work with length 1
* adjusted `...` in `plotMetricPerFov` to be passed only to `geom_line`

# spatialFDA 0.99.3
* adjust examples to build in under 10 min

# spatialFDA 0.99.2
* Rewrote vignette with smaller dataset
* Implemented quasi likelihood family in vignette to improve variance estimation
* Adjustd calcMetricPerFov.R to be able to calculate cross functions

# spatialFDA 0.99.1
* Adressed Bioconductor review
* Added unit tests
* Added functional boxplots
* Added diagnostic plots to vignette

# spatialFDA 0.99.0
* Initial Bioconductor submission.
