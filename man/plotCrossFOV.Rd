% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMetric.R
\name{plotCrossFOV}
\alias{plotCrossFOV}
\title{Creates a nXn plot of the cross metrics per sample}
\usage{
plotCrossFOV(sub_fov, theo, correction, x, image_id, ID = NULL)
}
\arguments{
\item{sub_fov}{a subset of the dataframe to the respective fov}

\item{theo}{logical; if the theoretical line should be plotted}

\item{correction}{the border correction to plot}

\item{x}{the x-axis variable to plot}

\item{image_id}{the ID of the image/fov}

\item{ID}{the (optional) ID for plotting combinations}
}
\value{
a ggplot object
}
\description{
Creates a nXn plot of the cross metrics per sample
}
