% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal.feature.R
\name{cal.feature}
\alias{cal.feature}
\title{Calculate meta features from RNA-seq and ATAC-seq}
\usage{
cal.feature(Smat, pcs, neigs = 50)
}
\arguments{
\item{Smat}{the similarity matrix.}

\item{pcs}{principal components from RNA-seq.}

\item{neigs}{number of Eigen vectors to calculate.}
}
\value{
Return meta features.
}
\description{
Calculate meta features from RNA-seq and ATAC-seq
}
\examples{
\dontrun{
feature_coassay<-cal.feature(Smat=Smat,
                            pcs=prcomp(t(expr))$x,
                            neigs=50)
}

}
