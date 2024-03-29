% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epiConv_funs.R
\name{run.epiConv.parallel}
\alias{run.epiConv.parallel}
\title{Calculate similarities between cells in parallel}
\usage{
run.epiConv.parallel(
  mat,
  row_sample = NULL,
  lib_size,
  nbootstrap,
  nsample,
  inf = (-8),
  bin = 1000,
  ncore = 5,
  ...
)
}
\arguments{
\item{mat}{Matrix object constains the peak by cell matrix.}

\item{row_sample}{indices of sampled cells.}

\item{lib_size}{library size of single cells.}

\item{nbootstrap}{number of bootstraps performed.}

\item{nsample}{number of peaks sampled in each bootstrap.}

\item{inf}{a small value used to replace -Inf.}

\item{bin}{Matrix are divided by bin to avoid very long vectors, It does not affect the results.}

\item{ncore}{number of threads.}

\item{...}{Arguments passed to big.matrix when creating the similarity matrix. You can specify a backup file through "backingfile" and "descriptorfile". If not specified, the function will create a in-memory backup file.}
}
\value{
Return a similarity matrix.
}
\description{
Calculate similarities between cells in parallel
}
\examples{
\dontrun{
Smat<-run.epiConv.parallel(mat=mat,
                           row_sample=1:ncol(mat),
                           lib_size=colSums(mat),
                           nbootstrap=15,
                           nsample=floor(nrow(mat)*0.2),
                           bin=5000,
                           inf=(-8),
                           ncore=5,
                           backingfile="backupfile",
                           descriptorfile="backupfile.descriptor")
}

}
