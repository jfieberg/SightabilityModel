\name{summary.sightest}
\alias{summary.sightest}
\alias{summary.sightest_ratio}
\title{
Summarize sightability estimator
}
\description{
Calculates confidence interval (based on asympototic [normal or log-normal assumption])
}
\usage{
 \method{summary}{sightest}(object,...)
}
\arguments{
  \item{object}{Sightability object, output from call to Sight.Est function.}
  \item{...}{arguments to be passed to or from other methods}
}
\value{
  \item{Nhat or Ratiohat}{Sightability population estimate}
  \item{lcl}{Lower confidence limit}
  \item{ucl}{Upper confidence limit}
}
\author{
John Fieberg and Carl James Schwarz
}

\seealso{
 \code{\link{Sight.Est}}, \code{\link{Sight.Est.Ratio}}
}
\keyword{summary} 
