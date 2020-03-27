\name{summary.sightest}
\alias{summary.sightest}
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
  \item{Nhat}{Sightability population estimate}
  \item{lcl}{Lower confidence limit}
  \item{ucl}{Upper confidence limit}
}
\author{
John Fieberg
}

\seealso{
 \code{\link{Sight.Est}}
}
\keyword{print} 
