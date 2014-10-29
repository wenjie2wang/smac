\name{ex1.data}
\alias{ex1.data}
\docType{data}
\title{Toy example data set used in smac}
\description{
 This is a simulated toy example data set to illustrate the functionality of smac.}

\usage{data(ex1.data)}
\format{
  Two matrices consisting of predictors, and their corresponding label vectors.
  \describe{
    \item{\code{ex1.x}}{The training matrix consisting of 2 predictors and 100 observations.}
    \item{\code{ex1.y}}{The label for the training data.}
    \item{\code{ex1.new.x}}{The testing matrix consisting of 2 predictors and 100 observations.}
    \item{\code{ex1.new.y}}{The label for the testing data.}
  }
}
\references{
C. Zhang and Y. Liu (2014). Multicategory Angle-based Large-margin Classification. \emph{Biometrika}, 101(3), 625-640. 
}

\examples{
data(ex1.data)
dim(ex1.data$ex1.x)
ex1.data$ex1.y
dim(ex1.data$ex1.new.x)
ex1.data$ex1.new.y
}
