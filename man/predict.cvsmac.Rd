\name{predict.cvsmac}
\alias{predict.cvsmac}
\title{A function that provides class label prediction and class conditional probability estimation for objects returned by the "cvsmac" function.}
\description{This function provides prediction on a test data set using the obtained classifier from a call of the cvsmac function.}
\usage{
\method{predict}{cvsmac}(object,new.x = NULL,...)
}

\arguments{
  \item{object}{An object returned by a call from the function "cvsmac".}
  \item{new.x}{The predictor matrix used for prediction. If not specified, the training matrix x will be used.}
  \item{...}{Not used.}
}

\value{
  \item{new.x}{The new testing matrix used for prediction.}
  \item{best.lambda}{The sequence of lambda values that have the smallest cross validation error in the training data set.}
  \item{best.beta0}{The beta0 values that correspond to the best.lambda sequence.}
  \item{best.beta}{The beta values that correspond to the best.lambda sequence.}
  \item{best.pred.y}{The predicted labels for new.x.}
  \item{best.pred.prob}{A list of data.frame containing predicted class conditional probabilities. Each data.frame corresponds to a value in the best.lambda object. Each row of the data.frames corresponds to an observation in new.x, and each column represents a class, with the column name as the class label.}
  \item{call}{The function call that returns this result.}
}

\references{
C. Zhang and Y. Liu (2014). Multicategory Angle-based Large-margin Classification. \emph{Biometrika}, 101(3), 625-640. 
}

\author{Chong Zhang, Guo Xian Yau and Yufeng Liu}

\seealso{
\code{\link{cvsmac}, \link{predict.smac}}
}

\examples{
data(ex1.data)
a=cvsmac(ex1.data$ex1.x,ex1.data$ex1.y,loss="p",nlambda=30)
predict(a,ex1.data$ex1.new.x)
}
