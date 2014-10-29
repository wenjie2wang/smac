\name{predict.smac}
\alias{predict.smac}
\title{A function that provides class label prediction and class conditional probability estimation for objects returned by the "smac" function.}
\description{This function provides prediction on a test data set using the obtained classification model from a call of the smac function.}
\usage{
\method{predict}{smac}(object,new.x=NULL,lambda=NULL,...)
}
\arguments{
  \item{object}{An object returned by the "smac" function.}
  \item{new.x}{The new predictor matrix. The number and order of predictors in new.x should be the same as those of x, which is used in the smac function. If not specified, the program uses the training x matrix as the prediction object.}
  \item{lambda}{The lambda values on which the user wants to predict. If not specified, the program will use the lambda values returned by the smac function.}
  \item{...}{Not used.}
 }

\value{
  \item{new.x}{The predictor matrix used for prediction.}
  \item{lambda}{The lambda values used for prediction.}
  \item{fitted.beta0}{The predicted intercept for each lambda value.}
  \item{fitted.beta}{A list of predicted matrices of parameters for each lambda value. The matrices in the list correspond to the lambda values in the returned lambda object in orders. Each row of a matrix corresponds to a predictor, and the name of the predictor is recorded as the row name. Note that a predictor does not affect the class label if and only if all elements in that corresponding row are 0.}
  \item{pred.y}{A list of predicted labels for each observation in new.x. A prediction is made for each lambda value in the returned lambda object.} 
  \item{pred.prob}{A list of data.frames containing predicted class conditional probabilities. Each data.frame object corresponds to a value in the lambda object. Each row of the data.frames corresponds to an observation in new.x, and each column represents a class, with the column name as the class label.}
   }
\references{
C. Zhang and Y. Liu (2014). Multicategory Angle-based Large-margin Classification. \emph{Biometrika}, 101(3), 625-640. 
}
\author{Chong Zhang, Guo Xian Yau and Yufeng Liu}

\seealso{\code{\link{smac}}} 

\examples{
data(ex1.data)
a=smac(ex1.data$ex1.x,ex1.data$ex1.y,loss="p",nlambda=30)
predict(a,ex1.data$ex1.new.x)
}
