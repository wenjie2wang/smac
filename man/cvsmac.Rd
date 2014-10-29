\name{cvsmac}
\alias{cvsmac}
\title{A cross validation function for smac.}
\description{This function is used to perform cross validation based on a training data set to select the best tuning parameter.}
\usage{
cvsmac(x,y,kfold=5,lambda=NULL,nlambda=100,lambda.min=NULL,seed=0,weight=NULL,...)
}

\arguments{
 \item{x}{The usage of this argument is the same as that in the "smac" function.}
 \item{y}{The usage of this argument is the same as that in the "smac" function.}
 \item{kfold}{The number of subsamples for cross validation. This number should be less than the sample size of the class with the fewest obervations. A warning will be given if kfold is too large. Default is 5.}
 \item{lambda}{The usage of this argument is the same as that in the "smac" function.}
 \item{nlambda}{The usage of this argument is the same as that in the "smac" function.}
 \item{lambda.min}{The usage of this argument is the same as that in the "smac" function.}
 \item{seed}{The seed for generating the random split of the training data set. Default is 0.}
 \item{weight}{The usage of this argument is the same as that in the "smac" function.}
 \item{...}{Other arguments used by function "smac" that are not specified.}
}

\value{
 \item{lambda}{The sequence of tuning parameters used in cross validation. Notice that the lambdas will be in a decreasing order.}
 \item{beta0}{The estimated intercepts with respect to each tuning parameter lambda.}
 \item{beta}{The estimated parameters of the predictors with respect to each tuning parameter lambda.}
 \item{error}{The total number of misclassifications with respect to each tuning parameter lambda. If weight is specified, each miscalssification is multiplied by its corresponding weight.}
 \item{best.lambda}{The sequence of lambda values that have the smallest cross validation error.}
 \item{best.beta0}{The intercepts that correspond to best.lambda.}
 \item{best.beta}{The parameters that correspond to best.lambda.}
 \item{model}{The fitted model with respect to each tuning parameter, using the entire training data set.}
 \item{min.error}{The minimum error in cross validation.}
}

\references{
C. Zhang and Y. Liu (2014). Multicategory Angle-based Large-margin Classification. \emph{Biometrika}, 101(3), 625-640. 
}

\author{Chong Zhang, Guo Xian Yau and Yufeng Liu}

\seealso{\code{\link{smac},\link{predict.cvsmac}}}

\examples{
data(ex1.data)
cvsmac(ex1.data$ex1.x,ex1.data$ex1.y,loss="p",nlambda=30)
}
