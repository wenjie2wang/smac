\name{smac}
\alias{smac}
\title{Classification function that provides solution path to Multicategory Angle-based large-margin Classifiers (MAC) with L1 penalty}
\description{A classifier that works under the structure of MAC (Zhang and Liu, 2014) with linear learning and the L1 penalty for variable selection.}
\usage{
smac(x,y,loss=c("logistic", "psvm", "boost"),weight=NULL,nlambda = 100,
lambda.min=ifelse(nobs < np, 0.05, 1e-03),lambda = NULL,
standardize = TRUE, epsilon = 1e-05)
}
\arguments{
  \item{x}{The x matrix/data.frame for the training dataset. Columns represent the covariates, and rows represent the instances. There should be no NA/NaN values in x.}
  \item{y}{The labels for the training dataset.}
  \item{loss}{The binary large margin classifer loss function to be used. By default the program uses the logistic loss. Exponential loss in boosting and squared loss in proximal support vector machines are also available. "l" or "logi" for logistic loss, "b" or "boost" for boosting loss, and "p" or "psvm" for squared loss in proximal support vector machines.}
  \item{weight}{The weight vector for each observation. By default the program uses equal weights for all observations.}
  \item{nlambda}{The number of lambda values in a solution path, if the user does not specify which lambdas to use. Default is 100.}
  \item{lambda.min}{In a classification problem where the user does not provide a list of lambda values, the program will automatically find the smallest lambda value that makes all the estimated parameters 0 as a starting lambda. Then the program will create a solution path for a list of lambda values based on the starting lambda (this starting lambda is in fact the largest lambda in the solution path). This option specifies how small the last lambda is compared to the starting lambda in terms of ratios. By default if the number of observations is larger than the number of parameters, the smallest lambda in the solution path is set to be 1/1,000 of the starting lambda, and is set to be 1/20 otherwise. The program then chooses nlambda's of lambda values between the starting lambda and the last lambda, based on an even split of log(lambda) values.}
  \item{lambda}{The user specified lambda values. If used, the options nlambda and lambda.min will be ignored.}
  \item{standardize}{Whether the input x should be standardized or not. Default is TRUE (standardize).}
  \item{epsilon}{Convergence threshold in coordinate descent circling algorithm. The smaller epsilon is, the more accurate the final model is, and the more time it takes for calculation. Default is 1e-5.}
}

\value{
  \item{All}{All arguments that are used are recorded.}
  \item{k}{Number of classes in the classification problems.}
  \item{x.name}{The column names of x.}
  \item{y.name}{The class names of y.}
  \item{lambda}{The lambda vector of all lambdas in the solution path.}
  \item{beta0}{A list of the intercepts of the classification function. Each vector in the list corresponds to the lambda in the solution path in order.}
  \item{beta}{A list of matrices containing the estimated parameters of the classification function. Each matrix in the list corresponds to the lambda value in the solution path in order. For one single matrix, the rows correspond to a specific predictor, whose name is recorded as the row name. Note that a predictor does not have a significant effect on the label if and only if all elements in its corresponding row are 0.}
  \item{loss}{The loss function used.}
  \item{way}{A numeric value specifying if the user provides the lambda values in the solution path (2), or not (1). This return is mainly used in the prediction function.}
  \item{call}{The call of smac.}
 }
\references{
C. Zhang and Y. Liu (2014). Multicategory Angle-based Large-margin Classification. \emph{Biometrika}, 101(3), 625-640. 
}
\author{Chong Zhang, Guo Xian Yau and Yufeng Liu}

\seealso{\code{\link{predict.smac}}} 

\examples{
data(ex1.data)
smac(ex1.data$ex1.x,ex1.data$ex1.y,loss="p",nlambda=30)
}
