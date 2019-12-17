\name{hgm.MLEso3}
\alias{hgm.MLEso3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    The function hgm.MLEso3 constructs the MLE for the Fisher
  distribution on SO(3).
}
\description{
    The function hgm.MLEso3 constructs the maximum likelihood estimator (MLE) for the Fisher
  distribution on SO(3) given the sample mean of a data set as input.
}
\usage{
hgm.MLEso3(Y,X = c(0.01, 0.005, 0),
           ord=21,method="H-BFGS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A 3 by 3 matrix of the sample mean.}
  \item{X}{(optional) A 1 by 3 vector giving the start point of the optimization algoritm. The default is \eqn{(0.01, 0.005, 0)}. See Details section below.}
  \item{ord}{(optional) The order at which to truncate the series expansion for the normalising constant.
  This is used to calculate the starting values of the normalising constant and its gradient. 
  The higher the order the more accurate the result. However, the computational cost increases. The default value is 20.}
  \item{method}{(optional) The optimization algorthim used. The options are 'H-BFGS', using the BFGS algorthim in conjunction with the holonomic gradient method, 'H-Newton', the holonomic Newton method, and 'Asymptotic', which utilises the formula for singular values close to unity. The 'H-BFGS' method is used by default.}
  }
  
\details{
  The MLE is defined as the argmax of the function 
  \deqn{l(X) = d_1*X_1 + d_2*X_2 + d_3*X_3 - log(c(X)).}
  Here X is the 1 by 3 vector of singular values of the matrix parameter, 
  d is the 1 by 3 vector of sign preserving singular values of the sample mean
  and c(X) is the normalising constant.
}
\value{
A list with the components:
\item{parameter}{The maximum likelihood estimator for the 3 by 3 parameter matrix Theta.}
\item{value}{  The maximum log-likelihood value.}

}
\references{
Michael F. Adamer, Andras C. Lorincz, Anna-Laura Sattelberger and Bernd Sturmfels:
Algebraic Analysis of Rotation Data.
ArXiv preprint 1912.00396, 2019,
\url{https://arxiv.org/abs/1912.00396}
}
\author{
Michael F. Adamer
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%%\code{\link{oxm.matrix_r2tfb}}
%}
\examples{
## ===========================================================================
## Example 1. Computing the MLE for the Asteroid Data of Sei et al.
## ===========================================================================
Y <- matrix(c(0.074,0.018,-0.001,0.012,0.003,0.949,0.016,-0.074,0.002),nrow=3)
hgm.MLEso3(Y)

## ===========================================================================
## Example 2. Computing the MLE for the Comet Data of Sei et al.
## ===========================================================================
Y <- matrix(c(0.257,0.158,0.079,0.044,-0.052,0.765,0.189,-0.146,0.004),nrow=3)
hgm.MLEso3(Y,ord=40,method="H-Newton")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Maximum likelihood estimation }
\keyword{ Directional statistics }
\keyword{ Fisher distribution on SO(3)}
\keyword{ Holonomic gradient method }
\keyword{ HGM }
