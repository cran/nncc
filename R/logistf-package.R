#' @title Firth's Bias-Reduced Logistic Regression
#' @name logistf-package
#' @docType package
#' @aliases {logistf}-package
#'
#' @description
#' Fit a logistic regression model using Firth's bias reduction method, equivalent to penalization of the log-likelihood by the Jeffreys prior.
#' Confidence intervals for regression coefficients can be computed by penalized profile likelihood. Firth's method was proposed as ideal
#' solution to the problem of separation in logistic regression. If needed, the bias reduction can be turned off such that ordinary
#' maximum likelihood logistic regression is obtained.
#'
#' @details
#' The package logistf provides a comprehensive tool to facilitate the application of Firth's modified score
#' procedure in logistic regression analysis. It was written on a PC with S-PLUS 4.0, later translated to S-PLUS 6.0, and to \R.
#'
#' Version 1.10 improves on previous versions by the possibility to include case weights and offsets, and better control
#' of the iterative fitting algorithm.
#'
#' Version 1.20 provides a major update in many respects:
#' \enumerate{
#' \item{Many S3Methods have been defined for objects of type logistf, including add1, drop1 and anova methods}
#' \item{New forward and backward functions allow for automated variable selection using penalized likelihood ratio tests}
#' \item{The core routines have been transferred to C code, and many improvements for speed have been done}
#' \item{Handling of multiple imputed data sets: the 'combination of likelihood profiles' (CLIP) method has been implemented, which builds on datasets that were imputed
#' by the package \code{mice}, but can also handle any imputed data.}
#' }
#'
#' The call of the main function of the library follows the structure of the standard functions as lm or glm, requiring a data.frame and a
#' formula for the model specification.  The resulting object belongs to the new class logistf, which includes penalized maximum likelihood
#' (`Firth-Logistic'- or `FL'-type) logistic regression parameters, standard errors, confidence limits, p-values, the value of the maximized
#' penalized log likelihood, the linear predictors, the number of iterations needed to arrive at the maximum and much more.  Furthermore,
#' specific methods for the resulting object are supplied. Additionally, a function to plot profiles of the penalized likelihood function and a
#' function to perform penalized likelihood ratio tests have been included.
#'
#' In explaining the details of the estimation process we follow mainly the description in Heinze & Ploner (2003). In general, maximum likelihood
#' estimates are often prone to small sample bias. To reduce this bias, Firth (1993) suggested to maximize the penalized log likelihood
#' \eqn{\log L(\beta)^* = \log L(\beta) + 1/2 \log |I(\beta)|}{log L(beta)* = log L(beta) + 1/2 log|I(beta)|}, where \eqn{I(\beta)}{I(beta)} is the
#' Fisher information matrix, i. e. minus the second derivative of the log likelihood. Applying this idea to logistic regression, the score
#' function \eqn{U(\beta)}{U(beta)} is replaced by the modified score function
#' \eqn{U(\beta)^* = U(\beta) + a}{U(beta)* = U(beta) + a}, where \eqn{a} has \eqn{r}th entry
#' \eqn{a_r = 0.5tr{I(\beta)^{-1} [dI(\beta)/d\beta_r]}, r = 1,...,k}{a_r = 0.5tr{I(beta)^{-1} [d I(beta)/d beta_r]}, r = 1,...,k}.
#' Heinze and Schemper (2002) give the explicit formulae for \eqn{I(\beta)}{I(beta)}
#' and \eqn{I(\beta)/d \beta_r}{d I(beta)/d beta_r}.
#'
#' In our programs estimation of \eqn{\beta}{beta} is based on a Newton-Raphson
#' algorithm. Parameter values are initialized usually with 0, but in
#' general the user can specify arbitrary starting values.
#'
#' With a starting value of \eqn{\beta^{(0)}}{beta^(0)}, the penalized maximum
#' likelihood estimate \eqn{\beta}{beta} is obtained iteratively:
#' \deqn{\beta^{(s+1)}= \beta^{(s)} + I(\beta^{(s)})^{-1} U(\beta^{(s)})^* }{beta^(s+1)= beta^(s) + I(beta^(s))^{-1} U(beta^(s))* }
#'
#' If the penalized log likelihood evaluated at \eqn{\beta^{(s+1)}} is less
#' than that evaluated at \eqn{\beta^{(s)}}{beta^(s)} , then (\eqn{\beta^{(s+1)}}{beta^(s+1)} is
#' recomputed by step-halving. For each entry \eqn{r} of \eqn{\beta}{beta} with
#' \eqn{r = 1,...,k} the absolute step size \eqn{|\beta_r^{(s+1)}-\beta_r^s|}{|beta_r^(s+1)-beta_r^s|}
#' is restricted to a maximal allowed value \code{maxstep}. These two means should avoid
#' numerical problems during estimation. The iterative process is continued
#' until the parameter estimates converge, i. e., until three criteria are met: the change in log likelihood is less than \code{lconv},
#' the maximum absolute element of the score vector is less than \code{gconv}, the maximum absolute change in beta is less than \code{xconv}.
#' \code{lconv, gconv, xconv} can be controlled by \code{control=logistf.control(lconv=...,}
#' \code{gconv=..., xconv=...)}.
#'
#' Computation of profile penalized likelihood confidence intervals for
#' parameters (\code{logistpl}) follows the algorithm of Venzon and
#' Moolgavkar (1988). For testing the hypothesis of \eqn{\gamma =
#' \gamma_0}{gamma = gamma_0}, let the likelihood ratio statistic
#'
#' \deqn{LR = 2 [ \log L(\gamma, \delta) - \log L(\gamma_0,\delta_{\gamma_0})^*]}{LR = 2 [ log L(gamma, delta) - log L(gamma_0,delta_{gamma_0})*]}
#'
#' where \eqn{(\gamma, \delta)}{(gamma, delta)}  is the joint penalized maximum likelihood estimate of \eqn{\beta=
#' (\gamma,\delta)}{beta=(gamma,delta)}, and \eqn{\delta_{\gamma_0}}{delta_{gamma_0}} is the penalized maximum
#' likelihood estimate of \eqn{\delta}{delta} when  \eqn{\gamma= \gamma_0}{gamma= gamma_0}. The
#' profile penalized likelihood confidence interval is the continuous set
#' of values \eqn{\gamma_0}{gamma_0} for which \eqn{LR} does not exceed the \eqn{(1 - \alpha)100}{(1 - alpha)100}th
#' percentile of the \eqn{\chi^2_1}{chi^2_1}-distribution. The
#' confidence limits can therefore be found iteratively by approximating
#' the penalized log likelihood function in a neighborhood of \eqn{\beta}{beta} by
#' the quadratic function
#' \deqn{ l(\beta+\delta) = l(\beta) + \delta'U^* - 0.5 \delta' I \delta }{ l(beta+delta) = l(beta) + delta'U* - 0.5 delta' I delta }
#' where \eqn{U^* = U(\beta)^*}{U^* = U(beta)*} and \eqn{-I = -I(\beta)}{-I = -I(beta)}.
#'
#' In some situations computation of profile penalized likelihood
#' confidence intervals may be time consuming since the iterative procedure
#' outlined above has to be repeated for the lower and for the upper
#' confidence limits of each of the k parameters. In other problems one may
#' not be interested in interval estimation, anyway. In such cases, the
#' user can request computation of Wald confidence intervals and P-values,
#' which are based on the normal approximation of the parameter estimates
#' and do not need any iterative estimation process. Standard errors
#' \eqn{\sigma_r, r = 1,...,k}{sigma_r, r = 1,...,k}, of the parameter estimates are computed as
#' the roots of the diagonal elements of the variance matrix \eqn{V(\beta) =
#' I(\beta)^{-1}}{V(beta) =  I(beta)^{-1}} . A \eqn{100(1 - \alpha)}{100(1 - alpha)} per cent Wald confidence interval for
#' parameter \eqn{\beta_r}{beta_r} is then defined as \eqn{[\beta_r +
#' \Psi_{\alpha/2}\sigma_r, \beta_r+\Psi_{1-\alpha/2}\sigma_r]}{[beta_r +
#' Psi_{alpha/2}sigma_r, beta_r+Psi_{1-alpha/2}sigma_r]} where
#' \eqn{\Psi_{\alpha}}{Psi_{alpha}} denotes the \eqn{\alpha}{alpha}-quantile of the standard normal
#' distribution function. The adequacy of Wald confidence intervals for
#' parameter estimates should be verified by plotting the profile penalized
#' log likelihood (PPL) function. A symmetric shape of the PPL function
#' allows use of Wald intervals, while an asymmetric shape demands profile
#' penalized likelihood intervals (\cite{Heinze & Schemper (2002)}).  Further documentation
#' can be found in \cite{Heinze & Ploner (2004)}.
#'
#' The latest version now also includes functions to work with multiply imputed data sets, such as generated by the \code{\link{mice}} package.
#' Results on individual fits can be pooled to obtain point and interval estimates, as well as profile likelihood confidence intervals and likelihood
#' profiles in general (Heinze, Ploner and Beyea, 2013).
#'
#' @author Georg Heinze <georg.heinze@meduniwien.ac.at> and Meinhard Ploner
#' @references
#' Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 80, 27--38.
#'
#' Heinze G, Schemper M (2002). A solution to the problem of
#' separation in logistic regression. \emph{Statistics in Medicine} 21: 2409-2419.
#'
#' Heinze G, Ploner M (2003). Fixing the nonconvergence bug in
#' logistic regression with SPLUS and SAS. \emph{Computer Methods and Programs in Biomedicine} 71: 181-187.
#'
#' Heinze G, Ploner M (2004). Technical Report 2/2004: A SAS-macro, S-PLUS library and
#' R package to perform logistic regression without convergence problems. Section of Clinical Biometrics, Department of
#' Medical Computer Sciences, Medical University of Vienna, Vienna, Austria.
#' \url{http://www.meduniwien.ac.at/user/georg.heinze/techreps/tr2_2004.pdf}
#'
#' Heinze G (2006). A comparative investigation of methods for logistic regression
#' with separated or nearly separated data. \emph{Statistics in Medicine} 25: 4216-4226.
#'
#' Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining profile likelihood information from logistic regressions. Statistics in Medicine, to appear.
#'
#' Venzon DJ, Moolgavkar AH (1988). A method for computing profile-likelihood
#' based confidence intervals. \emph{Applied Statistics} 37:87-94.
#'
#' @keywords models regression
#'
#' @importFrom stats add1 anova as.formula binomial coef density drop1 glm lm model.frame model.matrix model.offset model.response model.weights pchisq pnorm prcomp predict qchisq qnorm terms uniroot update vcov
#' @importFrom graphics abline axis grid legend lines mtext par plot points segments title
#' @importFrom utils capture.output head
#' @importFrom stats nobs
#' @importFrom mgcv uniquecombs
#' @importFrom mice complete
#' @importFrom formula.tools lhs.vars
#'
#' @useDynLib nncc
NULL
