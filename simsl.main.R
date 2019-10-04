#' Single-index models with a surface-link (main function)
#'
#' \code{simsl} is the wrapper function for fitting a single-index model with a surface-link (SIMSL).
#' The function estimates a linear combination (a single-index) of baseline covariates X, and models a nonlinear interactive structure between the single-index and a treatment variable defined on a continuum, via estimating a smooth link function on the index-treatment domain.
#'
#' SIMSL captures the effect of covariates via a single-index and their interaction with the treatment via a 2-dimensional smooth link function.
#' Interaction effects are determined by shapes of the link surface.
#' The SIMSL allows comparing different individual treatment levels and constructing individual treatment rules,
#' as functions of a biomarker signature (single-index), efficiently utilizing information on patient’s characteristics.
#' The resulting \code{simsl} object can be used to estimate an optimal dose rule for a new patient with baseline clinical information.
#'
#' @param y   a n-by-1 vector of treatment outcomes; y is assumed to follow an exponential family distribution; any distribution supported by \code{mgcv::gam}.
#' @param A  a n-by-1 vector of treatment variable; each element is assumed to take a value on a continuum.
#' @param X   a n-by-p matrix of baseline covarates.
#' @param family  specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; the defult is "gaussian"; can be any family supported by \code{mgcv::gam}.
#' @param mu.hat   a n-by-1 vector of the fitted main effect term of the model provided by the user; the defult is \code{NULL} and it is taken as a vector of zeros; the optimal choice for this vector is h(E(y|X)), where h is the canonical link function.
#' @param beta.ini  an initial solution of \code{beta.coef}; a p-by-1 vector; the defult is \code{NULL}.
#' @param beta.ini.gam  if \code{TRUE}, employ a \code{mgcv::gam} smooth function representation of the variable A effect when inializing \code{beta.coef}; otherwise use a linear model representation for the A effect at initialization.
#' @param ind.to.be.positive  for identifiability of the solution \code{beta.coef}, we restrict the jth component of \code{beta.coef} to be positive; by default \code{j=1}.
#' @param bs type of basis for representing the treatment-specific smooths; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines)
#' @param k  basis dimension; the same number (k) is used for all treatment groups, however, the smooths of different treatments have different roughness parameters.
#' @param knots  a list containing user specified knot values to be used for basis construction (for the treatment and the index variables, respectively).
#' @param sp  a vector of smoothing parameters associated with the 2-dimensional smooth
#' @param method  the smoothing parameter estimation method; "GCV.Cp" to use GCV for unknown scale parameter and Mallows' Cp/UBRE/AIC for known scale; any method supported by \code{mgcv::gam} can be used.
#' @param pen.order 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{beta.coef}.
#' @param lambda  a regularziation parameter associated with the penalized LS of \code{beta.coef}.
#' @param max.iter  an integer specifying the maximum number of iterations for \code{beta.coef} update.
#' @param eps.iter a value specifying the convergence criterion of algorithm.
#' @param center.X   if \code{TRUE}, center X to have zero mean.
#' @param scale.X    if \code{TRUE}, scale X to have unit variance.
#' @param si.main.effect   if \code{TRUE}, once the convergece in the estimates of \code{beta.coef} is reached, include the main effect associated with the fitted single-index (beta.coef'X) to the final surface-link estimate.
#' @param bootstrap if \code{TRUE}, compute bootstrap confidence intervals for the single-index coefficients, \code{beta.coef}; the default is \code{FALSE}.
#' @param boot.conf  a value specifying the confidence level of the bootstrap confidence intervals; the defult is \code{boot.conf = 0.95}.
#' @param nboot  when \code{bootstrap=TRUE}, a value specifying the number of bootstrap replications.
#' @param seed  when  \code{bootstrap=TRUE}, randomization seed used in bootstrap resampling.
#' @param trace.iter if \code{TRUE}, trace the estimation process and print the differences in \code{beta.coef}.
#'
#' @return a list of information of the fitted SIMSL including
#'  \item{beta.coef}{ the estimated single-index coefficients.} \item{g.fit}{a \code{mgcv:gam} object containing information about the estimated 2-dimensional link function.} \item{beta.ini}{the initial value used in the estimation of \code{beta.coef}} \item{beta.path}{solution path of \code{beta.coef} over the iterations} \item{d.beta}{records the change in \code{beta.coef} over the solution path, \code{beta.path}} \item{X.scale}{sd of pretreatment covariates X} \item{X.center}{mean of pretreatment covariates X} \item{A.range}{range of the observed treatment variable A} \item{p}{number of baseline covariates X} \item{n}{number of subjects} \item{boot.ci}{\code{boot.conf}-level bootstrap CIs (LB, UB) associated with \code{beta.coef}} \item{boot.mat}{a (nboot x p) matrix of bootstrap estimates of  \code{beta.coef}}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import mgcv stats
#' @seealso \code{pred.simsl},  \code{fit.simsl}
#' @export
#'
#' @examples
#'

#'set.seed(1234)
#'n <- 200
#'n.test <- 500
#'
#'## simulation 1
#'# generate training data
#'p <- 30
#'X <- matrix(runif(n*p,-1,1),ncol=p)
#'A <- runif(n,0,2)
#'f_opt <- 1 + 0.5*X[,2] + 0.5*X[,1]
#'mu <- 8 + 4*X[,1] - 2*X[,2] - 2*X[,3] - 25*((f_opt-A)^2)
#'y <- rnorm(length(mu),mu,1)
#'# fit SIMSL
#'simsl.obj <- simsl(y=y, A=A, X=X)
#'
#'# generate testing data
#'X.test <- matrix(runif(n.test*p,-1,1),ncol=p)
#'A.test <- runif(n.test,0,2)
#'f_opt.test <- 1 + 0.5*X.test[,2] + 0.5*X.test[,1]
#'pred <- pred.simsl(simsl.obj, newx= X.test)  # make prediction based on the estimated SIMSL
#'value <- mean(8 + 4*X.test[,1] - 2*X.test[,2] - 2*X.test[,3] - 25*((f_opt.test- pred$trt.rule)^2))
#'value  # the "value" of the estimated treatment rule; the "oracle" value is 8.

#'
#'## simulation 2
#'p <- 10
#'# generate training data
#'X = matrix(runif(n*p,-1,1),ncol=p)
#'A = runif(n,0,2)
#'f_opt = I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) +
#'  1.2*I(X[,1] < -0.5) + X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
#'mu =   8 + 4*cos(2*pi*X[,2]) - 2*X[,4] - 8*X[,5]^3 - 15*abs(f_opt-A)
#'y = rnorm(length(mu),mu,1)
#'Xq <- cbind(X, X^2)  # include a quadratic term
#'# fit SIMSL
#'simsl.obj <- simsl(y=y, A=A, X=Xq)
#'
#'# generate testing data
#'X.test = matrix(runif(n.test*p,-1,1),ncol=p)
#'A.test = runif(n.test,0,2)
#'f_opt.test = I(X.test[,1] > -0.5)*I(X.test[,1] < 0.5)*0.6 + 1.2*I(X.test[,1] > 0.5) +
#'  1.2*I(X.test[,1] < -0.5) + X.test[,4]^2 + 0.5*log(abs(X.test[,7])+1) - 0.6
#'Xq.test <- cbind(X.test, X.test^2)
#'pred <- pred.simsl(simsl.obj, newx= Xq.test)  # make prediction based on the estimated SIMSL
#'value <- mean(8 + 4*cos(2*pi*X.test[,2]) - 2*X.test[,4] - 8*X.test[,5]^3 -
#'               15*abs(f_opt.test-pred$trt.rule))
#'value  # the "value" of the estimated treatment rule; the "oracle" value is 8.
#'
#'
#'
#'
#'### air pollution data application
#'data(chicago); head(chicago)
#'chicago <- chicago[,-3][complete.cases(chicago[,-3]), ]
#'#plot(chicago$death)
#'#chicago$death[2856:2859]
#'chicago <- chicago[-c(2856:2859), ]  # get rid of the gross outliers in y
#'#plot(chicago$pm10median)
#'chicago <- chicago[-which.max(chicago$pm10median), ]  # get rid of the gross outliers in x
#'
#'## create lagged variables
#'lagard <- function(x,n.lag=5) {
#'  n <- length(x); X <- matrix(NA,n,n.lag)
#'  for (i in 1:n.lag) X[i:n,i] <- x[i:n-i+1]
#'  X
#'}
#'chicago$pm10 <- lagard(chicago$pm10median)
#'chicago <- chicago[complete.cases(chicago), ]
#'# create season varaible
#'chicago$time.day <- round(chicago$time %%  365)
#'
#'# fit SIMSL for modeling the season-by-pm10 interactions on their effects on outcomes
#'simsl.obj <- simsl(y = chicago$death, A = chicago$time.day, X=chicago[,7], bs= c("cc", "ps"),
#'                   beta.ini.gam = TRUE, family=poisson(), method = "REML")
#'simsl.obj$beta.coef  # the estimated single-index coefficients
#'summary(simsl.obj$g.fit)
#'#simsl.obj.boot <- simsl(y = chicago$death, A = chicago$time.day, X=chicago[,7],
#'#                        bs= c("cc", "ps"), family=poisson(), beta.ini.gam = TRUE,
#'#                        method = "REML", bootstrap = TRUE, nboot=5)  # nboot =500
#'#simsl.obj.boot$boot.ci
#'
#'
#'par(mfrow=c(1,3), mar=c(5.1,4.7,4.1,2.1))
#'additive.fit  <- mgcv::gam(chicago$death ~
#'                        s(simsl.obj$g.fit$model[,3], k=8, bs="ps") +
#'                        s(chicago$time.day, k=8, bs="cc"),
#'                        family = poisson(), method = "REML")
#'plot(additive.fit, shift= additive.fit$coefficients[1], select=2,
#'     ylab= "Linear predictor", xlab= "A", main = expression(paste("Individual A effect")))
#'plot(additive.fit, shift= additive.fit$coefficients[1], select = 1,
#'     xlab= expression(paste(beta*minute,"x")), ylab= " ",
#'     main = expression(paste("Individual ", beta*minute,"x effect")))
#'par(mar=c(2.1,2.5,4.1,2.1))
#'mgcv::vis.gam(simsl.obj$g.fit, view=c("A","single.index"), theta=-135, phi = 30,color="heat", se=1,
#'        ylab = "single-index", zlab = " ", main=expression(paste("Interaction surface ")))
#'
#'
#'
#'
#'### Warfarin data application
#'data(warfarin)
#'X <- warfarin$X
#'A <- warfarin$A
#'y <- -abs(warfarin$INR - 2.5)  # the target INR is 2.5
#'X[,1:3] <- scale(X[,1:3]) # standardize continuous variables
#'
#'## Estimate the main effect, using an additive model for continous variables and
#'## a linear model for the indicator variables
#'mu.fit <- mgcv::gam(y-mean(y)  ~ X[, 4:13] +
#'                    s(X[,1], k=5, bs="ps")+
#'                    s(X[,2], k=5, bs="ps") +
#'                    s(X[,3], k=5, bs="ps"), method="REML")
#'summary(mu.fit)
#'mu.hat <- predict(mu.fit)
#'# fit SIMSL (we do not scale/center X for the interpretabilty of the indicator variables in X).
#'simsl.obj <- simsl(y, A, X, mu.hat=mu.hat, scale.X = FALSE, center.X=FALSE, method="REML")
#'simsl.obj$beta.coef
#'## can also compute bootstrap CIs for the single-index coefficients (beta.coef)
#'#simsl.obj.boot <- simsl(y, A, X, mu.hat=mu.hat, scale.X=FALSE, center.X=FALSE,
#'#                        bootstrap = TRUE, nboot=5, method="REML")  # nboot = 500
#'#simsl.obj.boot$boot.ci
#'
#'
#'###
#'par(mfrow=c(1,3), mar=c(5.1,4.7,4.1,2.1))
#'additive.fit  <- mgcv::gam(y-mu.hat ~
#'                         s(A, k=8, bs="ps") +
#'                         s(simsl.obj$g.fit$model[,3], k=8, bs="ps"),
#'                       method = "REML" )
#'plot(additive.fit, shift= additive.fit$coefficients[1], select=1,
#'     ylab= "Y", main = expression(paste("Individual A effect")))
#'plot(additive.fit, shift= additive.fit$coefficients[1], select=2,
#'     xlab= expression(paste(beta*minute,"x")), ylab= " ",
#'     main = expression(paste("Individual ", beta*minute,"x effect")))
#'par(mar=c(2.1,2.5,4.1,2.1))
#'mgcv::vis.gam(simsl.obj$g.fit, view=c("A","single.index"), theta=55, phi = 30,color="heat", se=1,
#'        ylab = "single-index", zlab = "Y", main=expression(paste("Interaction surface ")))
#'

simsl <- function(y,  # a n-by-1 vector of treatment outcomes; y is assumed to follow an exponential family distribution; any distribution supported by \code{mgcv::gam}.
                  A,  # a n-by-1 vector of treatments; the ith element represents treatment dose applied to the ith individual.
                  X,  # a n-by-p matrix of covarates.
                  mu.hat = NULL,  # a n-by-1 vector for efficinecy augmentation provided by the user; the defult is \code{NULL}; the optimal choice for this vector is h(E(y|X)), where h is the canonical link function.
                  family = "gaussian",  # specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; the defult is "gaussian"; can be any family supported by \code{mgcv::gam}.
                  bs =c("ps", "ps"),  # type of basis (for the treatment and the index variables, respectively) for representing the 2-dimensional link function on the treatment-index domain; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines)
                  k = c(8, 8),    # basis dimension (for the treatment and the index variables, respectively).
                  knots= NULL,  # a list containing user specified knot values to be used for basis construction (for the treatment and the index variables, respectively).
                  sp = NULL,  # a vector of smoothing parameters associated with the smooths
                  method = "GCV.Cp",  # the smoothing parameter estimation method; "GCV.Cp" to use GCV for unknown scale parameter and Mallows' Cp/UBRE/AIC for known scale; any method supported by \code{mgcv::gam} can be used.
                  beta.ini = NULL,  # an initial solution of \code{beta.coef}; a p-by-1 vector; the defult is \code{NULL}.
                  beta.ini.gam = FALSE,
                  ind.to.be.positive = 1,  # or identifiability of the solution \code{beta.coef}, we restrict the jth component of \code{beta.coef} to be positive; by default \code{j=1}.
                  pen.order = 0,   # 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{beta.coef}.
                  lambda = 0,   # a regularziation parameter associated with the penalized LS of \code{beta.coef}.
                  max.iter = 30,   # an integer specifying the maximum number of iterations for \code{beta.coef} update.
                  eps.iter =10^{-2},  # a value specifying the convergence criterion of algorithm.
                  trace.iter = TRUE,  #  if \code{TRUE}, trace the estimation process and print the differences in \code{beta.coef}.
                  center.X=TRUE,  # if \code{TRUE}, center X to have zero mean.
                  scale.X=TRUE,   # if \code{TRUE}, scale X to have unit variance.
                  si.main.effect = TRUE,  # if \code{TRUE}, once the convergece in the estimates of \code{beta.coef} is reached, include the main effect of the estimated single-index (\code{beta.coef}'X) to the surface-link fit, using an unconstrained tensor product basis.
                  bootstrap = FALSE,  # if \code{TRUE}, compute bootstrap confidence intervals for the single-index coefficients, \code{beta.coef}; the default is \code{FALSE}.
                  nboot= 200,  # when \code{bootstrap=TRUE}, a value specifying the number of bootstrap replications.
                  boot.conf = 0.95,   # when \code{bootstrap=TRUE}, a value specifying the confidence level of bootstrap CIs.
                  seed= 1357)   # if \code{TRUE}, trace the estimation process and print the differences in \code{beta.coef}.
{

  simsl.obj <- fit.simsl(y=y, A=A, X=X, mu.hat=mu.hat, family=family,
                         bs =bs, k = k, knots=knots, sp=sp, method=method,
                         beta.ini = beta.ini, beta.ini.gam = beta.ini.gam,
                         ind.to.be.positive=ind.to.be.positive,
                         pen.order = pen.order, lambda = lambda,
                         max.iter = max.iter, trace.iter=trace.iter,
                         center.X= center.X, scale.X= scale.X,
                         si.main.effect= si.main.effect)
  boot.mat = boot.ci <- NULL
  if(bootstrap){
    set.seed(seed)
    indices <- 1:simsl.obj$n
    if(is.null(mu.hat)) mu.hat <- rep(0,simsl.obj$n)
    boot.mat <- matrix(0, nboot, simsl.obj$p)
    for(i in 1:nboot){
      boot.indices <- sample(indices, simsl.obj$n, replace = TRUE)
      boot.mat[i,] <- fit.simsl(y=y[boot.indices], A = A[boot.indices], X = X[boot.indices,],
                                mu.hat = mu.hat[boot.indices],
                                family=family, bs =bs, k = k, knots=knots, #sp = simsl.obj$g.fit$sp, #sp=sp, method=method,
                                beta.ini = beta.ini, beta.ini.gam= beta.ini.gam, ind.to.be.positive=ind.to.be.positive,
                                pen.order = pen.order, lambda = lambda, max.iter = max.iter, trace.iter=trace.iter,
                                center.X= center.X, scale.X= scale.X,
                                si.main.effect= si.main.effect)$beta.coef
    }

    boot.mat.abs <- abs(boot.mat)
    var.t0.abs <- apply(boot.mat.abs, 2, var)
    boot.ci <- cbind(simsl.obj$beta.coef-qnorm((1+boot.conf)/2)*sqrt(var.t0.abs),
                     simsl.obj$beta.coef+qnorm((1+boot.conf)/2)*sqrt(var.t0.abs))

    boot.ci <- cbind(simsl.obj$beta.coef, boot.ci, (boot.ci[,1] > 0 | boot.ci[,2] < 0) )
    colnames(boot.ci) <- c("coef", "LB", "UB", " ***")
    rownames(boot.ci) <- colnames(X)
  }
  simsl.obj$boot.mat <- boot.mat
  simsl.obj$boot.ci <- boot.ci

  return(simsl.obj)
}










#' Single-index models with a surface-link (workhorse function)
#'
#' \code{fit.simml} is the workhorse function for Single-index models with a surface-link (SIMSL).
#'
#' The function estimates a linear combination (a single-index) of covariates X, and captures a nonlinear interactive structure between the single-index and the treatment defined on a continuum via a smooth surface-link on the index-treatment domain.
#'
#' SIMSL captures the effect of covariates via a single-index and their interaction with the treatment via a 2-dimensional smooth link function.
#' Interaction effects are determined by shapes of the link function.
#' The model allows comparing different individual treatment levels and constructing individual treatment rules,
#' as functions of a biomarker signature (single-index), efficiently utilizing information on patient’s characteristics.
#' The resulting \code{simsl} object can be used to estimate an optimal dose rule for a new patient with pretreatment clinical information.
#'
#' @param y   a n-by-1 vector of treatment outcomes; y is assumed to follow an exponential family distribution; any distribution supported by \code{mgcv::gam}.
#' @param A  a n-by-1 vector of treatment variable; each element represents one of the L(>1) treatment conditions; e.g., c(1,2,1,1,1...); can be a factor-valued.
#' @param X   a n-by-p matrix of pre-treatment covarates.
#' @param family  specifies the distribution of y; e.g., "gaussian", "binomial", "poisson"; the defult is "gaussian"; can be any family supported by \code{mgcv::gam}.
#' @param mu.hat   a n-by-1 vector for efficinecy augmentation provided by the user; the defult is \code{NULL}; the optimal choice for this vector is h(E(y|X)), where h is the canonical link function.
#' @param beta.ini  an initial solution of \code{beta.coef}; a p-by-1 vector; the defult is \code{NULL}.
#' @param beta.ini.gam  if \code{TRUE}, employ a \code{mgcv::gam} smooth function representation of the variable A effect when inializing \code{beta.coef}; otherwise use a linear model representation for the A effect at initialization.
#' @param ind.to.be.positive  for identifiability of the solution \code{beta.coef}, we restrict the jth component of \code{beta.coef} to be positive; by default \code{j=1}.
#' @param knots  a list containing user specified knot values to be used for basis construction (for the treatment and the index variables, respectively).
#' @param sp  a vector of smoothing parameters associated with the 2-dimensional smooth
#' @param method  the smoothing parameter estimation method; "GCV.Cp" to use GCV for unknown scale parameter and Mallows' Cp/UBRE/AIC for known scale; any method supported by \code{mgcv::gam} can be used.
#' @param bs type of basis for representing the treatment-specific smooths; the defult is "ps" (p-splines); any basis supported by \code{mgcv::gam} can be used, e.g., "cr" (cubic regression splines)
#' @param k  basis dimension; the same number (k) is used for all treatment groups, however, the smooths of different treatments have different roughness parameters.
#' @param pen.order 0 indicates the ridge penalty; 1 indicates the 1st difference penalty; 2 indicates the 2nd difference penalty, used in a penalized least squares (LS) estimation of \code{beta.coef}.
#' @param lambda  a regularziation parameter associated with the penalized LS of \code{beta.coef}.
#' @param max.iter  an integer specifying the maximum number of iterations for \code{beta.coef} update.
#' @param eps.iter a value specifying the convergence criterion of algorithm.
#' @param center.X   if \code{TRUE}, center X to have zero mean.
#' @param scale.X    if \code{TRUE}, scale X to have unit variance.
#' @param si.main.effect   if \code{TRUE}, once the convergece in the estimates of \code{beta.coef} is reached, include the main effect associated with the fitted single-index (beta.coef'X) to the final surface-link estimate.
#' @param trace.iter if \code{TRUE}, trace the estimation process and print the differences in \code{beta.coef}.
#'
#' @return a list of information of the fitted SIMSL including
#'  \item{beta.coef}{ the estimated single-index coefficients.} \item{g.fit}{a \code{mgcv:gam} object containing information about the estimated 2-dimensional link function.} \item{beta.ini}{the initial value used in the estimation of \code{beta.coef}} \item{beta.path}{solution path of \code{beta.coef} over the iterations} \item{d.beta}{records the change in \code{beta.coef} over the solution path, \code{beta.path}} \item{X.scale}{sd of pretreatment covariates X} \item{X.center}{mean of pretreatment covariates X} \item{A.range}{range of the observed treatment variable A} \item{p}{number of baseline covariates X} \item{n}{number of subjects}
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @import mgcv
#' @seealso \code{pred.simsl},  \code{fit.simsl}
#' @export
#'
fit.simsl  <- function(y, A, X,
                       mu.hat = NULL,
                       family = "gaussian",
                       bs =c("ps", "ps"),
                       k = c(8, 8),
                       knots = NULL,   # a list consist of the knots for the variables A and X, respectively.
                       sp = NULL,      # a vector of smoothing parameters associated with the smooths
                       method = "GCV.Cp",
                       beta.ini = NULL,
                       beta.ini.gam = FALSE,
                       ind.to.be.positive =1,
                       pen.order = 0,
                       lambda = 0,
                       max.iter = 30,
                       eps.iter = 0.01,
                       trace.iter = TRUE,
                       scale.X = TRUE,
                       center.X = TRUE,
                       si.main.effect = TRUE)
{

  y <- as.vector(y)
  n <- length(y)
  p <- ncol(X)
  A.range <- range(A, na.rm = TRUE)  # the observed range of A

  ## Center and scale X
  Xc <- scale(X, center = center.X, scale = scale.X)
  X.center <- attr(Xc, "scaled:center")
  X.scale <- attr(Xc, "scaled:scale")

  ## If not provided by the user, the efficiency augmentation vector (corresponding to the X main effect) is set to be a zero vector.
  if(is.null(mu.hat)) mu.hat <- rep(0, n)

  ## Specify a penalty matrix associated with the penalized least squares for estimating beta.coef.
  D <- diag(p);  if(pen.order != 0)  for(j in 1:pen.order) D <- diff(D);
  Pen <- sqrt(lambda)*D

  ## Initialization for beta.coef
  if(is.null(beta.ini)){
    Ac <- A-mean(A)
    if(beta.ini.gam){
      tmp <- gam(y~ s(Ac, k=k[1], bs=bs[1]) + Xc + Ac:Xc, family = family)$coef
      tmp.ind <- grep(":Ac" ,names(tmp))
    }else{
      tmp <- glm(y~  Ac + Xc + Ac:Xc, family = family)$coef
      tmp.ind <- grep("Ac:" ,names(tmp))
    }
    beta.ini <- tmp[tmp.ind]
    beta.ini[which(is.na(beta.ini))] <- 0
  }

  beta.coef <- beta.ini/sqrt(sum(beta.ini^2))  # enforce unit L2 norm
  if(beta.coef[ind.to.be.positive] < 0) beta.coef <- -1*beta.coef       # for the (sign) identifiability
  # initialize the single-index
  single.index <- as.vector(Xc %*% beta.coef)
  # initialize the surface link (that consists of the A main effect + the A-by-X pure interaction effect)
  if(si.main.effect){
    g.fit <- gam(y ~ s(A, bs=bs[1], k=k[1]) + s(single.index, bs=bs[2], k=k[2]) +
                   ti(A, single.index, bs=bs, k=k), knots = knots, family =  family, gamma=1.4, sp=sp, method=method)
  }else{
    g.fit <- gam(y ~ s(A, bs=bs[1], k=k[1]) +
                   ti(A, single.index, bs=bs, k=k), knots = knots, family =  family, gamma=1.4, sp=sp, method=method)
  }

  ## Compute the adjusted responses (adjusted, for the nonlinearity associated with the GLM canonical link)
  # obtain the 1st derivative of the inverse canonical link, w.r.t. the "linear" predictor
  h.prime <- g.fit$family$mu.eta(predict.gam(g.fit, type="link"))
  adjusted.responses <- (y - predict.gam(g.fit, type = "response")) /h.prime + predict.gam(g.fit, type="link")
  # take the 1st deriavative of the treatment-specific smooths, w.r.t. the single.index.
  g.der <- der.link(g.fit)

  beta.path <- beta.coef
  d.beta <- NULL

  ## Start iteration
  for (it in 2:max.iter)
  {
    ## Update beta.coef and intercept through lsfit
    # adjusted responses, adjusted for the nonlinearity associated with the treatment-specific smooths
    y.star    <- adjusted.responses - predict.gam(g.fit, type="link")  + g.der*single.index - mu.hat
    # adjusetd covariates, adjusted for the nonlinearity of the the treatment smooths
    X.tilda   <- diag(g.der) %*% Xc
    nix       <- rep(0, nrow(D))
    X.p       <- rbind(X.tilda, Pen)
    y.p       <- c(y.star, nix)
    # perform a (penalized) WLS
    beta.fit   <- stats::lsfit(X.p, y.p, wt =c(g.fit$weights, (nix + 1)))
    # for the identifiability
    beta.new <- beta.fit$coef[-1]/sqrt(sum(beta.fit$coef[-1]^2))
    if(beta.new[ind.to.be.positive] < 0) beta.new <- -1*beta.new
    beta.path <- rbind(beta.path, beta.new)

    ## Check the convergence of alpha
    d.beta   <- c(d.beta, sum((beta.new-beta.coef)^2))

    if(trace.iter){
      cat("iter:", it, " "); cat(" difference in alpha: ", d.beta[(it-1)], "\n")
    }
    if (d.beta[(it-1)] < eps.iter)
      break
    beta.coef <- beta.new
    single.index <- as.vector(Xc %*% beta.coef)

    ## Given a single.index, estimate the link surface (subject to the interaction effect identifiability constraint).
    if(si.main.effect){
      g.fit <- gam(y ~ s(A, bs=bs[1], k=k[1]) + s(single.index, bs=bs[2], k=k[2]) +
                     ti(A, single.index, bs=bs, k=k), knots = knots, family =  family, gamma=1.4, sp=sp, method=method)
    }else{
      g.fit <- gam(y ~ s(A, bs=bs[1], k=k[1]) +
                     ti(A, single.index, bs=bs, k=k), knots = knots, family =  family, gamma=1.4, sp=sp, method=method)
    }

    ## compute the adjusted responses (adjusted, for the nonlinearity associated with the GLM canonical link)
    # obtain the 1st derivative of the inverse canonical link, w.r.t. the "linear" predictor
    h.prime <- g.fit$family$mu.eta(predict.gam(g.fit, type="link"))
    adjusted.responses <- (y - predict.gam(g.fit, type = "response")) /h.prime + predict.gam(g.fit, type="link")
    # take the 1st deriavative of the treatment-specific smooths, w.r.t. the single.index.
    g.der <- der.link(g.fit)
  }

  if(si.main.effect){
    g.fit <- gam(y ~ te(A, single.index, bs=bs, k=k), knots = knots, family =  family, gamma=1.4, sp=sp, method=method)
  }


  results <- list(beta.coef = round(beta.coef,3),
                  beta.ini = beta.ini, d.beta=d.beta, beta.path=beta.path,
                  g.fit= g.fit,
                  beta.fit=beta.fit,
                  X.scale=X.scale, X.center = X.center,
                  A.range=A.range, p=p, n=n, bs=bs, k=k)
  class(results) <- c("simsl", "list")

  return(results)
}


#' A subfunction used in estimation
#'
#' This function computes the 1st derivative of the surface-link function with respect to the argument associated with the pure interaction effect term of the smooth, using finite difference.
#'
#' @param g.fit  a \code{mgcv::gam} object
#' @param arg.number the argument of \code{g.fit} that is taken derivative with respect to. The default is \code{arg.number=2} (i.e., take deriviative with respect to the single-index).
#' @param eps a small finite difference used in numerical differentiation.
#' @seealso \code{fit.simsl}, \code{simsl}
#'
der.link <- function(g.fit, arg.number=2, eps=10^(-6))
{
  m.terms <- attr(stats::terms(g.fit), "term.labels")
  newD <- stats::model.frame(g.fit)[, m.terms, drop = FALSE]
  newDF <- data.frame(newD)  # needs to be a data frame for predict
  X0 <- predict.gam(g.fit, newDF, type = "lpmatrix")
  newDF[,arg.number]  <- newDF[,arg.number] + eps
  X1 <- predict.gam(g.fit, newDF, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xi <- Xp * 0
  want <- grep("A,single.index", colnames(X1))  # take only the pure interaction effect term
  Xi[, want] <- Xp[, want]
  g.der  <- as.vector(Xi %*% stats::coef(g.fit))  # the first derivative of eta
  return(g.der)
}


#' SIMSL prediction function
#'
#' This function makes predictions from an estimated SIMSL, given a (new) set of covariates.
#' The function returns a set of predicted outcomes given the treatment values in a dense grid of treatment levels for each individual, and a recommended treatment level (assuming a larger value of the outcome is better).
#'
#' @param simsl.obj  a \code{simsl} object
#' @param newx a (n-by-p) matrix of new values for the covariates X at which predictions are to be made.
#' @param newA a (n-by-L) matrix of new values for the treatment A at which predictions are to be made.
#' @param L when \code{newA=NULL}, a value specifying the length of the grid of A at which predictions are to be made.
#' @param type the type of prediction required; the default "response" is on the scale of the response variable; the alternative "link" is on the scale of the linear predictors.
#' @param maximize the default is \code{TRUE}, assuming a larger value of the outcome is better; if \code{FALSE}, a smaller value is assumed to be prefered.
#'
#' @return
#' \item{pred.new}{a (n-by-L) matrix of predicted values; each column represents a treatment option.}
#' \item{trt.rule}{a (n-by-1) vector of suggested treatment assignments}
#'
#'
#' @author Park, Petkova, Tarpey, Ogden
#' @seealso \code{simsl},\code{fit.simsl}
#' @export
#'
pred.simsl  <-  function(simsl.obj, newx, newA =NULL, L=30, type = "response", maximize=TRUE)
{
  #if(!inherits(simsl.obj, "simsl"))   # checks input
  #  stop("obj must be of class `simsl'")
  if(ncol(newx) != simsl.obj$p)
    stop("newx needs to be of p columns ")

  if(is.null(simsl.obj$X.scale)){
    if(is.null(simsl.obj$X.scale)){
      newx.scaled <- scale(newx, center = rep(0,simsl.obj$p), scale = rep(1,simsl.obj$p))
    }else{
      newx.scaled <- scale(newx, center = simsl.obj$X.center, scale = rep(1,simsl.obj$p))
    }
  }else{
    newx.scaled <- scale(newx, center = simsl.obj$X.center, scale = simsl.obj$X.scale)
  }
  single.index  <- newx.scaled %*% simsl.obj$beta.coef

  # compute treatment (A)-specific predicted outcomes
  if(is.null(newA)){
    A.grid <- seq(from =simsl.obj$A.range[1], to =simsl.obj$A.range[2], length.out =L)
    newA <-  matrix(A.grid, nrow(newx), L, byrow = TRUE)
  }else{
    newA <- as.matrix(newA)
    L <- ncol(newA)
  }

  pred.new <- matrix(0, nrow(newx), L)
  for(a in 1:L){
    newd <- data.frame(A= newA[,a], single.index=single.index)
    pred.new[ ,a] <- predict.gam(simsl.obj$g.fit, newd, type =type)
    rm(newd)
  }

  # compute optimal treatment assignment
  if(maximize){
    opt.trt.index <- apply(pred.new, 1, which.max)
  }else{
    opt.trt.index <- apply(pred.new, 1, which.min)
  }

  trt.rule <- rep(NA, nrow(pred.new))
  for(i in 1:nrow(pred.new)){
    trt.rule[i] <- newA[i, opt.trt.index[i]]
  }

  return(list(trt.rule = trt.rule, pred.new = pred.new))
}



######################################################################
## END OF THE FILE
######################################################################
