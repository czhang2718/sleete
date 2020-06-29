# preliminary functions

# truncation
trunc <- function(u, lo=-Inf, hi=Inf) pmax(lo,pmin(hi,u))

# logit = log-odds
logit <- function(p) log(p/(1-p))

# Agresti definition of h
h0 <- function(y1, y0) as.numeric(y1>y0)-as.numeric(y1<y0)

# Mann-Whitney definition of h
h1 <- function(y1, y0) as.numeric(y1>y0)+0.5*as.numeric(y1==y0)

# Kaplan-Meier estimate of a restricted (by tau) distribution
km.tau <- function(x, delta, tau=Inf) {
  dat <- Surv(x, delta)
  km.est <- survfit(dat~1)
  t.all <- c(km.est$time, Inf)
  s.all <- c(1, km.est$surv, 0)
  p.all <- s.all[-length(s.all)]-s.all[-1]
  t.is.small <- (t.all<tau)
  t.tau <- c(t.all[t.is.small], tau)
  p.tau <- p.all[t.is.small]
  pr.lt.tau <- sum(p.tau)
  p.tau <- c(p.tau, 1-pr.lt.tau)
  list(t=t.tau, p=p.tau)
}

# methods

# method for difference between two means or proportions
# point estimate
pt.est.mean.diff <- function(y, t) mean(y[t>0.5])-mean(y[t<0.5])
# influence function estimated from I, applied to J
inf.fct.mean.diff <- function(y, t, I=1:length(t), J=I, pi=NULL) {
  if (is.null(pi)) pi <- mean(t[I])
  mu1 <- mean(y[I][t[I]>0.5]); mu0 <- mean(y[I][t[I]<0.5])
  (t[J]*(y[J]-mu1)/pi)-((1-t[J])*(y[J]-mu0)/(1-pi))
}
# the method
mean.diff <- list(pt.est=pt.est.mean.diff, inf.fct.avail=TRUE, inf.fct=inf.fct.mean.diff)

# method for log-ratio of two means or proportions
# point estimate
pt.est.log.ratio <- function(y, t) log(mean(y[t>0.5])/mean(y[t<0.5]))
# influence function estimated from I, applied to J
inf.fct.log.ratio <- function(y, t, I=1:length(t), J=I, pi=NULL) {
  if (is.null(pi)) pi <- mean(t[I])
  mu1 <- mean(y[I][t[I]>0.5]); mu0 <- mean(y[I][t[I]<0.5])
  (t[J]*(y[J]-mu1)/(pi*mu1))-((1-t[J])*(y[J]-mu0)/((1-pi)*mu0))
}
# the method
log.ratio <- list(pt.est=pt.est.log.ratio, inf.fct.avail=TRUE, inf.fct=inf.fct.log.ratio)

# method for log-odds-ratio of two proportions
# point estimate
pt.est.log.or <- function(y, t) logit(mean(y[t>0.5]))-logit(mean(y[t<0.5]))
# influence function estimated from I, applied to J
inf.fct.log.or <- function(y, t, I=1:length(t), J=I, pi=NULL) {
  if (is.null(pi)) pi <- mean(t[I])
  p1 <- mean(y[I][t[I]>0.5]); p0 <- mean(y[I][t[I]<0.5])
  (t[J]*(y[J]-p1)/(pi*p1*(1-p1)))-((1-t[J])*(y[J]-p0)/((1-pi)*p0*(1-p0)))
}
# the method
log.odds.ratio <- list(pt.est=pt.est.log.or, inf.fct.avail=TRUE, inf.fct=inf.fct.log.or)

# method for Wilcoxon-Mann-Whitney effect based on an arbitrary h (default = h0)
# point estimate
pt.est.wmw <- function(y, t, h=h0) mean(outer(y[t>0.5], y[t<0.5], FUN=h))
# influence function estimated from I, applied to J
inf.fct.wmw <- function(y, t, I=1:length(t), J=I, pi=NULL, h=h0) {
  if (is.null(pi)) pi <- mean(t[I])
  theta <- pt.est.wmw(y[I],t[I],h=h)
  m <- length(J); inf <- numeric(m)
  for (k in 1:m) {
    if (t[J[k]]>0.5) {
      inf[k] <- (mean(h(y[J[k]],y[I]))-theta)/pi
    } else {
      inf[k] <- (mean(h(y[I],y[J[k]]))-theta)/(1-pi)
    }
  }
  inf
}
# the method
wmw <- list(pt.est=pt.est.wmw, inf.fct.avail=TRUE, inf.fct=inf.fct.wmw)

# method for Wilcoxon-Mann-Whitney effect basedon an arbitrary h (default = h0)
# estimated from randomly right-censored data restricted by tau
# point estimate
pt.est.wmw.cens <- function(y, t, tau=Inf, h=h0) {
  x <- y[,1]; delta <- y[,2]
  t.eq.1 <- (t>0.5)
  km.rst.1 <- km.tau(x[t.eq.1], delta[t.eq.1], tau=tau)
  t1 <- km.rst.1$t; p1 <- km.rst.1$p
  km.rst.0 <- km.tau(x[!t.eq.1], delta[!t.eq.1], tau=tau)
  t0 <- km.rst.0$t; p0 <- km.rst.0$p
  H <- outer(t1, t0, FUN=h)
  as.vector(t(p1)%*%H%*%p0)
}
# the method
wmw.cens <- list(pt.est=pt.est.wmw.cens, inf.fct.avail=FALSE)

# method for difference between two survival probabilities at time tau
# estimated from randomly right-censored data using the Kaplan-Meier approach
# point estimate
pt.est.surv.diff <- function(y, t, tau=Inf) {
  x <- y[,1]; delta <- y[,2]
  t.eq.1 <- (t>0.5)
  km.rst.1 <- km.tau(x[t.eq.1], delta[t.eq.1])
  s1 <- 1-sum(km.rst.1$p[km.rst.1$t<=tau])
  km.rst.0 <- km.tau(x[!t.eq.1], delta[!t.eq.1])
  s0 <- 1-sum(km.rst.0$p[km.rst.0$t<=tau])
  s1-s0
}
# the method
surv.diff <- list(pt.est=pt.est.surv.diff, inf.fct.avail=FALSE)

# method for difference between two mean restricted (by tau) survival times
# estimated from randomly right-censored data using the Kaplan-Meier approach
# point estimate
pt.est.mrst.diff <- function(y, t, tau=Inf) {
  x <- y[,1]; delta <- y[,2]
  t.eq.1 <- (t>0.5)
  km.rst.1 <- km.tau(x[t.eq.1], delta[t.eq.1], tau=tau)
  mu1 <- sum(km.rst.1$p*km.rst.1$t)
  km.rst.0 <- km.tau(x[!t.eq.1], delta[!t.eq.1], tau=tau)
  mu0 <- sum(km.rst.0$p*km.rst.0$t)
  mu1-mu0
}
# the method
mrst.diff <- list(pt.est=pt.est.mrst.diff, inf.fct.avail=FALSE)
rmst.diff <- mrst.diff

# method for log-hazard-ratio in Cox model
# point estimate
pt.est.haz.ratio <- function(y, t) {
  x <- y[,1]; delta <- y[,2]
  dat <- Surv(x, delta)
  mod <- coxph(dat~t)
  mod$coeff
}
# the method
log.haz.ratio <- list(pt.est=pt.est.haz.ratio, inf.fct.avail=FALSE)

# empirical influence function for a generic point estimator
# estimated from subjects in set I, then applied to subjects in set J
emp.inf.fct <- function(est, y, t, I=1:length(t), J=I) {
  y.is.matrix <- is.matrix(y)
  if (y.is.matrix) est0 <- est(y[I,], t[I])
  else est0 <- est(y[I], t[I])
  m <- length(I); n <- length(J); jack.est <- numeric(n)
  for (j in 1:n) {
    Ij <- c(I, J[j])
    if (y.is.matrix) jack.est[j] <- est(y[Ij,], t[Ij])
    else jack.est[j] <- est(y[Ij], t[Ij])
  }
  (m+1)*(jack.est-est0)
}

#'Estimates treatment effect in randomized clinical trials using a SuperLearner
#'
#'The sleete function uses a super learner to minimize the variance of an
#'augmented estimator of a specified treatment effect measure in a randomized
#'clinical trial. It returns a matrix of point estimates and standard errors for
#'the super learner as well as individual algorithms in the super learner
#'library, with or without sample splitting.
#'
#'Currently, there are eight built-in methods available for \code{method}. Four
#'of them are for fully observed univariate outcomes: \code{mean.diff} for the
#'difference between two means or proportions, \code{log.ratio} for the
#'log-ratio of two means or proportions, \code{log.odds.ratio} for the
#'log-odds-ratio of two proportions, and \code{wmw} for the
#'Wilcoxon-Mann-Whitney (WMW) effect (Zhang et al., 2019), the default version
#'of which is also known as the win-lose probability difference. There are four
#'other methods for right-censored survival outcomes: \code{wmw.cens} for the
#'WMW effect for restricted survival times, \code{surv.diff} for the difference
#'between two survival probabilities, \code{mrst.diff} (or \code{rmst.diff}) for
#'the difference in mean restricted survival time, and \code{log.haz.ratio} for
#'the log-hazard-ratio. The methods for right-censored survival outcomes are
#'implemented without an analytical influence function (i.e.,
#'\code{inf.fct.avail=FALSE}). Users can define their own methods under the same
#'guidelines. For illustration, the current definitions of the
#'\code{log.odds.ratio} and \code{wmw} methods are provided below as examples.
#'
#'@param y Outcome data represented as a vector (for a univariate outcome) or a
#'  matrix (for a right-censored survival outcome or multiple outcomes to be
#'  analyzed together). For a right-censored survival outcome, y is a matrix
#'  with two columns: observed time followed by event type (1 failure; 0
#'  censoring).
#'@param t A vector of 1s and 0s representing treatment assignment. The values 1
#'  and 0 represent the experimental and control treatments, respectively. The
#'  length of t should be equal to the number of subjects.
#'@param X A matrix of baseline covariates that may be related to y in one or
#'  both treatment groups. The number of rows in X should be equal to the number
#'  of subjects. The number of columns in X is the number of covariates. There
#'  is no need to have a column of 1s in X.
#'@param pi The probability of receiving the experimental treatment, usually
#'  known in a randomized clinical trial. If missing, will be replaced by the
#'  proportion of study subjects who were assigned to the experimental
#'  treatment.
#'@param bounds Known lower and upper bounds, if any, for the treatment effect
#'  measure to be estimated. For example, if the effect measure is a difference
#'  between two probabilities, the natural bounds are c(-1,1).
#'@param method A list of two mandatory components and one optional component
#'  specifying the (unadjusted) method for estimating the treatment effect of
#'  interest. The two mandatory components are \code{pt.est}, a function for
#'  obtaining a point estimate, and \code{inf.fct.avail}, an indicator for the
#'  availability of a function to compute the influence function of the point
#'  estimator analytically. If the value of \code{inf.fct.avail} is TRUE, one
#'  has to also supply a function named \code{inf.fct} to compute the influence
#'  function of the point estimator analytically. If the value of
#'  \code{inf.fct.avail} is FALSE, the function \code{inf.fct} is not needed and
#'  the empirical influence function (Zhang et al., 2020) will be computed.
#'  Following the \code{method} argument, there are optional arguments
#'  represented by dots (\code{...}). If specified, such optional arguments will
#'  be fed into the specified method. For instance, the \code{wmw} and
#'  \code{wmw.cens} methods involve a kernel function h. The default for h
#'  (named h0) and an illustrative alternative h1 are provided below as
#'  examples. The \code{wmw.cens}, \code{surv.diff} and \code{mrst.diff} methods
#'  require specifying a time point tau, which has no default value and must be
#'  supplied by the user.
#'@param SL.library A character vector of SuperLearner wrapper functions for the
#'  prediction algorithms that comprise the super learner library. A full list
#'  of wrapper functions included in the SuperLearner package can be found with
#'  \code{listWrappers()}.
#'@param cv The number of folds in the cross-validation for the super learner.
#'@param cf The number of folds in the sample splitting or cross-fitting
#'  procedure
#'@seealso See \code{\link[SuperLearner]{SuperLearner}} for details on
#'  \code{method}, \code{SL.library}, \code{cv}, \code{cf}, and \code{family}.
#'@return A matrix with two columns: point estimates of the treatment effect of
#'  interest and their standard errors. The number of rows is 2K+3, where K is
#'  the length of \code{SL.library}. The first row is for the unadjusted
#'  estimate as specified in the method argument. The next K+1 rows are for
#'  augmented estimates based on the individual algorithms in the super learner
#'  library (in the original order) followed by the super learner itself, all
#'  without sample splitting. The next K+1 rows are for augmented estimates
#'  based on the same set of algorithms (in the same order) with sample
#'  splitting. The standard error for the unadjusted estimate is based on the
#'  (analytical or empirical) influence function. The standard errors for the
#'  augmented estimates are cross-validated in the sample splitting procedure.
#'  Thus, the two sub-sets of augmented estimates (with and without sample
#'  splitting) have the same set of cross-validated standard errors.
#' @examples
#' # analysis of colon cancer data in the survival package
#' library(survival)
#' data(colon)
#' dim(colon); names(colon)
#' colon.data <- na.omit(subset(colon, subset=((etype==2)&(rx!="Lev")),
#' select <- c(rx, time, status, sex, age, obstruct, perfor,
#'   adhere, nodes, node4, surg, differ, extent)))
#' dim(colon.data)
#' attach(colon.data)
#' t = as.numeric(rx=="Lev+5FU")
#' y = cbind(time, status)
#' X = cbind(sex,age,obstruct,perfor,adhere,nodes,node4,surg,differ,extent)
#' detach()
#' pi = 0.5; tau = 5*365
#' sleete(y, t, X, pi=pi, method=surv.diff, bounds=c(-1,1), tau=tau)
#' sleete(y, t, X, pi=pi, method=mrst.diff, tau=tau)
#' sleete(y, t, X, pi=pi, method=wmw.cens, bounds=c(-1,1), tau=tau)
#' # the log-odds-ratio method
#' # logit = log-odds
#' logit = function(p) log(p/(1-p))
#' # point estimate
#' pt.est.log.or = function(y, t) logit(mean(y[t>0.5]))-logit(mean(y[t<0.5]))
#' # influence function estimated from subjects in set I
#' # then applied to subjects in set J
#' inf.fct.log.or = function(y, t, I=1:length(t), J=I, pi=NULL) {
#'   if (is.null(pi)) pi = mean(t[I])
#'   p1 = mean(y[I][t[I]>0.5]); p0 = mean(y[I][t[I]<0.5])
#'   (t[J]*(y[J]-p1)/(pi*p1*(1-p1)))-((1-t[J])*(y[J]-p0)/((1-pi)*p0*(1-p0)))
#' }
#' log.odds.ratio = list(pt.est=pt.est.log.or, inf.fct.avail=TRUE, inf.fct=inf.fct.log.or)
#'
#' # the wmw method with an arbitrary h (default = h0)
#' # Agresti definition of h
#' h0 = function(y1, y0) as.numeric(y1>y0)-as.numeric(y1<y0)
#' # Mann-Whitney definition of h
#' h1 = function(y1, y0) as.numeric(y1>y0)+0.5*as.numeric(y1==y0)
#' # point estimate
#' pt.est.wmw = function(y, t, h=h0) mean(outer(y[t>0.5], y[t<0.5], FUN=h))
#' # influence function estimated from subjects in set I
#' # then applied to subjects in set J
#' inf.fct.wmw = function(y, t, I=1:length(t), J=I, pi=NULL, h=h0) {
#'   if (is.null(pi)) pi = mean(t[I])
#'   theta = pt.est.wmw(y[I],t[I],h=h)
#'   m = length(J); inf = numeric(m)
#'   for (k in 1:m) {
#'     if (t[J[k]]>0.5) {
#'       inf[k] = (mean(h(y[J[k]],y[I]))-theta)/pi
#'     } else {
#'       inf[k] = (mean(h(y[I],y[J[k]]))-theta)/(1-pi)
#'     }
#'   }
#'   inf
#' }
#' wmw = list(pt.est=pt.est.wmw, inf.fct.avail=TRUE, inf.fct=inf.fct.wmw)
#'@export

sleete <- function(y, t, X, pi=mean(t), bounds=c(-Inf, Inf), method=mean.diff, ...,
        SL.library=c("SL.mean", "SL.glm", "SL.gam", "SL.rpart", "SL.randomForest"),
        cv=5, cf=5) {
  n <- length(t); K <- length(SL.library)
  est <- function(yy,tt) method$pt.est(yy,tt,...)
  if (method$inf.fct.avail) {
    inf.fct <- function(yy,tt,I=1:length(tt),J=I) method$inf.fct(yy,tt,I=I,J=J,pi=pi,...)
  } else {
    inf.fct <- function(yy,tt,I=1:length(tt),J=I) emp.inf.fct(est,yy,tt,I=I,J=J)
  }
  est.0 <- est(y,t)
  inf <- inf.fct(y,t)
  se.0 <- sd(inf)/sqrt(n)
  rsp <- inf/(t-pi); wt <- (t-pi)^2; X.df <- data.frame(X)
  SL.rst <- SuperLearner::SuperLearner(Y=rsp, X=X.df, newX=X.df,
    family=gaussian(), SL.library=SL.library, cvControl=list(V=cv), obsWeights=wt)
  b.hat <- cbind(SL.rst$library.predict, SL.rst$SL.predict)
  aug.mat <- (t-pi)*b.hat
  est.aug <- est.0-colMeans(aug.mat)
  se.aug <- sqrt(diag(var(inf-aug.mat))/n)
  inf.ss <- numeric(n)
  b.hat.ss <- matrix(0,n,K+1)
  fold <- sample(1:cf, n, replace=TRUE)
  for (v in 1:cf) {
    val <- (fold==v); trn <- !val
    I <- (1:n)[trn]; J <- (1:n)[val]
    inf.ss[val] = inf.fct(y,t,I=I,J=J)
    inf.trn <- inf.fct(y,t,I=I)
    rsp.trn <- inf.trn/(t[trn]-pi)
    SL.rst.trn <- SuperLearner::SuperLearner(Y=rsp.trn, X=X.df[trn,],
      newX=X.df[val,], family=gaussian(), SL.library=SL.library,
      cvControl=list(V=cv), obsWeights=wt[trn])
    b.hat.ss[val,] <- cbind(SL.rst.trn$library.predict, SL.rst.trn$SL.predict)
  }
  aug.mat.ss <- (t-pi)*b.hat.ss
  est.ss <- est.0-colMeans(aug.mat.ss)
  se.ss <- sqrt(diag(var(inf.ss-aug.mat.ss))/n)
  pe <- trunc(c(est.0, est.aug, est.ss), lo=bounds[1], hi=bounds[2])
  se <- c(se.0, se.ss, se.ss)
  rst <- cbind(pe, se)
  colnames(rst) <- c("Pt.Est", "Std.Err")
  SL.names <- c(SL.library, "SL")
  rownames(rst) <- c("Unadjusted", paste(SL.names, "w/o cross-fitting"),
        paste(SL.names, "w/ cross-fitting"))
  rst
}
