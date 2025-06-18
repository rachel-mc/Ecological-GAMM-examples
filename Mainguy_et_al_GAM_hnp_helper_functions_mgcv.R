## Required packages:

library(FSA)
library(hnp)
library(mgcv)

################################################################################
##################### "hnp" helper functions for "mgcv" ########################
################################################################################

## Up to now, "hnp" was unable to assess the adequacy of GAMs (and GAMMs) fitted 
## in "mgcv". Three helper functions are now available to allow such assessment,
## namely the dfun(), sfun(), and ffun() functions, which are all described below.
## Note that the dfun() function remains the same at all times.

## dfun(): a helper function to extract the deviance residuals of a model.

dfun <- function(obj) resid(obj, type = "deviance")

## sfun(): a helper function to perform empirically-derived simulations to 
## determine how the residuals should behave under the distributional assumptions  
## of the fitted GAM. This helper function is specific to the family being used.

## sfun() takes the following form:

sfun <- function(n, obj) {
  y <- ## specific to the considered family - see examples below
  return(y)
}

## ffun(): a helper function to refit the model 99 times (default value) for 
## each new simulated response. This function is also specific to the model being
## assessed. 

## Sample construction of the ffun() helper function:

## First define the fitted model, family, method, and data used.

model <- ## model name
family <- ## family used, such as poisson or nb (i.e., NB2)
method <- ## method used, such as "ML" or "REML"
data <- ## name of dataset

## Then, type the right-hand side of the fitted model formula next to 'resp ~ '

ffun <- function(resp) {
  gam(resp ~ ## indicate the right-hand side of the fitted model formula
      family = family,
      method = method,
      data = data)
}

# The hnp() syntax below can be used for all "mgcv" GAMs and GAMMs -------------

## The hnp() function can now be used to assess the adequacy of an "mgcv" object
## by relying on the three previously defined helper functions. Here, 
## a single "hnp" iteration will be performed based on 99 simulations (default).
## The argument how.many.out = TRUE will indicate the percentage of residuals
## found outside the simulated envelope of the half-normal plot. The corresponding
## half-normal plot will be shown with the argument plot = TRUE. The argument
## paint = TRUE will show the residuals found outside the envelope in red. The
## code below is just provided as an example and thus will not work if run.

hnp(model, ## same as the GAM or GAMM that was defined above
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

## The iterative process below can be used for all GAMs and GAMMs fitted in "mgcv".
## The specific dfun(), sfun(), and ffun() helper functions can be used to obtain 
## a mode and a mean percentage of residuals found outside the simulated
## envelope given the number of iterations (n) specified by the user.
## 10 "hnp" iterations are performed below. The set.seed() function is solely
## used for result reproducibility. This code generates a list of results
## saved as an object called "hnp_obj" and is applicable to any model (stored 
## as this alias) for which the dfun(), sfun(), and ffun() helper functions 
## have been previously correctly defined.

set.seed(2025)

n <- 10

hnp_obj <- list()
for (i in 1:n) {
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

## Calculate the percentage of residuals found outside the simulated envelope of 
## each "hnp" iteration and store it in an object referred to as hnp_summary

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

## Obtain the mode percentage (asymptotic value) from the associated density curve

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

## Descriptive statistics, including the mean percentage, using the Summarize() 
## function of the "FSA" package

Summarize(hnp_summary)

# GAUSSIAN DISTRIBUTION AND IDENTITY LINK  -------------------------------------

## Below, a fictive GAMM for which the deviance residuals are expected to
## approximate a normal distribution is provided as an example.

## Simulate continuous data using the normal distribution 

set.seed(2025)

X1 <- rnorm(n = 300)
X2 <- rnorm(n = 300)
X3 <- factor(sample(1:5, size = 300, replace = TRUE)) # factor with 5 levels for the random effect

Y <- X1^3 + 5 * X2 + rnorm(nlevels(X3))[as.numeric(X3)] + rnorm(n = 300, mean = 1, sd = 2) 
# assign random effect to each observation
# + stochastic error

CONTINUOUS <- data.frame(Y = Y,
                         X1 = X1,
                         X2 = X2,
                         X3 = X3)

## Fit a Gaussian GAMM

GAMM_GAUSSIAN <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                    family = gaussian(link = "identity"),
                    method = "REML",
                    data = CONTINUOUS) 

## Pre-define the information needed for ffun()

model <- GAMM_GAUSSIAN 
family <- gaussian(link = "identity")
method <- "REML"  
data <- CONTINUOUS

## Helper function to extract the deviance residuals

dfun <- function(obj) resid(obj, type = "deviance")

## Helper function to simulate how the residuals should behave given the
## distributional assumptions of the model being investigated. 
## The sfun() function is specific to the family used here (family = gaussian)

sfun <- function(n, obj) {
  y <- rnorm(nrow(data),
             mean = predict(model, type = "response"),
             sd = sqrt(model$sig2))
   return(y)
}

## Helper function to refit the model 99 times at each iteration 

ffun <- function(resp) {
  gam(resp ~ s(X1) + X2 + s(X3, bs = "re"),
      family = family,
      method = method,
      data = data)
}

## Asking for a single "hnp" iteration

set.seed(2025)

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

## Performing 10 "hnp" iterations to obtain a mode and a mean percentage of 
## residuals outside the simulated envelope

n <- 10

hnp_obj <- list()
for (i in 1:n) {
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

## Mode percentage of residuals outside the simulated envelope

round(return_max(hnp_summary), 2)

## Mean percentage of residuals outside the simulated envelope

Summarize(hnp_summary)


## GAMMA DISTRIBUTION AND LOG LINK ---------------------------------------------

## Simulate a new continuous response variable using the Gamma distribution 

set.seed(2025)

## Use the same X1, X2, and X3 as defined above, with Y requiring to be non-zero
## and only positive to allow the use of the gamma distribution and a log link.

shape <- 5
mu <- exp(1 + 1.6 * X1 - 0.8 * X2 + 2 * as.numeric(X3))  # log link
Y <- rgamma(n = 300, shape = shape, rate = shape / mu)

NON_ZERO_POSITIVE_CONTINUOUS <- data.frame(Y = Y,
                                           X1 = X1,
                                           X2 = X2,
                                           X3 = X3)

GAMM_GAMMA <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                  family = Gamma(link = "log"), 
                  method = "REML", 
                  data = NON_ZERO_POSITIVE_CONTINUOUS)

model <- GAMM_GAMMA
family <- Gamma(link = "log")
method <- "REML"  
data <- NON_ZERO_POSITIVE_CONTINUOUS

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n, obj) {
  y <- rgamma(nrow(data),
              shape = shape,
              rate = shape / predict(obj, type = "response"))
  return(y)
}

ffun <- function(resp) {
  gam(resp ~ s(X1) + X2 + s(X3, bs = "re"),
      family = family,
      method = method,
      data = data)
}

set.seed(2025)

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

n <- 10

hnp_obj <- list()
for (i in 1:n) {
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)


## POISSON DISTRIBUTION AND LOG LINK -------------------------------------------

## Simulate data from the type-II negative binomial (NB2) distribution to 
## illustrate that the Poisson distribution is inadequate for such overdispersed
## count data.

set.seed(2025)

mu <- exp(1 + 2* X1 - 1.5 * X2 + 0.38 * as.numeric(X3))  # log link
# same X1, X2, and X3 as above
Y <- rnbinom(n = 300, mu = mu, size = 3)

COUNT <- data.frame(Y = Y,
                    X1 = X1,
                    X2 = X2,
                    X3 = X3)

GAMM_POISSON <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                     family = poisson, ## same as family = poisson(link = "log")
                     method = "REML",
                     data = COUNT) 

model <- GAMM_POISSON 
family <- poisson
method <- "REML"  
data <- COUNT

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n,obj) {
  y <- rpois(nrow(data),
             lambda = predict(model, type = "response"))
  return(y)
  }

ffun <- function(resp) {
  gam(resp ~ s(X1) + X2 + s(X3, bs = "re"),
      family = family,
      method = method,
      data = data)
}

set.seed(2025)

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

n <- 10

hnp_obj <- list()
for (i in 1:n) {
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)


## TYPE-II NEGATIVE BINOMIAL (NB2) WITH LOG LINK -------------------------------

GAMM_NB2 <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                    family = nb, ## same as family = nb(theta = NULL, link = "log")
                    method = "REML",
                    data = COUNT) ## same COUNT data as in the Poisson example above

model <- GAMM_NB2 
family <- nb
method <- "REML"  
data <- COUNT

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n,obj) {
  y <- rnbinom(nrow(data),
               size = model$family$getTheta(TRUE),
               mu = predict(model, type = "response"))
  return(y)
}

ffun <- function(resp) {
  gam(resp ~ s(X1) + X2 + s(X3, bs = "re"), 
      family = family,
      method = method,
      data = data)
}

set.seed(2025)

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

n <- 10

hnp_obj <- list()
for (i in 1:n) {
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)

## The Poisson model was inadequate to model such overdispersed data. Based on 
## the mode and mean % above, this NB2 GAMM is adequate. 


## BINOMIAL DISTRIBUTION AND LOGIT LINK (analysis of discrete proportions) -----

## Generate new data

set.seed(2025)

X1 <- runif(n = 300) # random numbers between 0 and 1
X2 <- rnorm(n = 300)
X3 <- factor(sample(1:7, size = 300, replace = TRUE))
TOTAL <- sample(5:20, size = 300, replace = TRUE)

eta <- 0.4 + X1^3 + 1.5 * X2 - 0.6 * as.numeric(X3) # systematic component using
## X1, X2, and X3 defined above
p <- 1 / (1 + exp(-eta)) # probabilities
p_star <- pmin(pmax(p, 1e-6), 1 - 1e-6) # avoid probabilities of exactly 0 or 1

YES <- rbinom(n = 300, size = TOTAL, prob = p_star)

DISCRETE_PROP <- data.frame(YES = YES,
                            TOTAL = TOTAL,
                            X1 = X1,
                            X2 = X2,
                            X3 = X3)

GAMM_BINOMIAL_LOGIT <- gam(cbind(YES, TOTAL - YES) ~ s(X1) + X2 + s(X3, bs = "re"), 
                           family = binomial, ## same as binomial(link = "logit")
                           method = "REML",
                           data = DISCRETE_PROP) 

## For discrete proportions, size (i.e., the denominator) must also be provided
## within the simulation function. In this fictive dataset (DISCRETE_PROP), size
## makes up the TOTAL variable.

model <- GAMM_BINOMIAL_LOGIT
family <- binomial
method <- "REML"  
data <- DISCRETE_PROP

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n, obj) {
  p_hat <- predict(obj, type = "response")
  size <- obj$model[[1]][, 1] + obj$model[[1]][, 2] # successes and failures from cbind()
  y <- rbinom(n = n,
              size = size,
              prob = p_hat)
  return(y)
  }

ffun <- function(new_response) {
  data_copy <- data
  data_copy$new_response <- new_response
  
  gam(cbind(new_response, TOTAL - new_response) ~ s(X1) + X2 + s(X3, bs = "re"),
      family = family,
      method = method,
      data = data_copy)
  }

set.seed(2025)

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

n <- 10

hnp_obj <- list()
for (i in 1:n) {
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)
