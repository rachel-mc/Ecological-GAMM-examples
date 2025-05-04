## Required packages:

library(FSA)
library(hnp)
library(gamlss)

################################################################################
################### "hnp" HELPER FUNCTIONS FOR "gamlss" ########################
################################################################################

## Up to now, "hnp" was unable to assess the adequacy of most GA(M)Ms fitted in
## "gamlss". Three helper functions are required for this. These are dfun(),
## sfun(), and ffun() functions, which are described below. Note that the dfun()
## function remains the same regardless of the GA(M)M being assessed.

## dfun(): a helper function to extract the "gamlss" quantile residuals.

dfun <- function(obj) resid(obj) 

## sfun(): a helper function to perform empirically-derived simulations to determine
## how the residuals should behave under the distributional assumptions upon 
## which the GA(M)M was fitted. This helper function is specific to the family
## that was used.

sfun <- function(n, obj) { 
  ## FAMILY-SPECIFIC CODE: SEE EXAMPLES BELOW 
} 

## ffun(): a helper function to refit the model 99 times (default value).
## The contents of this function are also specific to the model being assessed.
## The right-hand side of the fitted model's formula is typed next to "resp ~"

## To create a specific ffun() helper function, the model, family, and data used
## must be defined first using these object names.
## See examples below.


# hnp() Syntax ------------------------------------------------------------


## The hnp() function can be used to assess the adequacy of a GA(M)M by
## relying on the three previously defined helper functions. The following
## syntax can be used for all GA(M)Ms. Here, a single hnp iteration will be
## performed based on 99 simulations (default). The argument how.many.out = TRUE
## indicates that the percentage of residuals found outside the simulated 
## envelope of the half-normal plot will be calculated and displayed. 
## The argument plot = TRUE will produce the corresponding half-normal plot.
## The argument paint = TRUE will print the residuals found outside the
## simulated envelope in red.

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

## This code will give an error unless the model, dfun, sfun, and ffun helper
## functions are previously stored. 

## Using the specific dfun(), sfun(), and ffun() helper functions that relate 
## to a given GA(M)M, the code and syntax below allow to obtain a mode and mean
## percentage of residuals found outside the simulated envelope based on the
## number of iterations (n) specified by the user. The set.seed() function is 
## used solely for result reproducibility. Here, n = 10 "hnp" iterations are
## performed. This code generates a list of results saved as an object called
## "hnp_obj", and is applicable to any model for which the dfun(), sfun(), and 
## ffun() helper functions have been previously defined.

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

## Calculate the percentage of residuals found outside the envelope for each 
## "hnp" iteration listed in hnp_obj and store it in the hnp_summary object

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

## Obtain the mode (asymptotic value) from a density curve

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
  }

round(return_max(hnp_summary), 2)

## Descriptive statistics, including the mean, using the Summarize() function of 
## the "FSA" package

Summarize(hnp_summary)


# Gaussian distribution and identity link ---------------------------------


## Below, a fictive GAMM is provided as an example with family = NO.

## Simulate data - same as for the "mgcv" Gaussian case

set.seed(2025)

X1 <- rnorm(n = 300)
X2 <- rnorm(n = 300)
X3 <- factor(sample(1:5, size = 300, replace = TRUE)) 

Y <- X1^3 + 5 * X2 + rnorm(nlevels(X3))[as.numeric(X3)] + rnorm(n = 300, mean = 1, sd = 2) 

CONTINUOUS <- data.frame(Y = Y,
                         X1 = X1,
                         X2 = X2,
                         X3 = X3)

GAMM_GAUSSIAN <- gamlss(Y ~ pb(X1) + X2 + random(X3), 
                        family = NO, ## same as family = NO(mu.link = "identity")
                        data = CONTINUOUS) 

## Indicate the model being assessed, the family, and the data that were used to
## build the investigated GA(M)M. 

model <- GAMM_GAUSSIAN
family = NO()
data = CONTINUOUS

dfun <- function(obj) resid(obj) 

sfun <- function(n, obj) { 
  mu <- obj$mu.fv 
  sig <- obj$sigma.fv 
  rNO(n, mu = mu, sigma = sig) 
  }

ffun <- function(resp) gamlss(resp ~ pb(X1) + X2 + random(X3),
                              family = family,
                              data = data) 

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

set.seed(2025)

n <- 10

hnp_obj <- list()
for (i in 1:n) { ## run twice
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


# Gamma distribution and log link -----------------------------------------


set.seed(2025)

# Use the same X1, X2, and X3 as defined above

mu <- exp(1 + 1.6 * X1 - 0.8 * X2 + 2 * as.numeric(X3))
Y <- rgamma(n = 300, shape = 2, rate = 2 / mu)

POSITIVE_CONTINUOUS <- data.frame(Y = Y,
                                  X1 = X1,
                                  X2 = X2,
                                  X3 = X3)

GAMM_GAMMA <- gamlss(Y ~ pb(X1) + X2 + random(X3), 
                     family = GA, ## same as family = GA(mu.link = "log")
                     data = POSITIVE_CONTINUOUS) 

model <- GAMM_GAMMA
family = GA
data = POSITIVE_CONTINUOUS

dfun <- function(obj) resid(obj) 

sfun <- function(n, obj) { 
  mu <- obj$mu.fv 
  sig <- obj$sigma.fv 
  rGA(n, mu = mu, sigma = sig) 
  }

ffun <- function(resp) gamlss(resp ~ pb(X1) + X2 + random(X3),
                              family = family,
                              data = data)

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)

set.seed(2025)

n <- 10

hnp_obj <- list()
for (i in 1:n) { ## This takes a few minutes to run.
  hnp_obj[[i]] <- hnp(model, 
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
  }

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)


# Poisson distribution (log link) and its extensions ----------------------


## A specific function was created to produce half normal plots for
## log-linear models that analyse counts. For instance,
## GA(M)Ms fitted with PO (Poisson), NBI (NB2), NBII (NB1), DPO (Double Poisson),
## Generalised Poisson (GPO), and Poisson inverse-Gaussian (PIG), can all be
## assessed for adequacy using the hnp_gamlss_count() function described below.
## Any 3- and even 4-parameter distributions that are available in "gamlss" can
## also be used here.

hnp_gamlss_count <- function(model, 
                             d_fun, 
                             ...) {
  
  fam <- model$family[1]
  random_generator <- get(paste0("r", fam))
  params <- model$parameters
  
  n <- length(model$y)
  
  mu_hat <- predict(model, type = "response")
  if("sigma" %in% params) sigma_hat <- predict(model, what = "sigma", type = "response")
  if("nu" %in% params) nu_hat <- predict(model, what = "nu", type = "response")
  if("tau" %in% params) tau_hat <- predict(model, what = "tau", type = "response")
  
  if(length(params) == 1) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat)
      return(y_new)
    }
  } else if(length(params) == 2) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat, sigma = sigma_hat)
      return(y_new)
    }
  } else if(length(params) == 3) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat, sigma = sigma_hat, nu = nu_hat)
      return(y_new)
    }
  } else if(length(parms) == 4) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat, sigma = sigma_hat, nu = nu_hat, tau = tau_hat)
      return(y_new)
    }
  }
  
  if (missing(d_fun)) d_fun <- function(obj) resid(obj) # default gamlss residuals
  
  full_data <- eval(model$call$data)
  
  f_fun <- function(y_new) {
    full_data$y <- y_new
    newfit <- update(model, formula = y ~ ., data = full_data)
    return(newfit)
  }
  
  hnp_results <- hnp(model,
                     newclass = TRUE,
                     diagfun = d_fun,
                     simfun = s_fun,
                     fitfun = f_fun,
                     ...)
  
  return(invisible(hnp_results))
}


## This hnp_gamlss_count() function can then be directly applied, as follows.

set.seed(2025)

lambda <- exp(1 + 2* X1 - 1.5 * X2 + 0.38 * as.numeric(X3))  
# Using the same X1, X2, and X3 as above
Y <- rpois(n = 300, lambda = lambda)

COUNT <- data.frame(Y = Y,
                    X1 = X1,
                    X2 = X2,
                    X3 = X3)

GAMM_POISSON <- gamlss(Y ~ pb(X1) + X2 + random(X3), 
                       family = PO, ## same as family = PO(mu.link = "log")
                       data = COUNT)

model <- GAMM_POISSON ## no need to define the family and data that were used

hnp_gamlss_count(model,
                 how.many.out = TRUE,
                 plot.sim = TRUE, 
                 paint = TRUE)


## To perform a number of "hnp" iterations (n, n > 1) to get a mode and mean
## percentage of residuals that lie outside the simulated envelope, the
## following code can be used.

set.seed(2025)

n <- 10

hnp_obj <- list()
for (i in 1:10) { ## This takes a few minutes to run.
  hnp_obj[[i]] <- hnp_gamlss_count(model, 
                                how.many.out = TRUE,
                                plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max(hnp_summary)

Summarize(hnp_summary)


# Binomial distribution and its extensions --------------------------------


## Binomial distribution (BI) with logit link

set.seed(2025)

TOTAL <- sample(5:20, n, replace = TRUE)
eta <- 0.4 + X1^3 + 2.1 * X2 - 0.6 * as.numeric(X3) # systematic component using X1, X2, and X3 as defined above
p <- 1 / (1 + exp(-eta)) # probabilities for the logit link
YES <- rbinom(n, size = TOTAL, prob = p)

DISCRETE_PROP <- data.frame(YES = YES,
                            TOTAL = TOTAL,
                            X1 = X1,
                            X2 = X2,
                            X3 = X3)

GAMM_BI <- gamlss(cbind(YES, TOTAL - YES) ~ pb(X1) + X2 + random(X3), 
                  family = BI, ## same as family = BI(mu.link = "logit")
                  data = DISCRETE_PROP)

model <- GAMM_BI
family <- BI
data <- DISCRETE_PROP

dfun <- function(obj) resid(obj)

sfun <- function(n, obj) {
  mu <- obj$mu.fv
  bd <- obj$bd
  rBI(n, bd = bd, mu = mu)
  }

ffun <- function(new_response) { 
  gamlss(cbind(new_response, TOTAL - new_response) ~ pb(X1) + X2 + random(X3),
           family = family,
           data = data)
}

hnp(model,
    newclass = TRUE,
    diagfun = dfun,
    simfun = sfun,
    fitfun = ffun,
    how.many.out = TRUE,
    paint = TRUE)

set.seed(2025)

n <- 10

hnp_obj <- list()
for (i in 1:10) { ## This takes a few minutes to run.
  hnp_obj[[i]] <- hnp_gamlss_count(model, 
                                   how.many.out = TRUE,
                                   plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max(hnp_summary)

Summarize(hnp_summary)

## Beta-binomial (BB), zero-inflated binomial (ZIBI), and zero-inflated
## beta-binomial (ZIBB) distributions all fitted with the logit link (default).

GAMM_BB <- gamlss(cbind(YES, TOTAL - YES) ~ pb(X1) + X2 + random(X3), 
                  family = BB,
                  data = DISCRETE_PROP)

GAMM_ZIBI <- gamlss(cbind(YES, TOTAL - YES) ~ pb(X1) + X2 + random(X3), 
                    family = ZIBI,
                    data = DISCRETE_PROP)

GAMM_ZIBB <- gamlss(cbind(YES, TOTAL - YES) ~ pb(X1) + X2 + random(X3), 
                    family = ZIBB,
                    data = DISCRETE_PROP)

## These BI extensions can all be handled directly by the hnp() function. However, 
## modelling low to moderate sample sizes of discrete proportions is prone 
## to analytical issues from data resampling. The tryCatch() function, a while loop,
## and additional code are used below to ensure that a result will be obtained.
## Here is an example for the adequacy assessment of GAMM_BB using "hnp".
## In this case, the arguments sim = 19 and conf = 1 are specified to reduce
## run time (see description below).

model <- GAMM_BB

run_complete <- FALSE
while (run_complete == F) {
  tryCatch( {
    hnp(model, how.many.out = TRUE, paint = TRUE,
        sim = 19, conf = 1)
    run_complete <- TRUE  
  },
  error = function(e) message("Trying the \"hnp\" iteration again...")
  )
}

## ZIBI and ZIBB modelling issue risks are increased due to the zero-
## inflated term present. As such, accommodating for errors may still lead to an
## excessively long computational process. To produce a result
## in a reasonable amount of time, the default number of simulations is reduced
## from 99 to 19 and the default level of confidence is modified from
## 95% (0.95) to 100% (1) (see Moral et al. 2017).

model <- GAMM_ZIBB

run_complete <- FALSE
while (run_complete == F) { 
  tryCatch( {
    hnp(model, how.many.out = TRUE, paint = TRUE,
        sim = 19, conf = 1)
    run_complete <- TRUE  
  },
  error = function(e) message("Trying the \"hnp\" iteration again...")
  )
}


## Performing two or more "hnp" iterations (here n = 3), using GAMM_BB as an
## example.

model <- GAMM_BB

set.seed(2025)

n <- 3

hnp_obj <- list()
run_complete <- rep(FALSE, n)
for (i in 1:n) {
  while (run_complete[i] == F) { ## This takes several minutes to run.
    tryCatch( {
      hnp_obj[[i]] <- hnp(model, how.many.out = TRUE, plot.sim = FALSE,
                          sim = 19, conf = 1) 
      run_complete[i] <- TRUE 
    }, 
    error = function(e) message("Trying \"hnp\" iteration ", i, " again...")
    )
  }
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100)

round(return_max(hnp_summary), 2)

Summarize(hnp_summary) 


# References --------------------------------------------------------------


## Moral RA, Hinde J, Demétrio CGB (2017) Half-normal plots and overdispersed
##    models in R: the hnp package. J Stat Softw 81:1–23. https://doi.org/10.18637/jss.v081.i10
