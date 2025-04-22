## Required packages:

library(FSA)
library(hnp)
library(gamlss)


################################################################################
################### "HNP" HELPER FUNCTIONS FOR "GAMLSS" ########################
################################################################################

## As of now, "hnp" is unable to assess the adequacy of most GAM(M)s fitted with
## "gamlss". Three helper functions are required for this. These are the dfun(),
## sfun(), and ffun() functions, which are described below. Note that the dfun()
## function remains the same regardless of the GAM(M) being assessed.


## dfun(): a helper function to extract the quantile residuals

dfun <- function(obj) resid(obj) 


## sfun(): A helper function to perform empirically-derived simulations to determine
## how the residuals should behave under the distributional assumptions upon 
## which the GAM(M) was fitted. This helper function is specific to the family
## that was used.

sfun <- function(n, obj) { 
  ## FAMILY-SPECIFIC CODE: SEE EXAMPLES BELOW 
} 


## ffun(): A helper function to refit the model 99 times (default value)
## It is also specific to the model being assessed. A fictive GAM is presented
## below to illustrate the parametrisation of the ffun() helper function with
## Y = response variable, X1 = predictor, and family = NO (normal, i.e. Gaussian).

GAM_NO <- gamlss(Y ~ pb(X1), 
             family = NO(),
             data = CONTINUOUS)


## To define the ffun() helper function that is specific to GAM_NO, the model,
## family, and data used, are first defined.

model <- GAM_NO
family <- NO()
data <- CONTINUOUS


## The ffun() helper function is parametrised by indicating the right-hand part
## of the model next to "resp ~" which here is only "pb(X1)".

ffun <- function(resp) {
  gam(resp ~ pb(X1),
      family = family,
      method = method,
      data = data)
}


############# THE HNP() FUNCTION BELOW CAN BE USED WITH ALL GAM(M)S ############

## The hnp() function can now be used to assess the adequacy of a GAM by
## relying on the three previous helper functions that were just defined. Here, 
## a single hnp iteration will be performed based on 99 simulations (default).
## The argument how.many.out = TRUE indicates that the percentage of residuals 
## found outside the simulated envelope of the half-normal plot will be calculated
## and displayed. The argument plot = TRUE will produce the corresponding half-
## normal plot. Using the argument paint = TRUE will print in red the residuals
## found outside the simulated envelope.

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)


## The code below allows, using the specific dfun(), sfun(), ffun() helper functions
## that relate to a given GAM(M) to obtain a mode and a mean for the percentage
## of residuals found outside the simulate envelope based on the number (n) of 
## iterations specified by the user. The set.seed() function is solely used for
## results reproducibility. Here, n = 10 "hnp" iterations are performed.This code
## generates a list of results saved as an object referred to as "hnp_obj", and is
## applicable to any model for which the dfun(), sfun(), and ffun() helper functions
## have been previously defined.

set.seed(2025)

n <- 10

hnp_obj <- list()
for(i in 1:n) {
  hnp_obj[[i]]<-hnp(model, 
                    newclass = TRUE,
                    diagfun = dfun,
                    simfun = sfun,
                    fitfun = ffun,
                    how.many.out = TRUE,
                    plot.sim = FALSE)
}


## Calculate the percentage of residuals found outside the envelope of each 
## "hnp" iteration (listed in hnp_obj) and store it in the hnp_summary object

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 


## Obtain the mode (asymptotic value) from a density curve

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)


## Descriptive statistics, including the mean, with the Summarize() function of 
## the "FSA" package

Summarize(hnp_summary)


################### GAUSSIAN DISTRIBUTION AND IDENTITY LINK ####################

## Indicate the model being assessed, the family, and the data that were used to
## build the GAM(M) being investigated. Below, a fictive GAMM is provided as an 
## example where family = NO for the Gaussian case.

GAMM_GAUSSIAN <- gamlss(Y ~ pb(X1) + X2 + random(X3), 
                        family = NO, ## same as family = NO(mu.link = "identity")
                        data = CONTINUOUS) 

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


###################### GAMMA DISTRIBUTION AND LOG LINK #########################


GAMM_GAMMA <- gamlss(Y ~ pb(X1) + X2 + random(X3), 
                     family = GA, ## same as family = GA(mu.link = "log")
                     data = CONTINUOUS) 

model <- GAMM_GAMMA
family = GA
data = CONTINUOUS

dfun <- function(obj) resid(obj) 

sfun <- function(n, obj) { 
  mu <- obj$mu.fv 
  sig <- obj$sigma.fv 
  rGA(n, mu = mu, sigma = sig) 
}

ffun <- function(resp) gamlss(resp ~ pb(X1) + X2 + random(X3),
                              family = family,
                              data = data)


############ POISSON DISTRIBUTION AND ITS EXTENSIONS (LOG LINK) ################

## A specific "hnp" function was created for log-linear models analysing counts for
## GAM(M)s fitted with PO (Poisson), NBI (NB2), NBII (NB1), DPO (double Poisson),
## generalised Poisson (GPO), and Poisson inverse-Gaussian (PIG), which can all be
## assessed for adequacy with the hnp_gamlss_count() function described below, as
## well as any 3- and even 4-parameter distributions that are available in "gamlss".

hnp_gamlss_count <- function(model, d_fun, ...) {
  fam <- model$family[1]
  random_generator <- get(paste0("r", fam))
  parms <- model$parameters
  
  n <- length(model$y)
  
  mu_hat <- predict(model, what = "mu", type = "response")
  if("sigma" %in% parms) {
    sigma_hat <- predict(model, what = "sigma", type = "response")
  }
  if("nu" %in% parms) {
    nu_hat <- predict(model, what = "nu", type = "response")
  }
  if("tau" %in% parms) {
    tau_hat <- predict(model, what = "tau", type = "response")
  }
  
  if(length(parms) == 1) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat)
      return(y_new)
    }
  } else if(length(parms) == 2) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat, sigma = sigma_hat)
      return(y_new)
    }
  } else if(length(parms) == 3) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat, sigma = sigma_hat,
                                nu = nu_hat)
      return(y_new)
    }
  } else if(length(parms) == 4) {
    s_fun <- function(n, obj) {
      y_new <- random_generator(n, mu = mu_hat, sigma = sigma_hat,
                                nu = nu_hat, tau = tau_hat)
      return(y_new)
    }
  }
  
  if(missing(d_fun)) {
    d_fun <- function(obj) resid(obj)
  }
  
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


## The hnp_gamlss_count() function can then directly be applied, as in this example.

GAMM_POISSON <- gamlss(Y ~ pb(X1) + X2 + random(X3), 
                       family = PO, ## same as family = PO(mu.link = "log")
                       data = COUNT)

model <- GAMM_POISSON ## no need to define the family and data that were used

hnp_gamlss_count(model,
                 how.many.out = TRUE,
                 plot.sim = TRUE, 
                 paint = TRUE)


## To perform a number (n) of "hnp" iterations to get a mode and a mean, the
## code below need to be used.

set.seed(2025)

n <- 10

hnp_obj <- list()
for(i in 1:10) {
  hnp_obj[[i]] <- hnp_gamlss_count(model, 
                                how.many.out = TRUE,
                                plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

return_max(hnp_summary)

Summarize(hnp_summary)


############### BINOMIAL DISTRIBUTION AND ITS EXTENSIONS #######################

## Binomial distribution (BI) with logit link

GAMM_BI <- gamlss(Y ~ pb(X1) + random(X2), 
                  family = BI, ## same as family = BI(mu.link = "logit")
                  data = DISCRETE_PROP)

model <- GAMM_BI
family = BI
data = DISCRETE_PROP

dfun <- function(obj) resid(obj)

sfun <- function(n, obj) {
  mu <- obj$mu.fv
  bd <- obj$bd
  rBI(n, bd = bd, mu = mu)
}

ffun <- function(new_response) 
    gamlss(cbind(new_response, TOTAL - new_response) ~ pb(X1) + random(X2),
           family = family,
           data = CHARR)

hnp(model,
    newclass = TRUE,
    diagfun = dfun,
    simfun = sfun,
    fitfun = ffun,
    how.many.out = TRUE,
    paint = TRUE)


## Beta-binomial distribution (BB), zero-inflated binomial (ZIBI) and zero-inflated
## beta-binomial (ZIBB), here all fitted with the logit link (default)

GAMM_BB <- gamlss(Y ~ pb(X1) + random(X2), 
                  family = BB,
                  data = DISCRETE_PROP)

GAMM_ZIBI <- gamlss(Y ~ pb(X1) + random(X2), 
                    family = ZIBI,
                    data = DISCRETE_PROP)

GAMM_ZIBB <- gamlss(Y ~ pb(X1) + random(X2), 
                    family = ZIBB,
                    data = DISCRETE_PROP)


## The hnp() function can be applied directly on BB, ZIBI, and ZIBB, given that 
## these BI extensions can all be handled by the hnp package. However, the modelling
## discrete proportions with low to moderate sample sizes is more prone to analytical
## issues from resampling the data, such that the tryCatch() function and additional
## coding are also used to ensure that a result will be obtained at the end of this
## process. Here is an example for the adequacy assessment of GAMM_BB with "hnp".

model <- GAMM_BB

works <- FALSE
while (works == F) {
  tryCatch({
    hnp(model,
        how.many.out = TRUE,
        paint = TRUE)
    works <- TRUE  
  }, error = function(e) message("Trying again...")
  )
}


## For the ZIBI and ZIBB, modelling issues risks are increased due to the zero-
## inflated term. As such, even allowing to accommodate errors during the process
## may lead to an excessively long computational process. To allow to reach a result
## in a reasonable ammount of time, reducing the default number of simulations from
## the default 99 to 19 while also modifying the default level of confidence from
## 95% to 100% (see Moral et al. 2017) will allow in most instances to reach a result.
## Below, the argument sim = 19 and conf = 1 are specified for the GAMM_ZIBB.

model <- GAMM_ZIBB

works <- FALSE
while (works == F) {
  tryCatch({
    hnp(model,
        how.many.out = TRUE,
        paint = TRUE,
        sim = 19,
        conf = 1)
    works <- TRUE  
  }, error = function(e) message("Trying again...")
  )
}


## To perform two or more "hnp" iterations (here n = 10), with GAMM_BB as an
## example.

model <- GAMM_BB

set.seed(2025)

n <- 10

hnp_obj <- list()
works <- rep(FALSE, n)
for (i in 1:n) {
  while (works[i] == F) {
    tryCatch({
      hnp_obj[[i]] <- hnp(model,
                       how.many.out = TRUE,
                       plot.sim = FALSE)
      works[i] <- TRUE 
    }, error = function(e) message("Trying iteration ", i, " again...")
    )
  }
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100)

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary) 


## To perform two or more "hnp" iterations (here n = 10) while also specifying
## the argument sim = 19 and conf = 1 (Moral et al. 2017).

set.seed(2025)

n <- 10

hnp_obj <- list()
works <- rep(FALSE, n)
for (i in 1:n) {
  while (works[i] == F) {
    tryCatch({
      hnp_obj[[i]] <- hnp(model,
                       how.many.out = TRUE,
                       plot.sim = FALSE,
                       sim = 19,
                       conf = 1)
      works[i] <- TRUE 
    }, error = function(e) message("Trying iteration ", i, " again...")
    )
  }
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100)

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)

