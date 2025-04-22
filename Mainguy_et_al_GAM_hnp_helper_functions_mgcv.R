## Required packages:

library(FSA)
library(hnp)
library(mgcv)


################################################################################
##################### "HNP" HELPER FUNCTIONS FOR "MGCV" ########################
################################################################################

## Up to now, "hnp" was unable to assess the adequacy of GAM(M)s fitted with "mgcv".
## Three helper functions are now made available for this. These are the dfun(),
## sfun(), and ffun() functions, which are described below. Note that the dfun()
## function remains the same regardless of the GAMM being assessed.

## dfun(): a helper function to extract the deviance residuals

dfun <- function(obj) resid(obj, type = "deviance")


## sfun(): A helper function to perform empirically-derived simulations to determine
## how the residuals should behave under the distributional assumptions upon 
## which the GAM(M) was fitted. This helper function is specific to the family
## that was used.


sfun <- function(n,obj){
  y <- ## SPECIFIC TO THE FAMILY CONSIDERED. SEE EXAMPLES BELOW.
  return(y)
  }


## ffun(): A helper function to refit the model 99 times (default value)
## The ffun() is also specific to the model being assessed. A fictive GAM is 
## presented below to illustrate the parametrisation of the ffun() helper
## function with Y = response variable and X1 = predictor, to which a smooth
## function, s(), is applied.

GAM_GAUSSIAN <- gam(Y ~ s(X1), 
                    family = gaussian ## same as family = gaussian(link = "identity"),
                    method = "GCV.Cp",
                    data = CONTINUOUS) 


## To define the ffun() helper function that is specific to GAM_GAUSSIAN, the model,
## family, method, and data used, are first defined.

model <- GAM_GAUSSIAN
family <- gaussian(link = "identity")
method <- "GCV.Cp"
data <- CONTINUOUS


## Then, the ffun() helper function is parametrised by indicating the right-hand
## part of the model next to "resp ~" which here is "s(X1)".

ffun <- function(resp) {
  gam(resp ~ s(X1),
      family = family,
      method = method,
      data = data)
      }


############ THE HNP() FUNCTION BELOW CAN BE USED WITH ALL GAM(M)S #############

## The hnp() function can now be used to assess the adequacy of a GAM by
## relying on the three previous helper functions that were just defined. Here, 
## a single hnp iteration will be performed based on 99 simulations (default).
## The argument how.many.out = TRUE will indicate the percentage of residuals found
## outside the simulated envelope of the half-normal plot. The corresponding half-
## normal plot will also be shown with the argument plot = TRUE. Using the argument
## paint = TRUE will show in red the residuals found outside the envelope.

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE)


###### THE ITERATIVE PROCESS DESCRIBED BELOW CAN BE USED WITH ALL GAM(M)S ######

## The specific dfun(), sfun(), ffun() helper functions can be used to obtain a
## mode and a mean for the percentage of residuals found outside the simulate
## envelope given the number (n) of iterations that are specified by the user. 
## The set.seed() function is solely used for results reproducibility, with 10 
## "hnp" iterations being performed below. The code generates a list of results
## saved under the object "hnp_obj" and is applicable to any model for which the
## dfun(), sfun(), and ffun() helper functions have been previously defined.

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
## "hnp" iteration and store it in the object hnp_summary

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

## Below, a fictive GAMM for which the deviance residuals are expected to approximate
## a normal distribution is provided as an example.

GAMM_GAUSSIAN <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                    family = gaussian(link = "identity"),
                    method = "REML",
                    data = CONTINUOUS) 

model <- GAMM_GAUSSIAN 
family <- gaussian(link = "identity")
method <- "REML"  
data <- CONTINUOUS


## helper function to extract the deviance residuals

dfun <- function(obj) resid(obj, type = "deviance")


## helper function to simulate from the data how the residuals should behave 
## given the distributional assumptions of the model being investigated. 
## The sfun() function is specific to the family used, with here family = gaussian

sfun <- function(n,obj){
   y <- rnorm(nrow(data),
       sd = sqrt(model$sig2),
       mean = predict(model, type = "response"))
   return(y)
}


## helper function to refit the model 99 times at each iteration 

ffun <- function(resp) {
  gam(resp ~ s(X1) + X2 + s(X3, bs = "re"),
      family = family,
      method = method,
      data = data)
}


###################### GAMMA DISTRIBUTION AND LOG LINK #########################

GAMM_GAMMA <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                     family = Gamma, ## same as Gamma(link = "log")
                     method = "REML",
                     data = CONTINUOUS) 

model <- GAMM_GAMMA
family <- Gamma(link = "log")
method <- "REML"  
data <- CONTINUOUS

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n, obj) simulate(obj)[,1]

ffun <- function(resp) gam(resp ~ s(X1) + X2 + s(X3, bs = "re"),
                           family = family,
                           method = method,
                           data = data)


##################### POISSON DISTRIBUTION AND LOG LINK ######################## 

GAMM_POISSON <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                     family = poisson, ## same as family = poisson(link = "log")
                     method = "REML",
                     data = COUNT) 

model <- GAMM_POISSON 
family <- poisson
method <- "REML"  
data <- COUNT

dfun<-function(obj) resid(obj, type = "deviance")

sfun<-function(n,obj){
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


############### TYPE-II NEGATIVE BINOMIAL (NB2) AND LOG LINK ###################

GAMM_NB2 <- gam(Y ~ s(X1) + X2 + s(X3, bs = "re"), 
                    family = nb, ## same as family = nb(theta = NULL, link = "log")
                    method = "REML",
                    data = COUNT) 

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


########## BINOMIAL DISTRIBUTION AND LOGIT, PROBIT, OR CLOGLOG LINK ############
########### ONLY FOR BINOMIAL GAM(M) ANALYSING DISCRETE PROPORTIONS ############

GAMM_BINOMIAL_LOGIT <- gam(cbind(YES, TOTAL - YES) ~ s(X1) + X2 + s(X3, bs = "re"), 
                           family = binomial, ## same as binomial(link = "logit")
                           method = "REML",
                           data = DISCRETE_PROP) 


## For discrete proportions, size must also be provided (i.e., denominator). In
## this fictive dataset (i.e., DISCRETE_PROP), it is referred to as the TOTAL variable.

model <- GAMM_BINOMIAL_LOGIT
family <- binomial
method <- "REML"  
data <- DISCRETE_PROP
size = DISCRETE_PROP$TOTAL

sfun <- function(n, obj) {
  p_hat <- predict(obj, type = "response")
  y <- rbinom(n = n, size = size, prob = p_hat)
  return(y)
}

dfun <- function(obj) {
  resid(obj, type = "deviance")
}

ffun <- function(new_response) {
  gam(cbind(new_response, TOTAL - new_response) ~ BIN50 + s(RIVER, bs = "re"),
      family = family,
      method = method,
      data = data)
}

