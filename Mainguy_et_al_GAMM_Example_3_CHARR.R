## Required packages:

library(FSA)
library(gamlss)
library(gratia)
library(hnp)
library(itsadug)
library(MASS)
library(mgcv)
library(mgcViz)
library(MuMIn)
library(performance)
library(tidyverse)
library(visreg)

################################################################################
########## EXAMPLE 3: SPAWNING PROBABILITY OF ARCTIC CHARR FEMALES #############
################################################################################

## Load the CHARR dataset as an example with discrete proportions 

CHARR <- read_delim(
       "https://raw.githubusercontent.com/rachel-mc/GAMM_examples/XXXXXXXXXXXX")


## Binomial (BI) GAMM built in "mgcv" with a complementary log-log (cloglog) link.
## A smooth function is applied to BIN50 to check if additional non-linearities
## may be present in the binomial regression curve.

CHARR$RIVER <- as.factor(CHARR$RIVER)

m_CHARR_MGCV_BI <- gam(cbind(YES, TOTAL - YES) ~ s(BIN50) + s(RIVER, bs = "re"),
                       family = binomial(link = "cloglog"),
                       method = "ML",
                       data = CHARR)

summary(m_CHARR_MGCV_BI)


## As edf = 1 for s(BIN50), this fixed predictor is refitted as a parametric
## component instead under the same name.

m_CHARR_MGCV_BI <- gam(cbind(YES, TOTAL - YES) ~ BIN50 + s(RIVER, bs = "re"),
                       family = binomial(link = "cloglog"),
                       method = "ML",
                       data = CHARR)

summary(m_CHARR_MGCV_BI)


## REML estimation is now used given that this model will not be compared to others

m_CHARR_MGCV_BI_REML <- gam(cbind(YES, TOTAL - YES) ~ BIN50 + s(RIVER, bs = "re"),
                            family = binomial(link = "cloglog"),
                            method = "REML",
                            data = CHARR)

summary(m_CHARR_MGCV_BI_REML)


## Estimate the L50 from the fixed effects where b0 is the intercept coefficient
## and b1 the coefficient associated to BIN50 (i.e., slope). The equation for the
## estimation of the L50 with a cloglog link is detailed in the Supplementary
## Material of Mainguy et al. (2024, Fisheries Research)

b0 <- coef(m_CHARR_MGCV_BI_REML)[[1]]
b1 <- coef(m_CHARR_MGCV_BI_REML)[[2]]
((-b0 - 0.3665129) / b1)


## Estimate the L50 with the dose.p() function of "MASS" to also obtain the SE
## as estimated from the Delta method.

dose.p(m_CHARR_MGCV_BI_REML)


## Fitting the same BI GAMM with "gamlss" with BIN50 included as a parametric
## component. Note that only fitting family = BI corresponds to fitting the
## canonical logit link, which is equivalent to family = BI(mu.link = "logit").

m_CHARR_GAMLSS_BI <- gamlss(cbind(YES, TOTAL - YES) ~ BIN50 + random(RIVER),
                            family = BI(mu.link = "cloglog"), 
                            data = CHARR)

summary(m_CHARR_GAMLSS_BI)


## Estimate the L50 from the fixed effects. 
## The dose.p() function of "MASS" cannot be used with a "gamlss" object.

b0 <- coef(m_CHARR_GAMLSS_BI)[[1]]
b1 <- coef(m_CHARR_GAMLSS_BI)[[2]]
((-b0 - 0.3665129) / b1)


## Now 4 candidate GAMMs are created with "gamlss" which will use the same data
## but to create overdispersion, the random RIVER intercept will be voluntarily
## omitted, with the binomial (BI) distribution expected to no longer be adequate,
## and with the beta-binomial (BB), zero-inflated BI (ZIBI) and zero-inflated BB
## (ZIBB) possibly being better adapted under such condition with a percentage of
## zeros too. Models are identified with NRS (No Random Structure)

m_CHARR_GAMLSS_NRS_BI <- gamlss(cbind(YES, TOTAL - YES) ~ BIN50,
                                family = BI(mu.link = "cloglog"),
                                data = CHARR)

summary(m_CHARR_GAMLSS_NRS_BI)

b0 <- coef(m_CHARR_GAMLSS_NRS_BI)[[1]]
b1 <- coef(m_CHARR_GAMLSS_NRS_BI)[[2]]
((-b0 - 0.3665129) / b1)


m_CHARR_GAMLSS_NRS_BB <- gamlss(cbind(YES, TOTAL - YES) ~ BIN50,
                                family = BB(mu.link = "cloglog"),
                                data = CHARR)

summary(m_CHARR_GAMLSS_NRS_BB)

b0 <- coef(m_CHARR_GAMLSS_NRS_BB)[[1]]
b1 <- coef(m_CHARR_GAMLSS_NRS_BB)[[2]]
((-b0 - 0.3665129) / b1)


m_CHARR_GAMLSS_NRS_ZIBI <- gamlss(cbind(YES, TOTAL - YES) ~ BIN50,
                                  family = ZIBI(mu.link = "cloglog"),
                                  data = CHARR)

summary(m_CHARR_GAMLSS_NRS_ZIBI)

b0 <- coef(m_CHARR_GAMLSS_NRS_ZIBI)[[1]]
b1 <- coef(m_CHARR_GAMLSS_NRS_ZIBI)[[2]]
((-b0 - 0.3665129) / b1)


m_CHARR_GAMLSS_NRS_ZIBB <- gamlss(cbind(YES, TOTAL - YES) ~ BIN50,
                                  family = ZIBB(mu.link = "cloglog"),
                                  data = CHARR)

summary(m_CHARR_GAMLSS_NRS_ZIBB)

b0 <- coef(m_CHARR_GAMLSS_NRS_ZIBB)[[1]]
b1 <- coef(m_CHARR_GAMLSS_NRS_ZIBB)[[2]]
((-b0 - 0.3665129) / b1)


## Despite the fact that all 4 GAMMs converged, the one fitted with ZIBI yielded
## an estimate for the zero-inflated term that is associated with an infinite SE,
## clearly indicating modelling issues leading to an unreliable that was thus
## discarded. The lower L50 estimate for ZIBB is also suspicious.


############################# MODEL ADEQUACY ###################################

## Diagnostic plots with the gam.check() function of "mgcv"

gam.check(m_CHARR_MGCV_BI_REML)


## The same diagnostic plots with the appraise() function of "gratia" with a 
## simulated, theoretically-based envelope relying on deviance residuals (default)

appraise(m_CHARR_MGCV_BI_REML,
         method = "simulate")


## Worm plot with the worm_plot() function of "gratia" with a simulated envelope

worm_plot(m_CHARR_MGCV_BI_REML,
          method = "simulate",
          n_simulate = 100)


## Adequacy assessment with "hnp" for "mgcv": helper functions and 1 iteration.

model <- m_CHARR_MGCV_BI_REML
family <- binomial(link = "cloglog")
method <- "REML"
data <- CHARR
size <- CHARR$TOTAL

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

hnp(model,
    newclass = TRUE,
    simfun = sfun,
    fitfun = ffun,
    diagfun = dfun, 
    how.many.out = TRUE,
    paint = TRUE)


## 10 iterations

set.seed(2025)
hnp_obj <- list()
    for(i in 1:10) {
      hnp_obj[[i]] <- hnp(
                   model, 
                   newclass = TRUE, 
                   diagfun = dfun, 
                   simfun = sfun,
                   fitfun = ffun, 
                   how.many.out = TRUE, 
                   plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100)

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)


## Calculate the mgcViz score (mode) and mean for a given model based on 
## 100 iterations. A 95% "uncertainty interval" (Ui) is also calculated.
## This can take several seconds to run

model <- m_CHARR_MGCV_BI_REML
predictor <- "BIN50"

set.seed(2025)

viz_fun <- function(model, predictor) {
  viz <- getViz(model, nsim = 100)
  plot <- check1D(viz, predictor) + l_gridCheck1D(n = 100)
  diag_plot <- plot$ggObj
  est <- diag_plot$layers[[1]]$data
  ll <- diag_plot$layers[[3]]$data$ll
  ul <- diag_plot$layers[[3]]$data$ul
  diag <- cbind(est, ll, ul)
  diagnostic <- mutate(diag,
                       ll_diff = y - ll,
                       ul_diff = ul - y,
                       test = ll_diff * ul_diff)
  n <- length(diagnostic$test)
  i_n <- sum(diagnostic$test >= 0)
  mgcViz_perc <- i_n / n * 100
  mgcViz_perc
}

summary_mgcViz <- replicate(100, viz_fun(model, predictor))

mgcViz_score <- mean(summary_mgcViz)
lower_ui <- (quantile(summary_mgcViz, probs = 0.025))[[1]]
upper_ui <- (quantile(summary_mgcViz, probs = 0.975))[[1]]

cbind(mgcViz_score, lower_ui, upper_ui)


## Diagnostic plots with the plot() function of "gamlss"

plot(m_CHARR_GAMLSS_BI)

plot(m_CHARR_GAMLSS_NRS_BI)

plot(m_CHARR_GAMLSS_NRS_BB)

plot(m_CHARR_GAMLSS_NRS_ZIBB)


## Worm plot with the wp() function of "gamlss"

wp(m_CHARR_GAMLSS_BI)

wp(m_CHARR_GAMLSS_NRS_BI)

wp(m_CHARR_GAMLSS_NRS_BB)

wp(m_CHARR_GAMLSS_NRS_ZIBB)


## Bucket plot with the bp() function of "gamlss"

bp(m_CHARR_GAMLSS_BI)

bp(m_CHARR_GAMLSS_NRS_BI)

bp(m_CHARR_GAMLSS_NRS_BB)

bp(m_CHARR_GAMLSS_NRS_ZIBB)


## Detrended-transformed Owen's Plot (DTOP) with the dtop() function of "gamlss"

dtop(m_CHARR_GAMLSS_BI)

dtop(m_CHARR_GAMLSS_NRS_BI)

dtop(m_CHARR_GAMLSS_NRS_BB)

dtop(m_CHARR_GAMLSS_NRS_ZIBB)


## Adequacy assessment with "hnp" for BI GAMM fitted with "gamlss" (1 iteration)
## The BB, ZIBI, and ZIBB can all be already directly handled by "hnp". The code
## to be used is shown further below.

## GAMM with binomial (BI) distribution and including a random RIVER intercept

model <- m_CHARR_GAMLSS_BI
family = BI(mu.link = "cloglog")
data = CHARR

dfun <- function(obj) resid(obj)

sfun <- function(n, obj) {
  mu <- obj$mu.fv
  bd <- obj$bd
  rBI(n, bd = bd, mu = mu)
}

ffun <- function(new_response) 
        gamlss(cbind(new_response, TOTAL - new_response) ~ BIN50 + random(RIVER),
               family = family,
               data = data)

hnp(model,
    newclass = TRUE,
    diagfun = dfun,
    simfun = sfun,
    fitfun = ffun,
    how.many.out = TRUE,
    paint = TRUE)

set.seed(2025)
hnp_obj <- list()
  for(i in 1:10) {
      hnp_obj[[i]] <- hnp(model,
                      newclass = TRUE,
                      diagfun = dfun,
                      simfun = sfun,
                      fitfun = ffun,
                      how.many.out = TRUE,
                      plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100) 

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)


## BI GAM but the random RIVER intercept is omitted

model <- m_CHARR_GAMLSS_NRS_BI
family = BI(mu.link = "cloglog")
data = CHARR

dfun <- function(obj) resid(obj)

sfun <- function(n, obj) {
  mu <- obj$mu.fv
  bd <- obj$bd
  rBI(n, bd = bd, mu = mu)
}

ffun <- function(new_response) 
        gamlss(cbind(new_response, TOTAL - new_response) ~ BIN50,
               family = family,
               data = data)

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
for(i in 1:n) {
  hnp_obj[[i]] <- hnp(model,
                   newclass = TRUE,
                   diagfun = dfun,
                   simfun = sfun,
                   fitfun = ffun,
                   how.many.out = TRUE,
                   plot.sim = FALSE)
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100) 

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)


## 1 "hnp" iteration for BB with the tryCatch() function such that "hnp" is able
## to continue to perform the 99 simulations despite modelling errors that may 
## occur in "gamlss" during this process

model <- m_CHARR_GAMLSS_NRS_BB

works <- FALSE
while (works == F) {
  tryCatch({
    hnp(model, how.many.out = TRUE, paint = TRUE)
    works <- TRUE  
  }, error = function(e) message("Trying again...")
  )
}


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

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100)

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
  }

round(return_max(hnp_summary), 2)

Summarize(hnp_summary) 


## Zero-inflated models are more prone to modelling errors, such that the default
## number of simulations per iteration is reduced from 99 to 19 and the default 
## level of confidence (conf) is increased from 95% to 100% (conf = 100). See
## Moral et al. (2017) for details. Here for ZIBB only as ZIBI was problematic.

model <- m_CHARR_GAMLSS_NRS_ZIBB 

works <- FALSE
while (works == F) {
  tryCatch({
    hnp(model, how.many.out = TRUE, paint = TRUE, sim = 19, conf = 1)
    works <- TRUE  
  }, error = function(e) message("Trying again...")
  )
}


set.seed(2025)

n <- 10

hnp_obj <- list()
works <- rep(FALSE, n)
  for (i in 1:n) {
    while (works[i] == F) {
      tryCatch({
        hnp_obj[[i]] <- hnp(m_CHARR_GAMLSS_NRS_ZIBB, 
                        how.many.out = TRUE, 
                        plot.sim = FALSE, 
                        sim = 19,
                        conf = 1)
        works[i] <- TRUE 
    }, error = function(e) message("Trying iteration ", i, " again...")
   )
  }
}

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100)

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary) 


## mgcViz score

model <- m_CHARR_MGCV_BI_REML
predictor <- "BIN50"

set.seed(2025)
viz_fun <- function(model, predictor) {
  viz <- getViz(model, nsim = 100)
  plot <- check1D(viz, predictor) + l_gridCheck1D(n = 100)
  diag_plot <- plot$ggObj
  est <- diag_plot$layers[[1]]$data
  ll <- diag_plot$layers[[3]]$data$ll
  ul <- diag_plot$layers[[3]]$data$ul
  diag <- cbind(est, ll, ul)
  diagnostic <- mutate(diag,
                       ll_diff = y - ll,
                       ul_diff = ul - y,
                       test = ll_diff * ul_diff)
  n <- length(diagnostic$test)
  i_n <- sum(diagnostic$test >= 0)
  mgcViz_perc <- i_n / n * 100
  mgcViz_perc
}
summary_mgcViz <- replicate(100, viz_fun(model, predictor))
return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}
mgcViz_score <- return_max(summary_mgcViz)
mgcViz_mean <- mean(summary_mgcViz)
lower_ui <- (quantile(summary_mgcViz, probs = 0.025))[[1]]
upper_ui <- (quantile(summary_mgcViz, probs = 0.975))[[1]]

cbind(mgcViz_score, mgcViz_mean, lower_ui, upper_ui)

## check concurvity

concurvity(m_CHARR_MGCV_BI_REML)


############################# MODEL SELECTION ##################################

## The implemented LR.test() function is used to compare the two GAMs.

LR.test(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_NRS_ZIBB)


## An information-theoretic approach based on AICc is used to compared the GAMs

model.sel(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_NRS_ZIBB)


## The test based on the Bayes Factor (BF) in "performance" 

test_bf(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_NRS_ZIBB)


##################### MODEL PREDICTIVE PERFORMANCE #############################

## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the final retained model relying on REML estimation

model <- m_CHARR_MGCV_BI_REML 

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
Total_df <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - Total_df) * (100 - D2))
cbind(D2, D2_adj)


## Estimate the pseudo-R2 of Nagelkerke (1991) for the BI, BB and ZIBB GAMM built
## with "gamlss"

Rsq(m_CHARR_GAMLSS_BI)

Rsq(m_CHARR_GAMLSS_NRS_BI)

Rsq(m_CHARR_GAMLSS_NRS_BB)

Rsq(m_CHARR_GAMLSS_NRS_ZIBB)


############################# MODEL PREDICTIONS ################################

## predictions with "mgcv" (Fig. 3a)

nd_CHARR <- data.frame(
  BIN50 = seq(175, 675, 1),
  RIVER = "VOLTZ")

model <- m_CHARR_MGCV_BI_REML

BIN50 <- nd_CHARR$BIN50

fitted <- predict(
  model, nd_CHARR, type = "link", exclude = "s(RIVER)", se.fit = TRUE)
pred <- 1-exp(-exp(fitted$fit))
lower <- 1-exp(-exp(fitted$fit - 1.96 * fitted$se.fit))
upper <- 1-exp(-exp(fitted$fit + 1.96 * fitted$se.fit))
cbind(BIN50, pred, lower, upper)


## predictions with "gamlss" (Fig. 3a)

nd_CHARR_ALL <- data.frame(
  BIN50 = seq(175, 675, 1))

nd_CHARR_AIPPARUSIK <- data.frame(
  BIN50 = seq(175, 675, 1),
  RIVER = "AIPPARUSIK")

nd_CHARR_TASIALLUJUAK <- data.frame(
  BIN50 = seq(175, 675, 1),
  RIVER = "TASIALLUJUAK")

nd_CHARR_VOLTZ <- data.frame(
  BIN50 = seq(175, 675, 1),
  RIVER = "VOLTZ")

model <- m_CHARR_GAMLSS_BI

fitted_AIPPARUSIK <- predict(model, 
                             what = "mu", 
                             newdata = nd_CHARR_AIPPARUSIK,
                             type = "response")

fitted_TASIALLUJUAK <- predict(model, 
                               what = "mu", 
                               newdata = nd_CHARR_TASIALLUJUAK,
                               type = "response")

fitted_VOLTZ <- predict(model, 
                        what = "mu",
                        newdata = nd_CHARR_VOLTZ,
                        type = "response")

fitted_ALL <- cbind(fitted_AIPPARUSIK, fitted_TASIALLUJUAK, fitted_VOLTZ)

fitted_ALL <- transform(fitted_ALL,
                        average = ((fitted_AIPPARUSIK + fitted_TASIALLUJUAK +
                                    fitted_VOLTZ) / 3))

pred_average <- fitted_ALL$average 

cbind(nd_CHARR_ALL, pred_average)


## Predictions with "gamlss" without the random RIVER intercept (Fig. 3b)
## Predicted values are obtained from the observed data, such that se.fit = TRUE
## can be used. The BB and ZIBB GAMs can have their predictions generated the
## same way.

## BI

model <- m_CHARR_GAMLSS_NRS_BI

fitted <- predict(model,
                  what = "mu",
                  type = "link",
                  se.fit = TRUE)

estimate <- 1 - exp(-exp(fitted$fit))
lower_ci <- 1 - exp(-exp(fitted$fit - 1.96 * fitted$se.fit))
upper_ci <- 1 - exp(-exp(fitted$fit + 1.96 * fitted$se.fit))
cbind(CHARR, estimate, lower_ci, upper_ci)

