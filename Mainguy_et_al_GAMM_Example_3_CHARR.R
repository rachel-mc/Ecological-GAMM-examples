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

file_name <- "https://raw.githubusercontent.com/rachel-mc/Ecological-GAMM-examples/main/CHARR.txt"
CHARR <- read.delim(file_name)


# Model Fitting -----------------------------------------------------------


## Binomial (BI) GAMM built in "mgcv" with a complementary log-log (cloglog) link.
## A smooth function is applied to BIN50 to check if additional non-linearities
## are present in the binomial regression curve.

CHARR$RIVER <- as.factor(CHARR$RIVER)

m_CHARR_MGCV_BI <- gam(cbind(YES, TOTAL - YES) ~ s(BIN50) + s(RIVER, bs = "re"),
                       family = binomial(link = "cloglog"),
                       method = "ML",
                       data = CHARR)

summary(m_CHARR_MGCV_BI)

## As edf = 1 for s(BIN50), this fixed predictor is refitted as a parametric
## component instead with the same model name.

m_CHARR_MGCV_BI <- gam(cbind(YES, TOTAL - YES) ~ BIN50 + s(RIVER, bs = "re"),
                       family = binomial(link = "cloglog"),
                       method = "ML",
                       data = CHARR)

summary(m_CHARR_MGCV_BI)

## REML estimation is now used given that this model will not be compared to others.

m_CHARR_MGCV_BI_REML <- gam(cbind(YES, TOTAL - YES) ~ BIN50 + s(RIVER, bs = "re"),
                            family = binomial(link = "cloglog"),
                            method = "REML",
                            data = CHARR)

summary(m_CHARR_MGCV_BI_REML)

## Estimate the L50 using the fixed effects where b0 is the intercept coefficient
## and b1 is the coefficient associated with BIN50 (i.e., slope). The equation to
## estimate the L50 for a cloglog link is detailed in the Supplementary
## Material of Mainguy et al. (2024, Fisheries Research).

b0 <- coef(m_CHARR_MGCV_BI_REML)[[1]]
b1 <- coef(m_CHARR_MGCV_BI_REML)[[2]]
((-b0 - 0.3665129) / b1)

## Estimate the L50 using the dose.p() function of "MASS" to also obtain the 
## standard error (SE), which is calculated using the Delta method.

dose.p(m_CHARR_MGCV_BI_REML)

## Fitting the same BI GAMM in "gamlss" with BIN50 included as a parametric
## component. Note that the default canonical link is logit, so family = BI 
## actually corresponds to family = BI(mu.link = "logit").

m_CHARR_GAMLSS_BI <- gamlss(cbind(YES, TOTAL - YES) ~ BIN50 + random(RIVER),
                            family = BI(mu.link = "cloglog"), 
                            data = CHARR)

summary(m_CHARR_GAMLSS_BI)

## Estimate the L50 using the fixed effects. 
## The dose.p() function of "MASS" cannot be used with a "gamlss" object.

b0 <- coef(m_CHARR_GAMLSS_BI)[[1]]
b1 <- coef(m_CHARR_GAMLSS_BI)[[2]]
((-b0 - 0.3665129) / b1)

## Now, 4 candidate GAMMs are created in "gamlss" which will use the same data as before
## but the random intercept RIVER will be voluntarily omitted to create
## overdispersion. The binomial (BI) distribution is expected to no longer be adequate.
## The beta-binomial (BB), zero-inflated BI (ZIBI) and zero-inflated BB
## (ZIBB) will possibly adapt better to such conditions that also have a percentage
## of zeros. Models are named with NRS (No Random Structure).

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
## an 'infinite' SE for the estimate of the zero-inflated term.
## This indicates modelling issues that would lead to unreliable inferences, 
## thus this model was discarded. The lower L50 estimate for ZIBB is also suspicious.


# Model Adequacy ----------------------------------------------------------


## Diagnostic plots using the gam.check() function of "mgcv".

gam.check(m_CHARR_MGCV_BI_REML)

## The same diagnostic plots using the appraise() function of "gratia" with a 
## simulated, theoretically-based envelope relying on deviance residuals (default).

appraise(m_CHARR_MGCV_BI_REML,
         method = "simulate")

## Worm plot using the worm_plot() function of "gratia" with a simulated envelope.

worm_plot(m_CHARR_MGCV_BI_REML,
          method = "simulate",
          n_simulate = 100)

## Adequacy assessment using "hnp" for "mgcv": helper functions and 1 iteration.

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

dfun <- function(obj) resid(obj, type = "deviance")

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
  hnp_obj[[i]] <- hnp(model, 
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

## Calculate the mgcViz score (mode and mean) for a given model based on 
## 100 iterations. A 95% "uncertainty interval" (ui) is also calculated.
## This can take several seconds to run.

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

summary_mgcViz <- replicate(100, viz_fun(model = m_CHARR_MGCV_BI_REML,
                                         predictor = "BIN50"))

mgcViz_score <- mean(summary_mgcViz)
lower_ui <- (quantile(summary_mgcViz, probs = 0.025))[[1]]
upper_ui <- (quantile(summary_mgcViz, probs = 0.975))[[1]]

cbind(mgcViz_score, lower_ui, upper_ui)

## Diagnostic plots using the plot() method of "gamlss"

plot(m_CHARR_GAMLSS_BI)

plot(m_CHARR_GAMLSS_NRS_BI)

plot(m_CHARR_GAMLSS_NRS_BB)

plot(m_CHARR_GAMLSS_NRS_ZIBB)

## Worm plot using the wp() function of "gamlss"

wp(m_CHARR_GAMLSS_BI)

wp(m_CHARR_GAMLSS_NRS_BI)

wp(m_CHARR_GAMLSS_NRS_BB)

wp(m_CHARR_GAMLSS_NRS_ZIBB)

## Bucket plot using the bp() function of "gamlss"

bp(m_CHARR_GAMLSS_BI)

bp(m_CHARR_GAMLSS_NRS_BI)

bp(m_CHARR_GAMLSS_NRS_BB)

bp(m_CHARR_GAMLSS_NRS_ZIBB)

## Detrended-transformed Owen's Plot (DTOP) using the dtop() function of "gamlss"

dtop(m_CHARR_GAMLSS_BI)

dtop(m_CHARR_GAMLSS_NRS_BI)

dtop(m_CHARR_GAMLSS_NRS_BB)

dtop(m_CHARR_GAMLSS_NRS_ZIBB)

## Adequacy assessment using "hnp" for BI GAMM fitted in "gamlss" (1 iteration).
## The BB, ZIBI, and ZIBB can all be directly handled by "hnp". 
## The code to be used is shown below.

## GAMM including the random intercept RIVER and fitted with the binomial (BI) distribution

model <- m_CHARR_GAMLSS_BI
family = BI(mu.link = "cloglog")
data = CHARR

dfun <- function(obj) resid(obj)

sfun <- function(n, obj) {
  mu <- obj$mu.fv
  bd <- obj$bd
  rBI(n, bd = bd, mu = mu)
  }

ffun <- function(new_response) {
  gamlss(cbind(new_response, TOTAL - new_response) ~ BIN50 + random(RIVER),
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

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)

## BI GAM (the random intercept RIVER is omitted)

model <- m_CHARR_GAMLSS_NRS_BI
family = BI(mu.link = "cloglog")
data = CHARR

# Same diagnostic and simulation function as above
# New fitting function
ffun <- function(new_response) {
  gamlss(cbind(new_response, TOTAL - new_response) ~ BIN50,
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
n <- 10 # can be altered depending on time constraints

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

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)

## 1 "hnp" iteration for BB using the tryCatch() function such that "hnp" is able
## to perform 99 simulations despite modelling errors that may occur in "gamlss"
## during this process.

model <- m_CHARR_GAMLSS_NRS_BB

run_complete <- FALSE
while (run_complete == F) {
  tryCatch( {
    hnp(model, how.many.out = TRUE, paint = TRUE)
    run_complete <- TRUE  
  },
  error = function(e) message("Trying the \"hnp\" iteration again...")
  )
}

set.seed(2025)

n <- 10

hnp_obj <- list()
run_complete <- rep(FALSE, n)
for (i in 1:n) {
  while (run_complete[i] == F) {
    tryCatch( {
      hnp_obj[[i]] <- hnp(model, how.many.out = TRUE, plot.sim = FALSE) 
      run_complete[i] <- TRUE 
      }, 
      error = function(e) message("Trying \"hnp\" iteration ", i, " again...")
    )
  }
  }

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100)

round(return_max(hnp_summary), 2)

Summarize(hnp_summary) 

## Zero-inflated models are more prone to modelling errors, so the default
## number of simulations per iteration is reduced from 99 to 19 and the default 
## level of confidence is increased from 95% to 100% (conf = 1) here. See
## Moral et al. (2017) for details. We use ZIBB only below as ZIBI was problematic.

model <- m_CHARR_GAMLSS_NRS_ZIBB 

run_complete <- FALSE
while (run_complete == F) {
  tryCatch( {
    hnp(model, how.many.out = TRUE, paint = TRUE, sim = 19, conf = 1)
    run_complete <- TRUE  
  },
  error = function(e) message("Trying the \"hnp\" iteration again...")
  )
}

set.seed(2025)

n <- 10

hnp_obj <- list()
run_complete <- rep(FALSE, n)
for (i in 1:n) {
  while (run_complete[i] == F) {
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

## mgcViz score

set.seed(2025)

summary_mgcViz <- replicate(100, viz_fun(model = m_CHARR_MGCV_BI_REML, predictor = "BIN50"))
mgcViz_score <- return_max(summary_mgcViz)
mgcViz_mean <- mean(summary_mgcViz)
lower_ui <- (quantile(summary_mgcViz, probs = 0.025))[[1]]
upper_ui <- (quantile(summary_mgcViz, probs = 0.975))[[1]]

cbind(mgcViz_score, mgcViz_mean, lower_ui, upper_ui)

## Check concurvity

concurvity(m_CHARR_MGCV_BI_REML)


# Model Selection ---------------------------------------------------------


## The implemented LR.test() function in "gamlss" is used to compare the two GAMs.

LR.test(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_NRS_ZIBB)

## An information-theoretic approach based on AICc to compare the GAMs.

model.sel(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_NRS_ZIBB)

## The test based on the Bayes Factor (BF) in "performance" to compare the GAMs. 

test_bf(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_NRS_ZIBB)


# Predictive Performance --------------------------------------------------


## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the BI GAMM in "mgcv" that relies on REML estimation.

model <- m_CHARR_MGCV_BI_REML 

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
Total_df <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - Total_df) * (100 - D2))
cbind(D2, D2_adj)

## Estimate the pseudo-R2 of Nagelkerke (1991) for the BI, BB, and ZIBB GAMMs built
## in "gamlss".

Rsq(m_CHARR_GAMLSS_BI)

Rsq(m_CHARR_GAMLSS_NRS_BI)

Rsq(m_CHARR_GAMLSS_NRS_BB)

Rsq(m_CHARR_GAMLSS_NRS_ZIBB)


# Model Predictions -------------------------------------------------------


## Predictions with "mgcv" (Fig. 3a).

# New data
nd_CHARR <- data.frame(BIN50 = seq(175, 675, 1),
                       RIVER = "VOLTZ")

model <- m_CHARR_MGCV_BI_REML

BIN50 <- nd_CHARR$BIN50

fitted <- predict(model, nd_CHARR, type = "link", exclude = "s(RIVER)", se.fit = TRUE)
pred <- 1 - exp(-exp(fitted$fit))
lower <- 1 - exp(-exp(fitted$fit - 1.96 * fitted$se.fit))
upper <- 1 - exp(-exp(fitted$fit + 1.96 * fitted$se.fit))
cbind(BIN50, pred, lower, upper)

## Predictions with "gamlss" (Fig. 3a).

# New data
nd_CHARR_ALL <- data.frame(BIN50 = seq(175, 675, 1))

nd_CHARR_AIPPARUSIK <- data.frame(BIN50 = seq(175, 675, 1),
                                  RIVER = "AIPPARUSIK")

nd_CHARR_TASIALLUJUAK <- data.frame(BIN50 = seq(175, 675, 1),
                                    RIVER = "TASIALLUJUAK")

nd_CHARR_VOLTZ <- data.frame(BIN50 = seq(175, 675, 1),
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

pred_average <- rowMeans(fitted_ALL)

cbind(nd_CHARR_ALL, pred_average)

## Predictions for "gamlss" models without the random intercept RIVER (Fig. 3b).
## Predicted values are obtained from the observed data so that se.fit = TRUE
## can be used. The BB and ZIBB GAMs can have their predictions generated in the
## same way, as follows.

## BI example
model <- m_CHARR_GAMLSS_NRS_BI

fitted <- predict(model,
                  what = "mu",
                  type = "link",
                  se.fit = TRUE)

estimate <- 1 - exp(-exp(fitted$fit))
lower_ci <- 1 - exp(-exp(fitted$fit - 1.96 * fitted$se.fit))
upper_ci <- 1 - exp(-exp(fitted$fit + 1.96 * fitted$se.fit))
cbind(CHARR, estimate, lower_ci, upper_ci)


# References --------------------------------------------------------------


## Mainguy J, Bélanger M, Ouellet-Cauchon G et al. (2024) Monitoring reproduction
##    in fish: assessing the adequacy of ogives and the predicted uncertainty of
##    their L50 estimates for more reliable biological inferences. 
##    Fish Res 269:106863. https://doi.org/10.1016/j.fishres.2023.106863

## Moral RA, Hinde J, Demétrio CGB (2017) Half-normal plots and overdispersed
##    models in R: the hnp package. J Stat Softw 81:1–23. https://doi.org/10.18637/jss.v081.i10

## Nagelkerke NJ (1991) A note on a general definition of the coefficient of determination.
##    Biometrika 78:691-692. https://doi.org/10.1093/biomet/78.3.691
