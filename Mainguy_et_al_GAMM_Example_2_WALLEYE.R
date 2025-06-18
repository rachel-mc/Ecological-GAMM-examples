## Required packages:

library(FSA)
library(gam.hp)
library(gratia)
library(hnp)
library(itsadug)
library(mgcv)
library(mgcViz)
library(MuMIn)
library(performance)
library(tidyverse)
library(visreg)

################################################################################
#### EXAMPLE 2: CPUE TEMPORAL TREND OF SMALL ST. LAWRENCE RIVER WALLEYES #######
################################################################################

## Load the WALLEYE dataset as an example with counts 

file_name <- 
"https://raw.githubusercontent.com/rachel-mc/Ecological-GAMM-examples/main/WALLEYE.txt"

WALLEYE <- read.delim(file_name)


## GENERALIZED ADDITIVE MIXED-EFFECTS MODELING (GAMM) --------------------------

## Poisson GAMM with the default settings of the smooth function, s(), for all
## predictors, i.e., thin plate regression splines (bs = "tp") and k = 10, with
## Maximum Likelihood (ML) estimation and AREA fitted as a random effect in "mgcv".

WALLEYE$AREA <- as.factor(WALLEYE$AREA)

m_WALLEYE_POISSON_FULL <- gam(N ~ s(YEAR) + s(CONDUCT) + s(D1AUG) + s(DEPTH)
                              + s(EFFORT) + s(TEMP) + s(TURBID) + s(AREA, bs = "re"),
                              family = poisson,
                              method = "ML",
                              data = WALLEYE)

summary(m_WALLEYE_POISSON_FULL)

## Check for likely overdispersion issues (see text) in the Poisson GAMM.

check_overdispersion(m_WALLEYE_POISSON_FULL)

## Check for zero-inflation issues (see text) in the Poisson GAMM.

check_zeroinflation(m_WALLEYE_POISSON_FULL)

## Attempt to handle both overdispersion and zero-inflation with the type-II 
## negative binomial (NB2) extension instead in "mgcv".

m_WALLEYE_NB2_FULL <- gam(N ~ s(YEAR) + s(CONDUCT) + s(D1AUG) + s(DEPTH) +
                          s(EFFORT) + s(TEMP) + s(TURBID) + s(AREA, bs = "re"),
                          family = nb,
                          method = "ML",
                          data = WALLEYE)

summary(m_WALLEYE_NB2_FULL)

## To first determine whether this NB2 GAMM fitted in "mgcv" could be reduced by
## removing some predictors that may explain little to no variation in the response
## variable, a variable-selection approach which more heavily penalizes each smooth
## term with the argument select = TRUE will be used.

m_WALLEYE_NB2_SELECT <- gam(N ~ s(YEAR) + s(CONDUCT) + s(D1AUG) + s(DEPTH) +
                            s(EFFORT) + s(TEMP) + s(TURBID) + s(AREA, bs = "re"),
                            family = nb,
                            method = "ML",
                            select = TRUE,
                            data = WALLEYE)

summary(m_WALLEYE_NB2_SELECT)

## The predictors s(CONDUCT), s(EFFORT), and s(TEMP) have their respective
## effective degrees of freedom (edf) shrunk toward zero (0), indicating that
## these three covariates fitted as smooth terms can be discarded. The NB2 GAMM
## is refitted without these predictors and referred to as m_WALLEYE_REDUCED_NB2.

m_WALLEYE_NB2_REDUCED <- gam(N ~ s(YEAR) + s(D1AUG) + s(DEPTH) + s(TURBID) +
                                  s(AREA, bs = "re"),
                                  family = nb,
                                  method = "ML",
                                  data = WALLEYE)

summary(m_WALLEYE_NB2_REDUCED)

## The smooth terms s(D1AUG) and s(DEPTH) now both exhibit edf = 1, indicating 
## that both predictors can be included as parametric components instead. A final
## NB2 GAMM is thus fitted with fewer smooth terms.

m_WALLEYE_NB2_FINAL <- gam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID) +
                           s(AREA, bs = "re"),
                           family = nb,
                           method = "ML",
                           data = WALLEYE)

summary(m_WALLEYE_NB2_FINAL)

## Check if overdispersion and zero-inflation were sufficiently accounted for in
## the retained final NB2 GAMM, which can only be done for the one built in 
## "mgcv" when the NB2 extension is used. Note that for the zero-inflation test,
## a p-value is provided when the NB2 extension is used.

check_overdispersion(m_WALLEYE_NB2_FINAL)

check_zeroinflation(m_WALLEYE_NB2_FINAL)


## MODEL ADEQUACY --------------------------------------------------------------

## Diagnostic plots with the gam.check() function of "mgcv" (final GAMM NB2)

gam.check(m_WALLEYE_NB2_FINAL)

## The same diagnostic plots using the appraise() function of "gratia" instead 
## with a simulated envelope, relying on deviance residuals.

appraise(m_WALLEYE_NB2_FINAL, method = "simulate")

## Worm plot using the worm_plot() function of "gratia" with a simulated envelope
## based on 1000 simulations instead of the 50 by default.

worm_plot(m_WALLEYE_NB2_FINAL, method = "simulate", n_simulate = 1000)

## Rootogram with the rg() function of "gratia"

rg <- rootogram(m_WALLEYE_NB2_FINAL)
draw(rg)

## Adequacy assessment based on half-normal plots with "hnp". This can take many
## seconds to a few minutes to run.

model <- m_WALLEYE_NB2_FINAL
family <- nb
method <- "ML"
data <- WALLEYE

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n, obj) {
  y <- rnbinom(nrow(data),
               size = model$family$getTheta(TRUE),
               mu = predict(model, type = "response"))
  return(y)
}

ffun <- function(resp) {
      gam(resp ~ s(YEAR) + D1AUG + DEPTH + s(TURBID) + s(AREA, bs = "re"),
          family = family,
          method = method,
          data = data)
}

## 1 run to produce 1 half-normal plot only
## This can take many seconds to run given the GAMM's complexity.

hnp(model,
    newclass = TRUE,
    diagfun = dfun,
    simfun = sfun,
    fitfun = ffun,
    how.many.out = TRUE,
    plot.sim = TRUE,
    paint = TRUE,
    ylab = "Deviance residuals")

## 10 "hnp" runs to get a mode and mean % of residuals outside the simulated
## envelope. This can take many minutes to run given the GAMM's level of complexity.

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

## Calculate the mgcViz score (mode and mean) based on 100 iterations, each relying
## on 100 model-based simulation, for a given model and its main explanatory variable.
## A 95% "uncertainty interval" (ui) is also calculated. This can take several 
## seconds to run.

set.seed(2025)

model <- m_WALLEYE_NB2_FINAL
predictor <- "YEAR"

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

## Concurvity and collinearity

## Check first for possible concurvity issues between s(YEAR), s(TURBIDITY), and
## s(AREA)

concurvity(m_WALLEYE_NB2_FINAL)

## Check for possible collinearity issues between the two parametric components
## of the model, i.e. D1AUG and DEPTH, with check_collinearity of "performance".
## A VIF value will also be associated to the smooth terms, but should be ignored.

check_collinearity(m_WALLEYE_NB2_FINAL)


## MODEL SELECTION -------------------------------------------------------------

## Using the compareML() function of "itsadug" shows that the final NB2 GAMM offers
## an equivalent fit to that of the initial full counterpart despite being much
## less complex, as captured by the Estimated degrees of freedom (Edf).

compareML(m_WALLEYE_NB2_FINAL, m_WALLEYE_NB2_FULL)

## Comparing these same two models under an information-theoretic approach
## indicates that the final NB2 GAMM is much more parsimonious than its full
## counterpart using the model.sel() function of "MuMIn" according to the AICc

model.sel(m_WALLEYE_NB2_FINAL, m_WALLEYE_NB2_FULL)

## The test based on the Bayes Factor (BF) in "performance" also favors the final
## model over the initial full NB2 GAMM built in "mgcv"

test_bf(m_WALLEYE_NB2_FULL, m_WALLEYE_NB2_FINAL)


## "IN-SAMPLE" PREDICTIVE PERFORMANCE ------------------------------------------

## First refit the final NB2 GAMM with REML instead of ML

m_WALLEYE_NB2_FINAL_REML <- gam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID)
                                + s(AREA, bs = "re"),
                                family = nb,
                                method = "REML",
                                data = WALLEYE)

summary(m_WALLEYE_NB2_FINAL_REML)

## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the final model relying on REML estimation

model <- m_WALLEYE_NB2_FINAL_REML

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
K <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - K) * (100 - D2))
cbind(D2, D2_adj)

## Estimate the contribution of the fixed effect YEAR and the random effect AREA
## to the estimated deviance explained (D2 above) with "gam.hp" from the model
## just defined above.

gam.hp(model)


## MODEL PREDICTIONS -----------------------------------------------------------

## Predictions with "mgcv" for the final NB2 GAMM refitted with REML (Fig. 2a)

## Create new data

nd_WALLEYE <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                         D1AUG = mean(WALLEYE$D1AUG),
                         DEPTH = mean(WALLEYE$DEPTH),
                         TURBID = mean(WALLEYE$TURBID),
                         AREA = "ALSP")

## Check the new data

headtail(nd_WALLEYE)

## Predictions

model <- m_WALLEYE_NB2_FINAL_REML

fitted <- predict(model,
                  nd_WALLEYE,
                  type = "link",
                  exclude = "s(AREA)",
                  se.fit = TRUE)

estimate <- exp(fitted$fit)
lower_ci <- exp(fitted$fit - 1.96 * fitted$se.fit)
upper_ci <- exp(fitted$fit + 1.96 * fitted$se.fit)

cbind(YEAR = nd_WALLEYE$YEAR, estimate, lower_ci, upper_ci)


## TEMPORAL AUTOCORRELATION ----------------------------------------------------

## Check for temporal autocorrelation in the final NB2 GAMM with REML ("mgcv")
## according to the test of Durbin and Watson (1950) implemented in "performance"
## which can detect such analytical issues.

check_autocorrelation(m_WALLEYE_NB2_FINAL_REML)

## Estimate the first-order autoregressive (AR1) coefficient, referred to as rho,
## with "itsadug" by refitting the final NB2 GAMM with the bam() function instead
## of gam() in "mgcv" using fast REML (fREML) estimation instead of REML. The 
## argument discrete = TRUE is also used (see van Rij et al. 2022).

m_WALLEYE_NB2_FINAL_fREML <- bam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID)
                                 + s(AREA, bs = "re"),
                                 family = nb, 
                                 method = "fREML", 
                                 discrete = TRUE, 
                                 data = WALLEYE)

summary(m_WALLEYE_NB2_FINAL_fREML)

(RHO <- start_value_rho(m_WALLEYE_NB2_FINAL_fREML, plot = TRUE))

## As RHO = 0.111, this indicates that approximately 11% of the variation
## in the response variable (N) in a given year is explained by that observed in
## the previous year. The same NB2 GAMM is refitted to now account for RHO and so
## AR1 is added to its name to identify it.

m_WALLEYE_NB2_FINAL_fREML_AR1 <- bam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID)
                                     + s(AREA, bs = "re"),
                                     family = nb,
                                     method = "fREML",
                                     discrete = TRUE,
                                     rho = RHO,
                                     data = WALLEYE)

summary(m_WALLEYE_NB2_FINAL_fREML_AR1)

## Visually inspect whether the addition of an autocorrelation matrix with RHO
## has reduced the temporal autocorrelation by looking at the ACF plot with "itsadug"

check_resid(m_WALLEYE_NB2_FINAL_fREML_AR1, select = 3)

## Use the compareML() function of "itsadug" to compare both models (with and 
## without RHO)

compareML(m_WALLEYE_NB2_FINAL_fREML, m_WALLEYE_NB2_FINAL_fREML_AR1)

## Compare both models under an information-theoretic approach based on AICc

model.sel(m_WALLEYE_NB2_FINAL_fREML, m_WALLEYE_NB2_FINAL_fREML_AR1)

## Compare both models using a test relying on the Bayes Factor (BF)

test_bf(m_WALLEYE_NB2_FINAL_fREML, m_WALLEYE_NB2_FINAL_fREML_AR1)

## Accounting for the AR1 term (using RHO) has statistically improved the model 
## fit, but when its predictions are plotted against those of the NB2 GAMM with 
## REML (i.e., not accounting for temporal autocorrelation), the differences 
## between the two regression curves are quite slight.


## REFERENCES (NOT CITED IN THE MAIN TEXT) -------------------------------------

## Consul PC (1989) Generalized Poisson distributions: properties and 
##    applications. Marcel Dekker, New York

## Efron B (1986) Double exponential families and their use in the generalized
##    linear regression. J Am Stat Assoc 81:709–721.

## Holla MS (1967) On a poisson-inverse gaussian distribution. Metrika 11:115–121.
