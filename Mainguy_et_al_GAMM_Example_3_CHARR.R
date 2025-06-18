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

file_name <- 
"https://raw.githubusercontent.com/rachel-mc/Ecological-GAMM-examples/main/CHARR.txt"
CHARR <- read.delim(file_name)


## GENERALIZED ADDITIVE MIXED-EFFECTS MODELING (GAMM) --------------------------

## Here the GAMMs will be identified with *_MGCV_* or *_GAMLSS_* given that both
## packages will be used in this last example.

## A binomial (BI) GAMM is fitted in "mgcv" with a complementary log-log (cloglog)
## link, which was previously used to model these data as Bernoulli trials (0/1)
## by Mainguy et al. (under review). A smooth function is applied to BIN50 to check
## if additional non-linearities are present in the binomial regression curve
## from modeling the same data as discrete proportions instead.

CHARR$RIVER <- as.factor(CHARR$RIVER)

m_CHARR_MGCV_BI <- gam(cbind(YES, TOTAL - YES) ~ s(BIN50) + s(RIVER, bs = "re"),
                       family = binomial(link = "cloglog"),
                       method = "ML",
                       data = CHARR)

summary(m_CHARR_MGCV_BI)

## As edf = 1 for s(BIN50), this fixed predictor is refitted as a parametric
## component instead using the same model name and now with REML estimation given
## that this model will not be compared to other candidate models.

m_CHARR_MGCV_BI_REML <- gam(cbind(YES, TOTAL - YES) ~ BIN50 + s(RIVER, bs = "re"),
                            family = binomial(link = "cloglog"),
                            method = "REML",
                            data = CHARR)

summary(m_CHARR_MGCV_BI_REML)

## Estimate the L50 using the fixed effects, where b0 is the intercept coefficient
## and b1 is the coefficient associated to BIN50 (i.e., slope). The equation used
## to estimate the L50 when a cloglog link is used is detailed in the Supplementary
## Material of Mainguy et al. (2024).

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

## Estimate the L50 using the fixed effects. The dose.p() function of "MASS" 
## cannot be used with a "gamlss" object.

b0 <- coef(m_CHARR_GAMLSS_BI)[[1]]
b1 <- coef(m_CHARR_GAMLSS_BI)[[2]]
((-b0 - 0.3665129) / b1)

## Now, four candidate GAMMs are created in "gamlss" that are using the same data,
## except that the random effect RIVER is voluntarily omitted to create overdispersion.
## The binomial (BI) distribution is therefore expected to no longer be adequate.
## The beta-binomial (BB), zero-inflated BI (ZIBI) and zero-inflated BB (ZIBB) will 
## possibly better model such overdispersed discrete proportions that are also
## potentially zero-inflated. Models are named with *_NRS_* = *No Random Structure*.

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

## Despite the fact that all four GAMMs converged, the one fitted with ZIBI yielded
## an 'infinite' SE for the estimate of the zero-inflated term using a logit link.
## This indicates modelling issues that would lead to unreliable inferences and 
## thus, this model was discarded. The lower L50 estimate for ZIBB was also less
## plausible given its related regression curve, such that it was discarded too.


## MODEL ADEQUACY --------------------------------------------------------------

## Diagnostic plots using the gam.check() function of "mgcv".

gam.check(m_CHARR_MGCV_BI_REML)

## The same diagnostic plots using the appraise() function of "gratia" with a 
## simulated, theoretically-based envelope relying on deviance residuals (default).

appraise(m_CHARR_MGCV_BI_REML, method = "simulate")

## Worm plot using the worm_plot() function of "gratia" with a simulated envelope.

worm_plot(m_CHARR_MGCV_BI_REML, method = "simulate", n_simulate = 1000)

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
    paint = TRUE,
    ylab = "Deviance residuals")

## 10 "hnp" iterations

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

## Calculate the mgcViz score (mode and mean %) for a given model based on 
## 100 iterations. A 95% "uncertainty interval" (ui) is also calculated.
## This can take several seconds to run.

set.seed(2025)

model <- m_CHARR_MGCV_BI_REML
predictor <- "BIN50"

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

## Check concurvity

concurvity(m_CHARR_MGCV_BI_REML)

## Diagnostic plots using the plot() method of "gamlss", here for the GAM fitted
## with the beta-binomial and no random structure (m_CHARR_GAMLSS_NRS_BB)

plot(m_CHARR_GAMLSS_NRS_BB)

## Worm plot using the wp() function of "gamlss", here for the GAM fitted with
## the beta-binomial (m_CHARR_GAMLSS_NRS_BB)

wp(m_CHARR_GAMLSS_NRS_BB)

## Bucket plot using the bp() function of "gamlss", here for the GAM fitted with
## the beta-binomial (m_CHARR_GAMLSS_NRS_BB)

bp(m_CHARR_GAMLSS_NRS_BB)

## Adequacy assessment using "hnp" for the GAM fitted with the beta-binomial in
## in "gamlss" when relying on a single "hnp" iteration. No helper functions are
## needed, with 19 simulations and a confidence level of 100% being used to allow
## to obtain an half-normal plot for adequacy assessment purposes (see Moral et
## al. 2017). Given that n < 30 (see Mainguy and Moral 2021), only a visual
## assessment should be performed given that each residual found outside of the 
## envelope in this example (n = 26) represent almost 4% of the sample. If the
## assessment fails, rerun the same code as the fitting procedure of "gamlss"
## when using the beta-binomial distribution with such small sample size is likely
## to cause analytical issues, i.e. the reason for relying on less simulations.

model <- m_CHARR_GAMLSS_NRS_BB

hnp(model,
    how.many.out = TRUE,
    plot = TRUE,
    paint = TRUE,
    sim = 19,
    conf = 1,
    ylab = "Deviance residuals")

## The binomial GAM fitted in "gamlss", which cannot currently be assessed with
## "hnp", is fitted in "mgcv" instead to then allow its assessment with this
## package using the same parametrization as above (sim = 19 and conf = 1) for
## comparison purposes of diagnostic half-normal plots on a same basis. Note that
## the random RIVER effect has been removed and that *_NRS_* has been added.

CHARR$RIVER <- as.factor(CHARR$RIVER)

m_CHARR_MGCV_NRS_BI <- gam(cbind(YES, TOTAL - YES) ~ s(BIN50),
                           family = binomial(link = "cloglog"),
                           method = "ML",
                           data = CHARR)

model <- m_CHARR_MGCV_NRS_BI
family <- binomial(link = "cloglog")
method <- "ML"
data <- CHARR
size <- CHARR$TOTAL

sfun <- function(n, obj) {
  p_hat <- predict(obj, type = "response")
  y <- rbinom(n = n, size = size, prob = p_hat)
  return(y)
}

dfun <- function(obj) resid(obj, type = "deviance")

ffun <- function(new_response) {
  gam(cbind(new_response, TOTAL - new_response) ~ BIN50,
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
    paint = TRUE,
    ylab = "Deviance residuals")

## The related half-normal plot suggest weak overdispersion. Because this simple
## binomial GAM has been fitted in "mgcv", this also allow to test for likely
## overdispersion with "DHARMa" which can also test for possible zero-inflation.
## Overdispersion is specifically assessed (as opposed to equidispersion) by also
## using the argument alternative = "greater" with the testDispersion() function.

simulationOuput <- simulateResiduals(fittedModel=m_CHARR_MGCV_NRS_BI)

testDispersion(simulationOuput, alternative = "greater")

## Using the same simulated randomized quantile residuals that were produced above,
## the testZeroInflation() function is used to assess a possible excess of observed
## zeros relative to those being predicted, which appears unlikely given that the
## zero-inflated binomial and zero-inflated beta-binomial were both discarded.

testZeroInflation(simulationOuput)

## With detected overdispersion, but no evidence for zero-inflation, the binomial
## GAM is not sufficiently adequate to model these slightly overdispersed discrete
## proportions, such that the adequate beta-binomial that handles the extra-binomial
## variation is the logical choice.


## MODEL SELECTION -------------------------------------------------------------

## The beta-binomial GAM was the best (and only) choice when the random RIVER 
## effect was removed, such that no step of model selection is required. Below,
## we show for demonstration purposes only, how the binomial GAMM making use of
## the random RIVER effect could be a better choice than the beta-binomial GAM.

## An information-theoretic approach based on AICc.

model.sel(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_BI)

## The test based on the Bayes Factor (BF) in "performance". 

test_bf(m_CHARR_GAMLSS_NRS_BB, m_CHARR_GAMLSS_BI)

## For both model selection approaches, the binomial GAMM using the random RIVER
## effect is identified as the one providing the best relative fit to the
## modeled (overdispersed) discrete proportions when compared to the
## beta-binomial GAM.


## "IN-SAMPLE" PREDICTIVE PERFORMANCE ------------------------------------------

## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the BI GAMM in "mgcv" that relies on REML estimation.

model <- m_CHARR_MGCV_BI_REML 

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
K <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - K) * (100 - D2))
cbind(D2, D2_adj)


## MODEL PREDICTIONS -----------------------------------------------------------

## New data for "mgcv"

nd_CHARR <- data.frame(BIN50 = seq(175, 675, by = 1), RIVER = "VOLTZ")

## Predictions with "mgcv" (Fig. 3a).

model <- m_CHARR_MGCV_BI_REML

fitted <- predict(model,
                  nd_CHARR,
                  type = "link",
                  exclude = "s(RIVER)",
                  se.fit = TRUE)

pred <- 1 - exp(-exp(fitted$fit))
lower <- 1 - exp(-exp(fitted$fit - 1.96 * fitted$se.fit))
upper <- 1 - exp(-exp(fitted$fit + 1.96 * fitted$se.fit))

bin50 <- nd_CHARR$BIN50
cbind(bin50, pred, lower, upper)

## Predictions with "gamlss" (Fig. 3a).

## New data for "gamlss"

nd_CHARR_ALL <- data.frame(BIN50 = seq(175, 675, by = 1))

nd_CHARR_AIPPARUSIK <- data.frame(BIN50 = seq(175, 675, by = 1), RIVER = "AIPPARUSIK")

nd_CHARR_TASIALLUJUAK <- data.frame(BIN50 = seq(175, 675, by = 1), RIVER = "TASIALLUJUAK")

nd_CHARR_VOLTZ <- data.frame(BIN50 = seq(175, 675, by = 1), RIVER = "VOLTZ")

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

## Predictions for "gamlss" models without the random RIVER effect (Fig. 3b).
## Predicted values are obtained from the observed data so that se.fit = TRUE
## can be used. The ZIBI and ZIBB GAM can have their predictions generated in the
## same way as the beta-binomial (BB) GAM below, but here were both discarded.

model <- m_CHARR_GAMLSS_NRS_BB

fitted <- predict(model,
                  what = "mu",
                  type = "link",
                  se.fit = TRUE)

estimate <- 1 - exp(-exp(fitted$fit))
lower_ci <- 1 - exp(-exp(fitted$fit - 1.96 * fitted$se.fit))
upper_ci <- 1 - exp(-exp(fitted$fit + 1.96 * fitted$se.fit))

cbind(CHARR, estimate, lower_ci, upper_ci)

## Note that the predicted values are the same for a given BIN50, i.e. regardless
## of the RIVER considered. Additional data manipulations are required to extract
## only one estimated probability with its related uncertainty to then produce the
## Fig. 3 presented in the main text (not shown). Using "mgcv" offers a clear
## advantage over "gamlss" to plot predictions.
