## Required packages:

library(DHARMa)
library(FSA)
library(gam.hp)
library(gratia)
library(hnp)
library(itsadug)
library(lmtest)
library(mgcv)
library(mgcViz)
library(moments)
library(MuMIn)
library(nortest)
library(nlme)
library(nlraa)
library(performance)
library(tidyverse)
library(visreg)

################################################################################
#### EXAMPLE 1: LENGTH-AT-AGE RELATIONSHIP IN AUPALUK FEMALE LAKE TROUT ########
################################################################################

## Load the TROUT dataset as an example with a continuous response variable 

file_name <- 
"https://raw.githubusercontent.com/rachel-mc/Ecological-GAMM-examples/main/TROUT.txt"

TROUT <- read.delim(file_name)


## NON-LINEAR MIXED-EFFECTS VON BERTALANFFY GROWTH MODEL -----------------------

## Create a function for the von Bertalanffy (VB) growth equation (see Ogle 2016
## and Ogle et al. 2025 for a description of the different growth parameters)

VB <- function(x, Linf, K, t0) Linf * (1 - exp(-K * (x - t0)))

## Estimate the starting values (sv) to be used in the nonlinear modeling with
## findGrowthStarts() function of "FSA". Asking for a plot to be printed allows
## to ensure that the global fit of the von Bertalanffy growth function to the
## length-at-age data seems to provide a reasonably good fit.

(sv_TROUT_VB <- findGrowthStarts(TL ~ AGE, 
                                 data = TROUT,
                                 type = "von Bertalanffy",
                                 plot = TRUE))

## Estimate the growth parameters for each LAKE (LAKE_PARMS) with the nlsList()
## function of "nlme". There are no missing values (NA's) in the dataset, yet
## the argument na.action = na.omit() is used to ensure that the nLslist() function
## will properly work.

LAKE_PARMS <- nlsList(TL ~ VB(AGE, Linf, K, t0) | LAKE,
                      start = c(Linf = sv_TROUT_VB[[1]],
                                K = sv_TROUT_VB[[2]],
                                t0 = sv_TROUT_VB[[3]]),
                      data = TROUT,
                      na.action = na.omit)

## Fit the non-linear mixed-effects von Bertalanffy growth model with the nlme()
## function of "nlme"

m_TROUT_VB <- nlme(TL ~ VB(AGE, Linf, K, t0),
                        fixed = Linf + K + t0 ~ 1,
                        random = Linf + K + t0 ~ 1 | LAKE,
                        start = fixef(LAKE_PARMS),
                        data = TROUT)

summary(m_TROUT_VB)

## The 95% confidence intervals for the three estimated growth parameters can be
## obtained as follows, if needed:

intervals(m_TROUT_VB, which = "fixed")

## The same approach was applied using instead the Gompertz and logistic growth
## functions (see Ogle et al. 2025 for equation-specific details), with the
## Gompertz growth model failing at converging, whereas the logistic one converged
## correctly but was found to exhibit a delta AIC of 4.32 relative to that of the
## von Bertalanffy growth model, such that the latter was preferred (code and
## analyses not shown). Finally, TL was predicted according to AGE from the 
## nonlinear mixed-effects von Bertalanffy growth model (Fig. 1a) with the code
## below. The related bootstrapped 95% CI is estimated with the boot_nlme() 
## function of "nlme". This step can takes several minutes to run. For details,
## see: https://cran.r-project.org/web/packages/nlraa/vignettes/Confidence-Bands.html

nd <- expand.grid(AGE = seq(1, 50, 0.1), LAKE = c("BRULE", "REDDOG", "VOLTZ"))

PRED_TL <- function(x) predict(x, newdata = nd, level = 0)

TROUT_VB_BOOT <- boot_nlme(m_TROUT_VB, PRED_TL)

TROUT_VB_BOOT_SUMMARY <- cbind(nd[,-4], summary_simulate(t(na.omit(TROUT_VB_BOOT$t))))

TROUT_VB_FINAL <- aggregate(cbind(Estimate, Est.Error, Q2.5, Q97.5) ~ AGE,
                                  data = TROUT_VB_BOOT_SUMMARY, FUN = "mean")

## Calculate the marginal R2 and conditional R2 with the R2M() function of the
## "nlraa" package (Miguez 2023)

R2M(m_TROUT_VB)


## GENERALIZED ADDITIVE MIXED-EFFECTS MODELING (GAMM) --------------------------

## Fitting the Gaussian distribution, as done for the nonlinear mixed-effects von 
## Bertalanffy growth model, to a GAMM with the default parametrization for s(),
## such that it is identified with *_TP* from using thin-plate regression splines.
## Note that s(AGE, bs = "tp", k = 10) is equivalent to just using s(AGE), whereas
## family = gaussian(link = "identity") is the same as using family = gaussian.

TROUT$LAKE <- as.factor(TROUT$LAKE)

m_TROUT_GAUSSIAN_TP <- gam(TL ~ s(AGE, bs = "tp", k = 10) + s(LAKE, bs = "re"),
                           family = gaussian(link = "identity"),
                           method = "ML",
                           data = TROUT)

summary(m_TROUT_GAUSSIAN_TP)

## Fitting an alternative Gaussian GAMM with an adaptive smooth *_AD* again using
## the ML method in "mgcv". The default k = 10 is increased to 20 to further penalize
## against wiggliness. Using k = 20 is possible given that the TROUT dataset contains
## 37 distinct AGE values.

m_TROUT_GAUSSIAN_AD <- gam(TL ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
                           family = gaussian(link = "identity"),
                           method = "ML",
                           data = TROUT)

summary(m_TROUT_GAUSSIAN_AD)

## Fitting the Gamma distribution with a log link (Forbes et al. 2010) instead as
## just done above with the Gaussian distribution. Note that using family = Gamma
## only will use an inverse link function as a default setting, such that explicitly
## indicating family = Gamma(link = "log") is required for fitting a log link 
## function instead, which is intended here. An adaptive smooth *_AD* is used again
## here given the bioligically unrealistic wiggling at older ages that thin-plate
## regression splines generate due to the sparsely distributed observations.

m_TROUT_GAMMA_AD <- gam(TL ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
                        family = Gamma(link = "log"),
                        method = "ML",
                        data = TROUT)

summary(m_TROUT_GAMMA_AD)


## MODEL ADEQUACY --------------------------------------------------------------

## Diagnostic plots with the gam.check() function of "mgcv" for the Gaussian
## GAMM with an adaptive smooth. Similar diagnostic plots are produced for the
## gamma GAMM, with homoscedasticity for instance only being applicable to the
## Gaussian case (which is not apparently respected from a visual inspection).

gam.check(m_TROUT_GAUSSIAN_AD)

## The same diagnostic plots with a simulated envelope relying on deviance 
## residuals (default) using the appraise() function of "gratia", again for the
## Gaussian GAMM with an adaptive smooth as an example.

appraise(m_TROUT_GAUSSIAN_AD, method = "simulate")

## Worm plot with a simulated envelope using the worm_plot() function of "gratia".
## The default coverage level for the reference band (level = 0.90) is used but 
## the default number of simulations used (50) is increased to 1000. Note that
## the worm plot and related reference band will slightly change at each
## iteration since simulations are used, similarly as for "hnp" (see below).

worm_plot(m_TROUT_GAUSSIAN_AD, method = "simulate", n_simulate = 1000)

## Looking at a series of worm plots for the alternative gamma GAMM with an
## adaptive smooth suggests a similar fit and possibly problematic nonrandom 
## pattern in the deviance residuals near zero for both distributions despite
## most being found in the confidence band.

worm_plot(m_TROUT_GAMMA_AD, method = "simulate", n_simulate = 1000)

## Using different normality tests: (i) Anderson and Darling (1954) with the 
## "nortest" package (Gross and Ligges 2015), (ii) Shapiro and Wilks (1965) with
## the default "stats" package, and (iii) Jarque and Bera (1980) with the "moments"
## package (Komsta and Novomestky 2022) for the Gaussian GAMM. These are only used
## to approximately assess whether the behavior of the deviance residuals in the 
## previous diagnostic plots appears to sufficiently correspond to a normal curve,
## which seems to be the case in most instances.

dev_resid_GAUSSIAN_AD <- residuals(m_TROUT_GAUSSIAN_AD, type = "deviance")

ad.test(dev_resid_GAUSSIAN_AD)

shapiro.test(dev_resid_GAUSSIAN_AD)

jarque.test(dev_resid_GAUSSIAN_AD)

## Using the testUniformity() function of "DHARMa" to assess whether the simulated
## quantile residuals from the Gaussian GAMM with an adaptive smooth are sufficiently
## uniformly distribution relative to the theoretical quantiles of the normal
## distribution as a more adequate test. Model-based simulated residuals are first
## produced and then be used for this specific test.

simulationOutput <- simulateResiduals(fittedModel = m_TROUT_GAUSSIAN_AD)

testUniformity(simulationOutput)

## The same simulated residuals as above (i.e., simulationOutput) are used to
## now assess whether they are sufficiently homogeneously distributed along the
## rank-transformed model predictions. Note that in this case, homogeneity of 
## variance is considered sufficient from relying on the test statistics, but
## an uneven pattern is nonetheless visually detectable in the associated 
## diagnostic plot given the low TL variation in young female lake trout that
## we can once more visually detect in the diagnostic plot produced by DHARMa.

testQuantiles(simulationOutput)

## An approximate assessment for possible heteroscedasticity for the Gaussian GAMM
## with an adaptive smooth according to the test of Harrison and McCabe (1979)
## which is performed with the "lmtest" package (Zeiles and Hothorn 2002).

hmctest(m_TROUT_GAUSSIAN_AD, point = 0.5, simulate.p = TRUE, nsim = 1000)

## Diagnostic half-normal plot of the deviance residuals using the hnp() 
## function of "hnp" and a newly-made available sfun() helper function for the
## Gaussian distribution. The right-hand part following "resp ~" of the ffun helper
## function should match that of the Gaussian GAMM built in "mgcv" that is being
## assessed for adequacy. This step can take several seconds to run.

model <- m_TROUT_GAUSSIAN_AD 
family <- gaussian(link = "identity")
method <- "ML"  
data <- TROUT

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n, obj) {
  y <- rnorm(nrow(data),
             mean = predict(model, type = "response"),
             sd = sqrt(model$sig2))
  return(y)
}

ffun <- function(resp) {
  gam(resp ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
      family = family,
      method = method,
      data = data)
}

hnp(model,
    newclass = TRUE, 
    diagfun = dfun, 
    simfun = sfun, 
    fitfun = ffun, 
    how.many.out = TRUE, 
    plot = TRUE, 
    paint = TRUE,
    ylab = "Deviance residuals")

## Performing 10 "hnp" iterations to obtain a mode and a mean percentage of 
## residuals outside the simulated envelope for the same Gaussian GAMM with an
## adaptive smooth. A seed is only used for result reproducibility.

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

hnp_summary <- sapply(hnp_obj, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

## Mode percentage of residuals outside the simulated envelope

round(return_max(hnp_summary), 2)

## Mean percentage of residuals outside the simulated envelope with other
## descriptive statistics using the Summarize() function of "FSA".

Summarize(hnp_summary)

## Diagnostic half-normal plot of the deviance residuals using the hnp() 
## function of "hnp" and a newly-made available sfun() helper function for the
## gamma distribution with a log link that uses an adaptive *_AD* spline.

model <- m_TROUT_GAMMA_AD
family = Gamma(link = "log")
method = "ML"
data = TROUT

dfun <- function(obj) resid(obj, type = "deviance")

sfun <- function(n, obj) simulate(obj)[,1]

ffun <- function(resp) gam(resp ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
                           family = family,
                           method = method,
                           data = data)

hnp(model,
    newclass = TRUE,
    diagfun = dfun,
    simfun = sfun,
    fitfun = ffun,
    how.many.out = TRUE,
    plot = TRUE,
    paint = TRUE,
    ylab = "Deviance residuals")

## Perform 10 consecutive "hnp" runs to obtain a mode and a mean percentage of 
## the points outside the simulated envelope from the same gamma GAMM being 
## assessed for adequacy. The set.seed() function is solely used for reproducibility. 
## This step can take several minutes to run

set.seed(2025)

hfun <- list()
for(i in 1:10) {
  hfun[[i]] <- hnp(model,
                   newclass = TRUE,
                   diagfun = dfun,
                   simfun = sfun,
                   fitfun = ffun,
                   how.many.out = TRUE,
                   plot.sim = FALSE)
}

hnp_summary <- sapply(hfun, function(x) x$out / x$total * 100) 

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)

## Diagnostic plots using the getviz(), check1D(), and l_gridCheck1D() functions 
## of "mgcViz" with 100 model-derived simulations (nsim) for an "mgcv" object, here
## for the gamma GAMM with an adaptive smooth. Only the model and predictor must
## be defined for the code to work.

model <- m_TROUT_GAMMA_AD
predictor <- "AGE"

set.seed(2025)

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

## Visualize the mgcViz diagnostic plot

diag_plot

## Obtain the percentage of the binned residual means found within the
## simulated default 80% confidence limits shown on the mgcViz diagnostic plot

mgcViz_perc

## Calculate the mgcViz score (mean) for a given model based on 100 iterations,
## each relying on 100 model-based simulations. An associated 95% "uncertainty 
## interval" (ui) is also calculated. Here again for the gamma GAMM with an
## adaptive smooth. This can take several seconds to run.

model <- m_TROUT_GAMMA_AD
predictor <- "AGE"

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

## Check for potential concurvity issues. Smooth terms should have an estimate
## index < 0.5 (Kovács, 2024)

concurvity(m_TROUT_GAMMA_AD)


## MODEL SELECTION -------------------------------------------------------------

## Compare Maximum Likelihood scores while accounting for model complexity (Edf) 
## for the two candidates GAMMs using the compareML() function of "itsadug". The
## GAMMs being compared must share the same random structure.

compareML(m_TROUT_GAUSSIAN_AD, m_TROUT_GAMMA_AD)

## Compare the two candidate gamma GAMMs fitted in "mgcv" based on AICc under
## an information-theoretic approach using the model.sel() function of "MuMIn".

model.sel(m_TROUT_GAUSSIAN_AD, m_TROUT_GAMMA_AD)

## Compare the candidate gamma GAMMs fitted in "mgcv" according to a test based
## on the Bayes Factor (BF) using the test_bf() function of "performance". Note
## that the first model listed is used as the reference level.

test_bf(m_TROUT_GAUSSIAN_AD, m_TROUT_GAMMA_AD)

## Compare the non-linear mixed-effects von Bertalanffy growth model with the
## "mgcv" Gamma GAMM using an adaptive smooth using the AIC() function.

AIC(m_TROUT_VB, m_TROUT_GAMMA_AD)

## Refit the gamma GAMM using an adaptive smooth with REML estimation in "mgcv"

m_TROUT_GAMMA_AD_REML <- gam(TL ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
                             family = Gamma(link = "log"),
                             method = "REML",
                             data = TROUT)

summary(m_TROUT_GAMMA_AD_REML)


## "IN-SAMPLE" PREDICTIVE PERFORMANCE ------------------------------------------

## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the preferred gamma GAMM using an adaptive smooth relying on REML estimation.
## Note here that K (i.e., number of parameters) corresponds to the sum of df 
## and edf as estimated by "mgcv" from using the logLik() function.

model <- m_TROUT_GAMMA_AD_REML

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
K <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - K) * (100 - D2))
cbind(D2, D2_adj)

## Estimate the contribution of the fixed AGE effect and the random LAKE effect
## to the estimated D2 of the gamma GAMM using REML with "gam.hp"

gam.hp(m_TROUT_GAMMA_AD_REML)


## MODEL PREDICTIONS -----------------------------------------------------------

## Visualization of the partial residuals of the fixed AGE effect and the random
## LAKE effect (i.e., factor) at the linear predictor scale with "visreg" for 
## the gamma GAMM with an adaptive smooth using REML. The predictions can also
## be visualized at the response scale, with both predictors being sequentially
## shown each time.

visreg(m_TROUT_GAMMA_AD_REML, scale = "linear")

visreg(m_TROUT_GAMMA_AD_REML, scale = "response")

## Predictions from the initial Gaussian GAMM using thin-plate regression splines.
## The argument exclude = "s(LAKE)" is used to provide predictions for the fixed
## effects only. As such, providing only one of the different considered levels 
## for the random effect is sufficient in the new data upon which the predictions
## will be made.

model <- m_TROUT_GAUSSIAN_TP

nd_TROUT <- data.frame(AGE = seq(1, 50, by = 0.1), LAKE = "BRULE")

fitted <- predict(model,
                  nd_TROUT,
                  type = "link",
                  exclude = "s(LAKE)",
                  se.fit = TRUE)

estimate <- fitted$fit
lower_ci <- fitted$fit - 1.96 * fitted$se.fit
upper_ci <- fitted$fit + 1.96 * fitted$se.fit

age <- nd_TROUT$AGE
cbind(age, estimate, lower_ci, upper_ci)

## Predictions from the best-retained gamma GAMM with REML estimation (Fig. 1c).
## Note that the predictions are at the log link scale, such that exponentiation
## of the predictions is required for their visualization at the response scale.

model <- m_TROUT_GAMMA_AD_REML

nd_TROUT <- data.frame(AGE = seq(1, 50, by = 0.1),
                       LAKE = "BRULE")

fitted <- predict(model,
                  nd_TROUT,
                  type = "link",
                  exclude = "s(LAKE)",
                  se.fit = TRUE)

estimate <- exp(fitted$fit)
lower_ci <- exp(fitted$fit - 1.96 * fitted$se.fit)
upper_ci <- exp(fitted$fit + 1.96 * fitted$se.fit)

age <- nd_TROUT$AGE
cbind(age, estimate, lower_ci, upper_ci)


## REFERENCES (NOT CITED IN THE MAIN TEXT) -------------------------------------

## Anderson TW, Darling DA (1954) A test of goodness of fit. J Am Stat Assoc 49:
##    765–769

## Forbes C, Evans M, Hastings N et al. (2010) Gamma distribution. In: Forbes C,
##    Evans M, Hastings N et al (eds) Statistical distributions, 4th edn. John 
##    Wiley & Sons, New Jersey, pp 109–113.

## Gross J, Ligges U (2015) nortest: Tests for normality. R package version 1.0-4.
##    https://CRAN.R-project.org/package=nortest

## Jarque CM, Bera AK (1980) Efficient tests for normality, homoscedasticity and
##    serial independence of regression residuals. Econ Lett 6:255–259

## Komsta L, Novomestky F (2022) moments: Moments, cumulants, skewness, kurtosis
##    and related tests. R package version 0.14.1.
##    https://CRAN.R-project.org/package=moments

## Miguez F (2023) nlraa: Nonlinear regression for agricultural applications.
##    R package version 1.9.7. https://CRAN.R-project.org/package=nlraa

## Shapiro SS, Wilk MB (1965) An analysis of variance test for normality (complete
##    samples). Biometrika 52:591-611

## Zeileis A, Hothorn T (2002) Diagnostic checking in regression relationships. R
##    News 2:7-10
