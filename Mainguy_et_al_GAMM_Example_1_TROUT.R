## Required packages:

library(FSA)
library(gam.hp)
library(gamlss)
library(gratia)
library(hnp)
library(itsadug)
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

file_name <- "https://raw.githubusercontent.com/rachel-mc/Ecological-GAMM-examples/main/TROUT.txt"
TROUT <- read.delim(file_name)


# Non-linear mixed-effects logistic growth model --------------------------


## Create a function for the logistic (LG) growth equation

LG <- function(x, Linf, K, t0) Linf / (1 + exp(-K * (x - t0)))

## Estimate the initial parameters (IP) to be used as starting values with the
## vbStarts() function of "FSA". Although the growth parameter K for the logistic
## equation differs from that of the von Bertalanffy and Gompertz equations, its
## estimation with the vbStarts() function suffices here as a starting value

IP_TROUT <- vbStarts(formula = TL ~ AGE, data = TROUT)

## Estimate the growth parameters for each LAKE (LAKE_PARMS) with the nlsList()
## function of "nlme". There are no missing values (NAs) in the dataset, yet
## the argument na.action = na.omit() is used to ensure that the nLslist() function
## will properly work

LAKE_PARMS <- nlsList(TL ~ LG(AGE, Linf, K, t0) | LAKE,
                      data = TROUT,
                      start = c(Linf = IP_TROUT$Linf,
                                K = IP_TROUT$K,
                                t0 = IP_TROUT$t0),
                      na.action = na.omit)

## Fit the (global) non-linear mixed-effects logistic growth model with the
## nlme() function of "nlme"

m_TROUT_NLME_LOGISTIC <- nlme(TL ~ LG(AGE, Linf, K, t0),
                              fixed = Linf + K + t0 ~ 1,
                              random = Linf + K + t0 ~ 1 | LAKE,
                              start = fixef(LAKE_PARMS),
                              data = TROUT)

summary(m_TROUT_NLME_LOGISTIC)

## The 95% confidence intervals for the three estimated growth parameters can be
## obtained as follows:

intervals(m_TROUT_NLME_LOGISTIC, which = "fixed")

## Predict total length (TL) according to AGE from the logistic growth model (Fig. 1a). 
## The related bootstrapped 95% CI is estimated with the boot_nlme() function of "nlme".
## See https://cran.r-project.org/web/packages/nlraa/vignettes/Confidence-Bands.html
## This step can takes several minutes to run.

nd <- expand.grid(AGE = seq(1, 50, 0.1),
                  LAKE = c("BRULE", "REDDOG", "VOLTZ"))

PRED_TL <- function(x) predict(x, newdata = nd, level = 0)

m_TROUT_NLME_LOGISTIC_BOOT <- boot_nlme(m_TROUT_NLME_LOGISTIC, PRED_TL)

m_TROUT_NLME_LOGISTIC_BOOT_SUMMARY <- cbind(nd[,-4],
                                            summary_simulate(t(na.omit(m_TROUT_NLME_LOGISTIC_BOOT$t))))

m_TROUT_NLME_LOGISTIC_FINAL <- aggregate(cbind(Estimate, Est.Error, Q2.5, Q97.5) ~ AGE,
                                         data = m_TROUT_NLME_LOGISTIC_BOOT_SUMMARY,
                                         FUN = mean)

## Calculate the marginal and conditional R2 with the R2M() function of the
## "nlraa" package (Miguez 2023)

R2M(m_TROUT_NLME_LOGISTIC)


# Generalized Additive Mixed-Effects modelling (GAMM) ----------------------------


## Fitting the Gamma distribution with a log link (Forbes et al. 2010) using the 
## gam() function of the "mgcv" framework according to the default parameterization
## for s(), i.e., using thin plate regression splines (tprs) and 10 basis dimensions.
## Note here that s(AGE) is equivalent to using s(AGE, bs = "tp", k = 10).
## The Maximum Likelihood (ML) estimation method is used. 
## The related predictions are presented in Fig. 1b. 
## IMPORTANT: fitting family = Gamma only will use of an inverse link function, such that
## explicitly indicating family = Gamma(link = "log") is required for a log link function.

TROUT$LAKE <- as.factor(TROUT$LAKE)

m_TROUT_MGCV_GAMMA <- gam(TL ~ s(AGE) + s(LAKE, bs = "re"),
                          family = Gamma(link = "log"),
                          method = "ML",
                          data = TROUT)

summary(m_TROUT_MGCV_GAMMA)

## Fitting an alternative Gamma GAMM with an adaptive (ad) smooth and using the ML method 
## in "mgcv". The default k = 10 is increased to 20 to further penalise against 
## wiggliness. Using k = 20 is possible given that the TROUT dataset contains 37
## distinct AGE values. The related predictions are presented in Fig. 1c.

m_TROUT_MGCV_GAMMA_AD <- gam(TL ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
                             family = Gamma(link = "log"),
                             method = "ML",
                             data = TROUT)

summary(m_TROUT_MGCV_GAMMA_AD)

## Fitting a Gamma GAMM using the recommended P-spline, pb(), in "gamlss".
## The related predictions are presented in Fig. 1b.

m_TROUT_GAMLSS_GAMMA <- gamlss(TL ~ pb(AGE) + random(LAKE),
                               family = GA,
                               data = TROUT)

summary(m_TROUT_GAMLSS_GAMMA)

## Fitting a Gamma GAMM using instead a monotonic P-spline, pbm(), in "gamlss".
## The related predictions are presented in Fig. 1c.

m_TROUT_GAMLSS_GAMMA_PBM <- gamlss(TL ~ pbm(AGE) + random(LAKE),
                                   family = GA,
                                   data = TROUT)

summary(m_TROUT_GAMLSS_GAMMA_PBM)


# Model Adequacy ----------------------------------------------------------


## Diagnostic plots with the gam.check() function of "mgcv"

gam.check(m_TROUT_MGCV_GAMMA_AD)

## The same diagnostic plots with a simulated, theoretically-based envelope
## relying on deviance residuals (default) using the appraise() function of "gratia"

appraise(m_TROUT_MGCV_GAMMA_AD, 
         method = "simulate")

## Worm plot with a simulated envelope using the worm_plot() function of "gratia".
## The default coverage level for the reference band (level = 0.90) is used but 
## the default number of simulations used (50) is increased to 100. Note that
## the worm plot and related reference band will slightly change at each
## iteration since simulations are used.

worm_plot(m_TROUT_MGCV_GAMMA_AD,
          method = "simulate",
          n_simulate = 100)

## Although residual normality is not one of the Gamma distributional assumptions,
## the deviance residuals are extracted and checked for demonstration purposes only, i.e., as
## if the Gaussian distribution was used instead. 
## The normality test of Anderson and Darling (1954) with the "nortest" package (Gross and Ligges 2015),
## that of Shapiro and Wilks (1965) with the default "stats" package, and that of 
## Jarque and Bera (1980) with the "moments" package (Komsta and Novomestky 2022)
## are all applied

DEV_RESID_MGCV <- residuals(m_TROUT_MGCV_GAMMA_AD,
                            type = "deviance")

ad.test(DEV_RESID_MGCV)

shapiro.test(DEV_RESID_MGCV)

jarque.test(DEV_RESID_MGCV)

## Diagnostic half-normal plot of the deviance residuals using the hnp() 
## function of "hnp" and a newly-made available sfun helper function for the
## Gamma distribution. The right-hand part following "resp ~" of the ffun helper
## function should match that of the Gamma GAMM built in "mgcv" that is being
## assessed for adequacy. This step can take several seconds to run.

model <- m_TROUT_MGCV_GAMMA_AD
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
    paint = TRUE)

## Perform 10 consecutive hnp runs to obtain a mode and a mean percentage of the points
## outside the simulated envelope from the same Gamma GAMM (defined model above) 
## being assessed for adequacy. 
## The set.seed() function is solely used for reproducibility. 
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

## Find the asymptotic value (i.e., mode) of a density curve describing these
## 10 "hnp" percentage values

return_max <- function(numvec) {
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
  }

round(return_max(hnp_summary), 2)

## Descriptive statistics (including the mean) of the percentages of residuals
## found outside the simulated envelope for 10 "hnp" runs 

Summarize(hnp_summary)

## Diagnostic half-normal plot of the (normalised) randomised quantile residuals 
## using the hnp() function of "hnp" and a newly-made available sfun helper 
## function for a Gamma model fitted in "gamlss" (with monotonic P-splines).
## The right-hand part following "resp ~" of the ffun helper function should 
## match that of the Gamma GAMM fitted in "gamlss" that is being assessed for
## adequacy. This step can take several seconds to run.

model <- m_TROUT_GAMLSS_GAMMA_PBM
family = GA
data = TROUT

dfun <- function(obj) resid(obj) 

sfun <- function(n, obj) { 
  mu <- obj$mu.fv 
  sig <- obj$sigma.fv 
  rGA(n, mu = mu, sigma = sig) 
  } 

ffun <- function(resp) gamlss(resp ~ pbm(AGE) + random(LAKE),
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

## Perform 10 consecutive hnp runs to obtain a modal and mean percentage.
## This step can takes several minutes to run

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

hnp_summary <- sapply(hfun, function(x) x$out/x$total*100) 

round(return_max(hnp_summary), 2)

Summarize(hnp_summary)

## Diagnostic plots using the getviz(), check1D(), and l_gridCheck1D() functions 
## of "mgcViz" with 100 model-derived simulations (nsim) for an "mgcv" object

model <- m_TROUT_MGCV_GAMMA_AD
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

## Calculate the mgcViz score (mean) for a given model based on 
## 100 iterations. An associated 95% "uncertainty interval" is also calculated.
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

summary_mgcViz <- replicate(100, viz_fun(model = m_TROUT_MGCV_GAMMA_AD, predictor = "AGE"))

mgcViz_score <- mean(summary_mgcViz)
lower_ui <- (quantile(summary_mgcViz, probs = 0.025))[[1]]
upper_ui <- (quantile(summary_mgcViz, probs = 0.975))[[1]]

cbind(mgcViz_score, lower_ui, upper_ui)

## Check for potential concurvity issues
## Smooth terms should have an estimate index < 0.5 (Kovács, 2024)

concurvity(m_TROUT_MGCV_GAMMA_AD)

## Diagnostic plots using the plot() method of "gamlss". As this assessment
## is based on quantile residuals, these should approximately follow a Normal
## distribution and exhibit even variation along the fitted values

plot(m_TROUT_GAMLSS_GAMMA_PBM)

## Worm plot using the wp() function of "gamlss"

wp(m_TROUT_GAMLSS_GAMMA_PBM)

## Bucket plot using the bp() function of "gamlss"

bp(m_TROUT_GAMLSS_GAMMA_PBM)

## Perform the same normality tests as for "mgcv" above

dev_resid_gamlss <- resid(m_TROUT_GAMLSS_GAMMA_PBM)

ad.test(dev_resid_gamlss)

shapiro.test(dev_resid_gamlss)

jarque.test(dev_resid_gamlss)


# Model Selection ---------------------------------------------------------


## Compare Maximum Likelihood scores while accounting for model complexity (Edf) 
## for "mgcv" GAMMs using the compareML() function of "itsadug".
## GAMMs being compared must share the same random structure.

compareML(m_TROUT_MGCV_GAMMA, m_TROUT_MGCV_GAMMA_AD)

## Compare the two Gamma GAMMs fit in either "mgcv" or "gamlss" using AICc scores
## from the model.sel() function of "MuMIn". As the methods used (ML vs. RS) are
## not the same, no comparison is made between the two GAM-oriented packages.

model.sel(m_TROUT_MGCV_GAMMA, m_TROUT_MGCV_GAMMA_AD)

model.sel(m_TROUT_GAMLSS_GAMMA, m_TROUT_GAMLSS_GAMMA_PBM)

## Compare the two Gamma GAMMs obtained in "mgcv" or "gamlss" according to
## a test based on the Bayes Factor (BF) using the test_bf() function of
## "performance"

test_bf(m_TROUT_MGCV_GAMMA, m_TROUT_MGCV_GAMMA_AD)

test_bf(m_TROUT_GAMLSS_GAMMA, m_TROUT_GAMLSS_GAMMA_PBM)

## Compare the non-linear mixed-effects Logistic growth model with the retained 
## "mgcv" Gamma GAMM with an adaptive smoother using the AIC() and BIC() functions.
## Both of these models used ML estimation, as opposed to "gamlss" that uses the
## RS algorithm to maximise the penalised likelihood.

AIC(m_TROUT_NLME_LOGISTIC, m_TROUT_MGCV_GAMMA_AD)

BIC(m_TROUT_NLME_LOGISTIC, m_TROUT_MGCV_GAMMA_AD)

## Refit the Gamma GAMM with an adaptive smooth and REML estimation in "mgcv"

m_TROUT_MGCV_GAMMA_AD_REML <- gam(TL ~ s(AGE, bs = "ad", k = 20) + s(LAKE, bs = "re"),
                                  family = Gamma(link = "log"),
                                  method = "REML",
                                  data = TROUT)

summary(m_TROUT_MGCV_GAMMA_AD_REML)


# Model Predictive Performance --------------------------------------------


## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the final retained model that relies on REML estimation

model <- m_TROUT_MGCV_GAMMA_AD_REML

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
Total_df <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - Total_df) * (100 - D2))
cbind(D2, D2_adj)

## Estimate the contribution of the fixed effect AGE and the random intercept LAKE
## to the estimated D2 (above) with "gam.hp"

gam.hp(m_TROUT_MGCV_GAMMA_AD_REML)

## Estimate the pseudo-R2 of Nagelkerke (1991) for the Gamma GAMM built with
## "gamlss" to describe a monotonously increasing relationship

Rsq(m_TROUT_GAMLSS_GAMMA_PBM)


# Model Predictions -------------------------------------------------------


## Visualise the partial residuals of the fixed effect AGE and the random
## intercept LAKE (factor) at both the linear predictor and response scales with
## "visreg" for an "mgcv" object. For "gamlss", only the central tendency would
## be shown.

visreg(m_TROUT_MGCV_GAMMA_AD_REML, "AGE", scale = "linear")

visreg(m_TROUT_MGCV_GAMMA_AD_REML, "AGE", scale = "response")

visreg(m_TROUT_MGCV_GAMMA_AD_REML, "LAKE", scale = "linear")

visreg(m_TROUT_MGCV_GAMMA_AD_REML, "LAKE", scale = "response")

## Predictions from the best-retained Gamma GAMM with REML (Fig. 1c). The argument
## exclude = "s(LAKE)" is use to provide predictions for the fixed effects only.
## Therefore, only one level of the random effect is needed in the new data upon
## which the predictions of the Gamma GAMM are made.

model <- m_TROUT_MGCV_GAMMA_AD_REML

nd_TROUT <- data.frame(AGE = seq(1, 50, 0.1),
                       LAKE = "BRULE")

fitted <- predict(model,
                  nd_TROUT,
                  type = "link",
                  exclude = "s(LAKE)",
                  se.fit = TRUE)

estimate <- exp(fitted$fit)
lower_ci <- exp(fitted$fit - 2 * fitted$se.fit)
upper_ci <- exp(fitted$fit + 2 * fitted$se.fit)

cbind(nd_TROUT$AGE, estimate, lower_ci, upper_ci)

## Obtain predictions from the best-retained Gamma GAMM in "gamlss" (Fig. 1c).
## Note that using new data does not allow the se.fit = TRUE argument, so no
## uncertainty band can be estimated.

nd_TROUT_ALL <- data.frame(
  AGE = seq(1, 50, 0.1))

nd_TROUT_BRULE <- data.frame(
  AGE = seq(1, 50, 0.1),
  LAKE = "BRULE")

nd_TROUT_REDDOG <- data.frame(
  AGE = seq(1, 50, 0.1),
  LAKE = "REDDOG")

nd_TROUT_VOLTZ <- data.frame(
  AGE = seq(1, 50, 0.1),
  LAKE = "VOLTZ")

model <- m_TROUT_GAMLSS_GAMMA_PBM

fitted_BRULE <- predict(model,
                        what = "mu", # default
                        newdata = nd_TROUT_BRULE,
                        type = "response")

fitted_REDDOG <- predict(model,
                         what = "mu",
                         newdata = nd_TROUT_REDDOG,
                         type = "response")

fitted_VOLTZ <- predict(model,
                        what = "mu",
                        newdata = nd_TROUT_VOLTZ,
                        type = "response")

fitted_ALL <- cbind(fitted_BRULE, fitted_REDDOG, fitted_VOLTZ)

## Calculate the mean of the fitted values across the sampling locations (i.e., the 3 LAKEs) 
## to get the predicted central tendency, referred to as pred_average below.

pred_average <- rowMeans(fitted_ALL)

cbind(nd_TROUT_ALL, pred_average)


# References --------------------------------------------------------------


## Anderson TW, Darling DA (1954) A test of goodness of fit. J Am Stat Assoc 49:
##    765–769

## Chhikara RS, Folks JL (1988) The inverse Gaussian distribution. Theory: 
##    Methodology, and Applications. Marcel Dekker, New York.

## Forbes C, Evans M, Hastings N et al. (2010) Gamma distribution. In: Forbes C,
##    Evans M, Hastings N et al (eds) Statistical distributions, 4th edn. John 
##    Wiley & Sons, New Jersey, pp 109–113.

## Gross J, Ligges U (2015) nortest: Tests for Normality. R package version 1.0-4.
##    https://CRAN.R-project.org/package=nortest

## Jarque CM, Bera AK (1980) Efficient tests for normality, homoscedasticity and
##    serial independence of regression residuals. Econ Lett 6:255–259

## Jørgensen B (1987) Exponential dispersion models. J R Stat Soc B 49:127–162

## Komsta L, Novomestky F (2022) moments: Moments, Cumulants, Skewness, Kurtosis
##    and Related Tests. R package version 0.14.1.
##    https://CRAN.R-project.org/package=moments

## Miguez F (2023) nlraa: Nonlinear Regression for Agricultural Applications.
##    R package version 1.9.7. https://CRAN.R-project.org/package=nlraa

## Shapiro SS, Wilk MB (1965) An analysis of variance test for normality (complete
##    samples). Biometrika 52:591-611
