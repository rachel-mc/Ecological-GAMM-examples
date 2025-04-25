## Required packages:

library(FSA)
library(gam.hp)
library(gamlss)
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

file_name <- "https://raw.githubusercontent.com/rachel-mc/Ecological-GAMM-examples/main/WALLEYE.txt"
WALLEYE <- read.delim(file_name)


## Poisson GAMM with the default smooth function, s(), for all predictors, thus
## thin plate regression splines (bs = "tp") and k = 10, with ML estimation and 
## AREA as a random intercept in "mgcv"

WALLEYE$AREA <- factor(WALLEYE$AREA)

m_WALLEYE_MGCV_POISSON_FULL <- gam(N ~ s(YEAR) + s(CONDUCT) + s(D1AUG) + s(DEPTH)
                                   + s(EFFORT) + s(TEMP) + s(TURBID) + s(AREA, bs = "re"),
                                   family = poisson,
                                   method = "ML",
                                   data = WALLEYE)

summary(m_WALLEYE_MGCV_POISSON_FULL)


## Equivalent Poisson GAMM in "gamlss" with P-splines, pb(), using the RS
## algorithm which maximises the penalised likelihood instead.
## Note that as opposed to the Poisson GAMM built with "mgcv" above, the one
## built with "gamlss" did not converge with the n = 20 iterations (n.cyc)
## used as default. Below, up to n.cyc = 200 are allowed in an attempt to reach
## convergence. This can take several seconds to run.

m_WALLEYE_GAMLSS_POISSON_FULL <- gamlss(N ~ pb(YEAR) + pb(CONDUCT) + pb(D1AUG) +
                                        pb(DEPTH) + pb(EFFORT) + pb(TEMP) + pb(TURBID)
                                        + random(AREA),
                                        family = PO,
                                        data = WALLEYE,
                                        control = gamlss.control(n.cyc = 200))

summary(m_WALLEYE_GAMLSS_POISSON_FULL)


## Check for likely overdispersion issues from both Poisson GAMMs

check_overdispersion(m_WALLEYE_MGCV_POISSON_FULL)

check_overdispersion(m_WALLEYE_GAMLSS_POISSON_FULL)


## Check for zero-inflation issues from both Poisson GAMMs

check_zeroinflation(m_WALLEYE_MGCV_POISSON_FULL)

check_zeroinflation(m_WALLEYE_GAMLSS_POISSON_FULL)


## Attempt to handle both overdispersion and zero-inflation with the type-II 
## negative binomial (NB2) extension instead with both packages. Note that the
## NB2 is referred to as NBI with "gamlss", whereas the NB1 is referred to as NBII.

WALLEYE$AREA <- factor(WALLEYE$AREA)

m_WALLEYE_MGCV_NB2_FULL <- gam(N ~ s(YEAR) + s(CONDUCT) + s(D1AUG) + s(DEPTH) +
                               s(EFFORT) + s(TEMP) + s(TURBID) + s(AREA, bs = "re"),
                               family = nb,
                               method = "ML",
                               data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_FULL)


m_WALLEYE_GAMLSS_NB2_FULL <- gamlss(N ~ pb(YEAR) + pb(CONDUCT) + pb(D1AUG) + pb(DEPTH)
                                    + pb(EFFORT) + pb(TEMP) + pb(TURBID) + random(AREA),
                                    family = NBI,
                                    data = WALLEYE)

summary(m_WALLEYE_GAMLSS_NB2_FULL)


## To first determine whether each of these NB2 GAMMs either fitted with "mgcv"
## or "gamlss" could be reduced by removing some predictors that may explain 
## little to no variation in the response variable, two different approaches
## that are specific to each package will be used.


## Variable selection by more heavily penalising each smooth term in "mgcv" with
## the argument select = TRUE

m_WALLEYE_MGCV_NB2_SELECT <- gam(N ~ s(YEAR) + s(CONDUCT) + s(D1AUG) + s(DEPTH) +
                                 s(EFFORT) + s(TEMP) + s(TURBID) + s(AREA, bs = "re"),
                                 family = nb,
                                 method = "ML",
                                 select = TRUE,
                                 data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_SELECT)


## The predictors s(CONDUCT), s(EFFORT), and s(TEMP) have their respective
## effective degrees of freedom (edf) shrunk toward zero (0), indicating that
## these three covariates fitted as smooth terms can be discarded. The NB2 GAMM
## is refitted without these predictors as m_WALLEYE_MGCV_REDUCED_NB2

m_WALLEYE_MGCV_NB2_REDUCED <- gam(N ~ s(YEAR) + s(D1AUG) + s(DEPTH) + s(TURBID) +
                                  s(AREA, bs = "re"),
                                  family = nb,
                                  method = "ML",
                                  data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_REDUCED)


## The smooth terms s(D1AUG) and s(DEPTH) now both exhibit edf = 1, indicating
## that both predictor can be included as parametric components instead. A final
## NB2 GAMM is thus fitted with fewer smooth terms

m_WALLEYE_MGCV_NB2_FINAL <- gam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID) +
                                s(AREA, bs = "re"),
                                family = nb,
                                method = "ML",
                                data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_FINAL)


## For the FULL NB2 GAMM that was fitted in "gamlss", the stepGAIC() function 
## implemented in this package will be used instead to perform a stepwise 
## selection. This can take a few minutes to run.

stepGAIC(m_WALLEYE_GAMLSS_NB2_FULL)


## The reduced "gamlss" NB2 GAMM has discarded the same smooth termps as "mgcv",
## in addition of pb(D1AUG). To allows comparisons on a more similar basis between
## the two packages, pb(D1AUG) is kept but converted to a parametric component
## in the final GAMM built with "gamlss", with the same done to pb(DEPHT) refitted
## as DEPTH, as for the final NB2 GAMM built in "mgcv".

m_WALLEYE_GAMLSS_NB2_FINAL <- gamlss(N ~ pb(YEAR) + D1AUG + DEPTH + pb(TURBID)
                                     + random(AREA),
                                     family = NBI,
                                     data = WALLEYE)

summary(m_WALLEYE_GAMLSS_NB2_FINAL)


## check if overdispersion and zero-inflation were sufficiently accounted for in
## the retained final NB2 GAMM, which can only be done for the one built in 
## "mgcv" when the NB2 extension is used.

check_overdispersion(m_WALLEYE_MGCV_NB2_FINAL)

check_zeroinflation(m_WALLEYE_MGCV_NB2_FINAL)


## Because "gamlss" can accommodate other extensions of the Poisson distribution,
## including the type-I negative binomial (NB1; Hilbe 2014), the double Poisson
## (DPO; Efron 1986), a parametrisation of the generalised Poisson (GPO; Consul 1989),
## and the Poisson inverse-Gaussian (PIG; Holla 1967) among others, this four
## extensions are also fitted.

m_WALLEYE_GAMLSS_NB1_FINAL <- gamlss(N ~ pb(YEAR) + D1AUG + DEPTH + pb(TURBID)
                                     + random(AREA),
                                     family = NBII, 
                                     data = WALLEYE)

summary(m_WALLEYE_GAMLSS_NB1_FINAL)


m_WALLEYE_GAMLSS_DPO_FINAL <- gamlss(N ~ pb(YEAR) + D1AUG + DEPTH + pb(TURBID)
                                     + random(AREA),
                                     family = DPO,
                                     data = WALLEYE)

summary(m_WALLEYE_GAMLSS_DPO_FINAL)


m_WALLEYE_GAMLSS_GPO_FINAL <- gamlss(N ~ pb(YEAR) + D1AUG + DEPTH + pb(TURBID)
                                     + random(AREA),
                                     family = GPO,
                                     data = WALLEYE)

summary(m_WALLEYE_GAMLSS_GPO_FINAL)  


m_WALLEYE_GAMLSS_PIG_FINAL <- gamlss(N ~ pb(YEAR) + D1AUG + DEPTH + pb(TURBID)
                                     + random(AREA),
                                     family = PIG,
                                     data = WALLEYE)

summary(m_WALLEYE_GAMLSS_PIG_FINAL)


## The NB1 and DPO failed to converge with the default 20 n.cyc, as for the Poisson, 
## and were discarded. The GPO and PIG converged, but were found to be inadequate 
## from using "hnp" when following the same adequacy assessment approach as for 
## the NB2 below. As such, the NB2 was retained as the best distribution family.


############################# MODEL ADEQUACY ###################################

## Diagnostic plots with the gam.check() function of "mgcv"

gam.check(m_WALLEYE_MGCV_NB2_FINAL)


## The same diagnostic plots with the appraise() function of "gratia" with a 
## simulated, theoretically-based envelope relying on deviance residuals (default)

appraise(m_WALLEYE_MGCV_NB2_FINAL, method = "simulate")


## Worm plot with the worm_plot() function of "gratia" with a simulated envelope

worm_plot(m_WALLEYE_MGCV_NB2_FINAL, method = "simulate", n_simulate = 100)


## Rootogramm with the rg() function of "gratia"

rg <- rootogram(m_WALLEYE_MGCV_NB2_FINAL)
draw(rg)


## Adequacy assessment based on half-normal plot with "hnp"
## This can take many seconds to a few minutes to run

model <- m_WALLEYE_MGCV_NB2_FINAL
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

## 1 simulation only to produce 1 half-normal plot
## This can take many seconds to run given the GAMM's complexity
hnp(model,
    newclass = TRUE,
    diagfun = dfun,
    simfun = sfun,
    fitfun = ffun,
    how.many.out = TRUE,
    plot.sim = TRUE,
    paint = TRUE,
    ylab = "deviance residuals")


## 10 simulations
## This can take many minutes to run given the GAMM's complexity
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

round(return_max(hnp_summary),2)

Summarize(hnp_summary)


## Calculate the mgcViz score (mode) and mean for a given model based on 
## 100 iterations. A 95% "uncertainty interval" (Ui) is also calculated.
## This can take several seconds to run

model <- m_WALLEYE_MGCV_NB2_FINAL
predictor <- "YEAR"

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

plot(m_WALLEYE_GAMLSS_NB2_FINAL)


## Worm plot with the wp() function of "gamlss"

wp(m_WALLEYE_GAMLSS_NB2_FINAL)


## Bucket plot with the bp() function of "gamlss"

bp(m_WALLEYE_GAMLSS_NB2_FINAL)


## Detrended Transformed Owen's Plot (DTOP) with the dtop() function of "gamlss"

dtop(m_WALLEYE_GAMLSS_NB2_FINAL)


## Helper functions contained in the hnp_gamlss_count() function to allow the
## adequacy assessment of GAMM fitted in "gamlss" with different exponential
## families such as the NB2, NB1, DPO, GPO, PIG, and others.

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


## 1 "hnp" run
## This can take many minutes as assessing a gamlss object is computationnally-
## intensive. If too long to perform, use argument sim = 19 and conf = 1 (see
## Moral et al. 2017)

model <- m_WALLEYE_GAMLSS_NB2_FINAL

hnp_gamlss_count(model,
                 how.many.out = TRUE,
                 plot.sim = TRUE, 
                 paint = TRUE)


## 10 "hnp" runs
## This can take many minutes as assessing a gamlss object is especially 
## computationally intensive. If too long to perform, use argument sim = 19 and
## conf = 1 (see Moral et al. 2017)

set.seed(2025)
hnp_obj <- list()
  for(i in 1:10) {
  hnp_obj[[i]] <- hnp_gamlss_count(model, 
                                how.many.out = TRUE, 
                                plot.sim = FALSE) 
}

hnp_summary <- sapply(hnp_obj, function(x) x$out/x$total*100) 

return_max <- function(numvec){
  dens <- density(numvec)
  return(dens$x[which.max(dens$y)][1])
}

round(return_max(hnp_summary),2)

Summarize(hnp_summary)


############################# MODEL SELECTION ##################################

## Using the compareML() function of "itsadug" show that the final NB2 GAMM offers
## an equivalent fit to that of the intial full counterpart despite being much
## less complex, as captured by the Estimated degrees of freedome (Edf).

compareML(m_WALLEYE_MGCV_NB2_FINAL,
          m_WALLEYE_MGCV_NB2_FULL)


## Comparing the same two models as above under an information-theoretic approach
## indicates that the final NB2 GAMM is much more parsimonious than its full
## counterpart using the model.sel() function of "MuMIn" according to the AICc

model.sel(m_WALLEYE_MGCV_NB2_FINAL,
          m_WALLEYE_MGCV_NB2_FULL)


## The test based on the Bayes Factor (BF) in performance also favours the final
## rather than the initial full NB2 GAMM built in "mgcv"

test_bf(m_WALLEYE_MGCV_NB2_FINAL,
        m_WALLEYE_MGCV_NB2_FULL)


## For the "gamlss" final NB2 GAMM, given that it is nested in the full NB2 GAMM,
## the implemented LR.test() function is used to compare the two NB2 GAMMs.

LR.test(m_WALLEYE_GAMLSS_NB2_FINAL,
        m_WALLEYE_GAMLSS_NB2_FULL)


## Comparing the two NB2 GAMMs under an information-theoretic approach

model.sel(m_WALLEYE_GAMLSS_NB2_FINAL,
          m_WALLEYE_GAMLSS_NB2_FULL)


## The test based on the Bayes Factor (BF) provides the same interpretation

test_bf(m_WALLEYE_GAMLSS_NB2_FINAL,
        m_WALLEYE_GAMLSS_NB2_FULL)


##################### MODEL PREDICTIVE PERFORMANCE #############################

## Estimate the deviance explained (D2) and its adjusted version (D2_adj) for
## the final retained model relying on REML estimation

model <- m_WALLEYE_MGCV_NB2_FINAL 

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
Total_df <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - Total_df) * (100 - D2))
cbind(D2, D2_adj)


## Estimate the contribution of the fixed AGE effect and random LAKE factor
## to the estimated D2 (above) with "gam.hp"

gam.hp(m_WALLEYE_MGCV_NB2_FINAL)


## Estimate the pseudo-R2 of Nagelkerke (1991) for the Gamma GAMM built with
## "gamlss" to describe a monotonously increasing relationship

Rsq(m_WALLEYE_GAMLSS_NB2_FINAL)


############################# MODEL PREDICTIONS ################################

## predictions with "mgcv" after refitting the final NB2 GAMM with REML (Fig. 2a)

m_WALLEYE_MGCV_NB2_FINAL_REML <- gam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID)
                                     + s(AREA, bs = "re"),
                                     family = nb,
                                     method = "REML",
                                     data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_FINAL_REML)


## predictions

nd_WALLEYE <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                         D1AUG = mean(WALLEYE$D1AUG),
                         DEPTH = mean(WALLEYE$DEPTH),
                         TURBID = mean(WALLEYE$TURBID),
                         AREA = "ALSP")


## check the new data that were just created
headtail(nd_WALLEYE)

model <- m_WALLEYE_MGCV_NB2_FINAL_REML
YEAR <- nd_WALLEYE$YEAR

fitted <- predict(model, nd_WALLEYE, type = "link", exclude = "s(AREA)", se.fit = TRUE)
estimate <- exp(fitted$fit)
lower_ci <- exp(fitted$fit - 1.96 * fitted$se.fit)
upper_ci <- exp(fitted$fit + 1.96 * fitted$se.fit)

cbind(YEAR, estimate, lower_ci, upper_ci)


## predictions with "gamlss" (Fig. 2a)

nd_WALLEYE_all <- data.frame(YEAR = seq(2001, 2021, by = 0.1))

nd_WALLEYE_ALSP <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "ALSP")

nd_WALLEYE_BEBA <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "BEBA")

nd_WALLEYE_LDDM <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "LDDM")

nd_WALLEYE_LSFR <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "LSFR")

nd_WALLEYE_LSLO <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "LSLO")

nd_WALLEYE_LSPN <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "LSPN")

nd_WALLEYE_LSPS <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "LSPS")

nd_WALLEYE_MTSL <- data.frame(YEAR = seq(2001, 2021, by = 0.1),
                              D1AUG = mean(WALLEYE$D1AUG),
                              DEPTH = mean(WALLEYE$DEPTH),
                              TURBID = mean(WALLEYE$TURBID),
                              AREA = "MTSL")

model <- m_WALLEYE_GAMLSS_NB2_FINAL

fitted_ALSP <- predict(model,
                       what = "mu",
                       newdata = nd_WALLEYE_ALSP,
                       type = "response")

fitted_BEBA <- predict(model,
                       what = "mu",
                       newdata = nd_WALLEYE_BEBA,
                       type = "response")

fitted_LDDM <- predict(model,
                       what = "mu",
                       newdata = nd_WALLEYE_LDDM,
                       type = "response")

fitted_LSFR <- predict(model,
                       what = "mu",
                       newdata = nd_WALLEYE_LSFR,
                       type = "response")

fitted_LSLO <- predict(model, 
                       what = "mu", 
                       newdata = nd_WALLEYE_LSLO,
                       type = "response")

fitted_LSPN <- predict(model,
                       what = "mu", 
                       newdata = nd_WALLEYE_LSPN, 
                       type = "response")

fitted_LSPS <- predict(model,
                       what = "mu", 
                       newdata = nd_WALLEYE_LSPS,
                       type = "response")

fitted_MTSL <- predict(model, 
                       what = "mu", 
                       newdata = nd_WALLEYE_MTSL,
                       type = "response")

fitted_ALL <- cbind(fitted_ALSP,
                    fitted_BEBA,
                    fitted_LDDM, 
                    fitted_LSFR,
                    fitted_LSLO,
                    fitted_LSPN,
                    fitted_LSPS,
                    fitted_MTSL)

## Calculate the mean from these 8 sampling locations (i.e., AREA) to get the
## predicted central tendency referred to as fitted_ALL below

fitted_ALL <- transform(fitted_ALL, 
                        average = ((fitted_ALSP + fitted_BEBA + fitted_LDDM +
                                    fitted_LSFR + fitted_LSLO + fitted_LSPN + 
                                    fitted_LSPS + fitted_MTSL) / 8))

pred_average <- fitted_ALL$average 

cbind(nd_WALLEYE_all,
      pred_average)


########################## TEMPORAL AUTOCORRELATION ############################

## Check for temporal autocorrelation for the final NB2 GAMM with REML ("mgcv")
## and for the final NB2 GAMM built in "gamlss" according to test of Durbin and
## Watson (1950) implemented in "performance" which detect such analytical issues.

check_autocorrelation(m_WALLEYE_MGCV_NB2_FINAL_REML)

check_autocorrelation(m_WALLEYE_GAMLSS_NB2_FINAL)


## Estimate the first-oder autoregressive (AR1) coefficient, referred to as rho
## with "itsadug" from refitting the final NB2 GAMM with fast REML (fREML) instead
## or REML with the bam() instead of the gam() function in "mgcv". The argument
## discrete = TRUE is also used to accelerate the process (see van Rij et al. 2022
## and Wood 2017).

m_WALLEYE_MGCV_NB2_FINAL_fREML <- bam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID)
                                      + s(AREA, bs = "re"),
                                      family = nb, 
                                      method = "fREML", 
                                      discrete = TRUE, 
                                      data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_FINAL_fREML)

RHO <- start_value_rho(m_WALLEYE_MGCV_NB2_FINAL_fREML,
                       plot = TRUE)
RHO


## As rho (RHO) = 0.111, this indicates that approximately 11% of the variation
## in the response variable (N) in a given year is explained by that observed in
## the previous year. The same NB2 GAMM is refitted to now account for rho and
## AR1 is added to identify it.

m_WALLEYE_MGCV_NB2_FINAL_fREML_AR1 <- bam(N ~ s(YEAR) + D1AUG + DEPTH + s(TURBID)
                                          + s(AREA, bs = "re"),
                                          family = nb,
                                          method = "fREML",
                                          discrete = TRUE,
                                          rho = RHO,
                                          data = WALLEYE)

summary(m_WALLEYE_MGCV_NB2_FINAL_fREML_AR1)

## Visually inspect whether the addition of an autocorrelation matrix with RHO
## as allowed to reduce the temporal autocorrelation by looking again at the 
## ACF plot with "itasdug"

check_resid(m_WALLEYE_MGCV_NB2_FINAL_fREML_AR1,
            select = 3)


## Use the compareML() function of itsadug to compare both models (with and without
## RHO being accounted for)

compareML(m_WALLEYE_MGCV_NB2_FINAL_fREML,
          m_WALLEYE_MGCV_NB2_FINAL_fREML_AR1)


## Compare both NB2 GAMMs under an information-theoretic approach using AICc

model.sel(m_WALLEYE_MGCV_NB2_FINAL_fREML,
          m_WALLEYE_MGCV_NB2_FINAL_fREML_AR1)


## Compare both NB2 GAMMs with the test relying on the Bayes Factor (BF)

test_bf(m_WALLEYE_MGCV_NB2_FINAL_fREML,
        m_WALLEYE_MGCV_NB2_FINAL_fREML_AR1)


## Accounting for the AR1 term (rho) has statistically improved model fit, but
## when its predictions are plotted against those of the NB2 GAMM with REML that
## does not account for temporal autocorrelation, the differences between the
## two regression curves are quite slight.


################################ REFERENCES ####################################

## Consul PC (1989) Generalized Poisson distributions: properties and 
##    applications. Marcel Dekker, New York

## Efron B (1986) Double exponential families and their use in the generalized
##    linear regression. J Am Stat Assoc 81:709–721.

## Holla MS (1967) On a poisson-inverse gaussian distribution. Metrika 11:115–121.
##    https://doi.org/10.1007/BF02613581