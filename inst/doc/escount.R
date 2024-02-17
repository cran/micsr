## -----------------------------------------------------------------------------
library("micsr")

## ----warning = FALSE----------------------------------------------------------
trips_2s <- escount(trips | car ~ workschl + size + dist + smsa + fulltime + distnod +
                        realinc + weekend + car | . - car - weekend + adults,
                    data = trips)
names(trips_2s)

## -----------------------------------------------------------------------------
summary(trips_2s)

## ----tripsest, warning = FALSE------------------------------------------------
trips_pois <- glm(trips ~ workschl + size + dist + smsa + fulltime + distnod +
                      realinc + weekend + car, data = trips, family = poisson)
trips_ml <- update(trips_2s, method = "ml")

## ----tripsres, warning = FALSE, echo = FALSE----------------------------------
if (requireNamespace("modelsummary")){
    library("modelsummary")
    msummary(list("Poisson" = trips_pois, "2 steps" = trips_2s, "ML" = trips_ml),
             label = "tab:tripsres",
             title = "Estimation results for the trip demand model")
}

## ----drinksest, warning = FALSE, cache = TRUE---------------------------------
kt_pois <- glm(drinks ~ advice + income + age + educ + race + marital + empstatus +
                   region, data = drinks, family = poisson)
kt_ml <- escount(drinks | advice ~ advice + income + age + educ + race + marital +
                     empstatus + region | income + age + educ + race + marital +
                     empstatus + region + medicare + medicaid + champus + hlthins +
                     regmed + dri + limits + diabete + hearthcond + stroke,
                 data = drinks, method = "ml")
kt_2s <- update(kt_ml, method = "twosteps")

## ----drinksres, warning = FALSE, eval = TRUE, echo = FALSE--------------------
msummary(list("Poisson" = kt_pois, "2 steps" = kt_2s, "ML" = kt_ml),
         gof_omit = "AIC|BIC|Log.Lik.",
         label = "tab:drinksres",
         coef_omit = "region|marital|empstatus|age",
         title = "Poisson and endogenous switching models for alcohol demand")

