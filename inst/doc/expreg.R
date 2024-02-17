## ----message = FALSE----------------------------------------------------------
library("micsr")
library("dplyr")
cigmales <- cigmales %>%
    mutate(age2 = age ^ 2, educ2 = educ ^ 2,
           age3 = age ^ 3, educ3 = educ ^ 3,
           educage = educ * age)                                 
pois_cig <- glm(cigarettes ~ habit + price + restaurant + income + age + age2 +
                    educ + educ2 + famsize + race, data = cigmales,
                family = quasipoisson)

## -----------------------------------------------------------------------------
iv_cig <- expreg(cigarettes ~ habit + price + restaurant + income + age + age2 +
                       educ + educ2 + famsize + race | . - habit + age3 + educ3 +
                       educage + lagprice + reslgth, data = cigmales,
                   twosteps = FALSE)
gmm_cig <- update(iv_cig, twosteps = TRUE)

## ----echo = FALSE-------------------------------------------------------------
if (requireNamespace("modelsummary")){
    modelsummary::msummary(list(ML = pois_cig, IV = iv_cig, GMM = gmm_cig),
                           fmt = 4, statistic = 'statistic', coef_omit = "(Intercept)")
}

## -----------------------------------------------------------------------------
sargan(gmm_cig)

## -----------------------------------------------------------------------------
ml_bwt <- glm(birthwt ~ cigarettes + parity + race + sex, data = birthwt,
              family = quasipoisson)
iv_bwt <- expreg(birthwt ~ cigarettes + parity + race + sex |
                      . - cigarettes + edmother + edfather + faminc + cigtax,
                  data = birthwt, method = "iv")
gmm_bwt <- update(iv_bwt, method = "gmm")

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  if (requireNamespace("modelsummary")){
#      modelsummary::msummary(list(ML = ml_bwt, IV = iv_bwt, GMM = gmm_bwt),
#                             fmt = 4, statistic = 'statistic', coef_omit = "(Intercept)")
#  }

## -----------------------------------------------------------------------------
sargan(gmm_bwt)

