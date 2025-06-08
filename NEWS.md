# micsr 0.1-2

* dependences on part of the `tidyverse` and on `knitr` are removed

* `maximize` now has a "newton" method

* `lnl_ivldv` now have an opposite argument

* the `x` and `y` arguments of `scoretest` are replaced by `object`
  and `...` 

* ancil is replaced by instr for the coefficients of the ivldv function

* commit

* the vcov method now has error msg when the estimator is not available

* value is now a matrix with three columns (model, saturated, null)

* logLik has a sum argument (if FALSE, the individual contributions
  are returned)

* last commit

* the analytical hessian for ordreg with censored observations is fixed

* a testing infrastructure is added, using the `tinytest` package

* for models estimated by maximum likelihood, the gradient and the
  hessian are checked if `check_gradient = TRUE` and the result is a
  `check_gradient` elemnt

* weights are added to `weibreg`, `ordreg`, `poisreg`

* in `ordreg`, `binomreg` and `poisreg`, the formula is replaced by terms in the
  result

* `dmills` and `d2mills` replaced by `mills(, deriv = 1/2)`

* `grep` and `invert` argument added to `select_coef`

* new data set `random_group`

* two new vignettes `lm_function` and `micsr_objects`

* `select_coef` is now exported and has a `fixed` argument (if `TRUE`,
  constant parameters are selected), `coef` and `vcov` also have this
  argument

* gaze.lm.RStests replaces gaze.LMtestlist due to a naming change in
  the spdep package

* lhs `y1 | y2` replaced by `y1 + y2` in `pscore` and `escount`

* `"twosteps"` replaced by `"twostep"` in `tobit1`

* `var_ replace` by `.var` in `pscore`

* new `quad_form` function

* `gaze` method for `mlogit` objects

* for `rsq`, `type` = w, lm, lr changed to wald, score and lr, same
  for names in `tests` for `binomreg`, `ordreg` and `poisreg`

* argument `.vcov` of `stder` is changed to `vcov`

* `gaze.ivreg` has now `signif.stars = FALSE`
 
* a `cmtest` method is added for `weibreg` objects

* new data sets `unemp_duration` and `recall`

* `ordreg` now has a `"cloglog"` link and a `Surv` object can be
  provided as the response for right-censored data
  
* `weibrug` is a new function to fit Weibull models; the AFT and PH
  parametrization are provided and the Gamma-Weibull model is obtained
  with `mixing = TRUE`

# micsr 0.1-1

* initial version of micsr on cran
