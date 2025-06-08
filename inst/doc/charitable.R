## -----------------------------------------------------------------------------
library(micsr)


## -----------------------------------------------------------------------------
head(charitable, 5)


## -----------------------------------------------------------------------------
charitable$logdon <- with(charitable, log(donation) - log(25))


## -----------------------------------------------------------------------------
#| message: false
char_form <- logdon ~ log(donparents) + log(income) +
    education + religion + married + south
if (requireNamespace("AER")){
    library("AER")
    ml_aer <- tobit(char_form, data = charitable)
}
if (requireNamespace("censReg")){
    library("censReg")
    ml_creg <- censReg(char_form, data = charitable)
}
ml <- tobit1(char_form, data = charitable)


## -----------------------------------------------------------------------------
scls <- update(ml, method = "trimmed")
ols <- update(ml, method = "lm")


## -----------------------------------------------------------------------------
cmtest(ml)


## -----------------------------------------------------------------------------
#| include: false
cmtest(ml_aer)
cmtest(ml_creg)


## -----------------------------------------------------------------------------
cmtest(ml, test = "heterosc")


## -----------------------------------------------------------------------------
cmtest(ml, test = "normality", opg = TRUE)
cmtest(ml, test = "heterosc", opg = TRUE)


## -----------------------------------------------------------------------------
cmtest(ml, test = "skewness")
cmtest(ml, test = "kurtosis")


## -----------------------------------------------------------------------------
#| label: fig-histnorm
#| fig-cap: "Empirical distribution of the response and normal approximation"
#| echo: false
#| message: false
if (requireNamespace("ggplot2")){
    library(ggplot2)
    mu <- mean(subset(charitable, logdon > 0)$logdon)
    std <- sd(subset(charitable, logdon > 0)$logdon)
    ggplot(subset(charitable, logdon > 0), aes(logdon)) +
        geom_histogram(aes(y = after_stat(density)), color = "black",
                       fill = "white", bins = 10) +
        geom_function(fun = dnorm, args = list(mean = mu,
                                               sd = std)) +
        labs(x = "log of charitable giving", y = NULL)
}

