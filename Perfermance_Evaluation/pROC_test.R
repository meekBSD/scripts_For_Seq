library(pROC)
data(aSAH)

roc(outcome ~ s100b, aSAH)
a <-roc(outcome ~ s100b, aSAH)

a$call
  ## roc.formula(formula = outcome ~ s100b, data = aSAH)
ci(a)
  ## 95% CI: 0.6301-0.8326 (DeLong)
ci.se(a)

spec.ci <- ci.sp(a, specificities=seq(0, 100, 5))

dat<- data.frame(pred=c(rep("P",20), rep("N",7)), ref=rep(5,27))
b <- roc(pred ~ ref, dat)
  ## Setting levels: control = N, case = P
  ## Setting direction: controls < cases

sens.ci <- ci.se(b, sensitivities=seq(0, 100, 1))
sens.ci

dat<- data.frame(pred=c(rep("P",20), rep("N",7)), ref=c(rep(5,6), rep(5.9,10), rep(2.3,11)))
b <- roc(pred ~ ref, dat)
  ## Setting levels: control = N, case = P
  ## Setting direction: controls < cases
sens.ci <- ci.se(b, sensitivities=seq(0, 100, 1))
sens.ci

## other reference : https://www.rdocumentation.org/packages/GenBinomApps/versions/1.0-2/topics/clopper.pearson.ci
