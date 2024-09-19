library(causl)
library(data.table)
library(grf)
library(devtools)
#install_github("vshirvaikar/rrcf", subdir = "r-package/rrcf")
library(rrcf)
library(ggplot2)
library(dplyr)

trials = 100
ns = c(2500, 5000, 7500, 10000)
rhos = c(0.25, 0.5, 0.75)
pvals = c()
mapes = c()

for(trial in 1:trials){
  trial.start = Sys.time()
  for(n in ns){
    for(rho in rhos){
      run.start = Sys.time()
      data = copula_obs(n, rho, trial)
      split <- sample(seq_len(n), size = n*0.8)
      x = data[, !c("A", "Y", "CRTE"), with=FALSE]
      x.train <- x[split, ]
      x.test <- x[-split, ]
      y.train = data$Y[split]
      y.test = data$Y[-split]
      t.train = data$A[split]
      t.test = data$A[-split]
      crte.test = data$CRTE[-split]
      
      forest.grf <- causal_forest(x.train, y.train, t.train,
                                  num.trees=500, seed=1234)
      forest.glm = rr_causal_forest(x.train, y.train, t.train,
                                    rct=FALSE, num.trees=500, seed=1234)
      pred.grf = rr_predict(forest.grf, x.test)
      pred.glm = rr_predict(forest.glm, x.test)
      mape.grf = mean(abs((crte.test-pred.grf)/crte.test))
      mape.glm = mean(abs((crte.test-pred.glm)/crte.test))
      
      anova.data = data.frame(cbind(y.test, t.test, x.test))
      model.base = glm(y.test ~ ., family = poisson, data = anova.data)
      model.grf = glm(y.test ~ ., family = poisson, 
                      data = cbind(anova.data, t.test*log(pred.grf)))
      model.glm = glm(y.test ~ ., family = poisson, 
                      data = cbind(anova.data, t.test*log(pred.glm)))
      anova.grf = anova(model.base, model.grf)
      anova.glm = anova(model.base, model.glm)
      p.grf = 1 - pchisq(anova.grf$Deviance[2], df = 1)
      p.glm = 1 - pchisq(anova.glm$Deviance[2], df = 1)
      #p.grf = rr_test_calibration(forest.grf, x.test, y.test, t.test)
      #p.glm = rr_test_calibration(forest.glm, x.test, y.test, t.test)
      
      pvals = rbind(pvals, c(n, rho, p.grf, p.glm))
      mapes = rbind(mapes, c(n, rho, mape.grf, mape.glm))
      
      vi.grf = variable_importance(forest.grf)
      vi.glm = variable_importance(forest.glm)
      vi.output = data.frame(vi.grf, vi.glm)
      row.names(vi.output) = colnames(x.train)
      
      print(paste0("Trial ", trial, ", n = ", n, ", rho = ", rho))
      print(c(p.grf, p.glm))
      print(c(mape.grf, mape.glm))
      print(Sys.time()-run.start)
      path = paste0("T", sprintf("%03d", trial), "n", sprintf("%05d", n), 
                    "rho", sprintf("%02d", rho*100))
      save.image(paste0("~/Desktop/GRF/TrialEnvsObs/", path, ".Rdata"))
    }
  }
  print(Sys.time()-trial.start)
}
