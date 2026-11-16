trials = 100
ns = c(2500, 5000, 7500, 10000)
rhos = c(0, 0.25, 0.5, 0.75)
trees = 500

pvals = c()
mapes = c()
vimps = c()
features = NULL

save_dir = file.path("outputs/obs")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
message(paste0("Saving trial outputs to: ", normalizePath(save_dir, mustWork = FALSE)))

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
                                  num.trees=trees, seed=1234)
      forest.glm = rr_causal_forest(x.train, y.train, t.train,
                                    rct=FALSE, num.trees=trees, seed=1234)
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

      vi.grf = variable_importance(forest.grf)
      vi.glm = variable_importance(forest.glm)
      features = colnames(x.train)
      
      pvals = rbind(pvals, c(n, rho, p.grf, p.glm))
      mapes = rbind(mapes, c(n, rho, mape.grf, mape.glm))
      for(i in seq_along(features)){
        vimps = rbind(vimps, c(trial, n, rho, i, vi.grf[i], vi.glm[i]))
      }
      
      print(paste0("Trial ", trial, ", n = ", n, ", rho = ", rho))
      print(c(p.grf, p.glm))
      print(Sys.time()-run.start)
      path = paste0("T", sprintf("%03d", trial), "n", sprintf("%05d", n), 
                    "rho", sprintf("%02d", rho*100))
      save.image(file.path(save_dir, paste0(path, ".Rdata")))
    }
  }
}
