library(causl)
library(data.table)

copula_obs <- function(n, rho = 0, seed = 111){
  fam <- list(c(1,3,4),c(5,5,5,1,2),c(5),c(1,3,4))  
  forms <- list(c(X1 ~ 1, X2 ~ X1, X3~ X1), 
                list(A ~ X1 + X3 + C2, C1 ~ 1, C2 ~ C1, C3 ~ C1:C2, C4 ~ C1), 
                Y ~ A + C1 + I(sin(C4)) + A:C1 + A:C2 + A:I(C3>0) + A:I(C4^2), 
                ~ C2) 
  pars <- list(X1 = list(beta=0, phi=1),
               X2 = list(beta=c(0.1,0.2), phi=1), 
               X3 = list(beta=c(0.1,0.1), phi=1), 
               A = list(beta=c(-1,2,2,-1)),
               C1 = list(beta=c(0)),
               C2 = list(beta=c(-2,1)),
               C3 = list(beta=c(0,0.1), phi = 1),
               C4 = list(beta=c(0,0.1), phi=0.1, par2=20),
               Y = list(beta=c(-2, -0.2, 0.3, 0.4, -rho, -rho, rho, rho)), 
               cop = list(Y=list(X1=list(beta=c(0,0.5)),
                                 X2=list(beta=c(-0.5,0)),
                                 X3=list(beta=c(0.5,0))))) 
  link <- list(c("identity","log","logit"),
               c("logit","logit","logit","identity","identity"),
               "log")
  
  set.seed(seed)
  dat <- as.data.table(rfrugalParam(n, formulas=forms, pars=pars, 
                                    family=fam, link=link))
  
  # construct a model matrix and calculate the true CATE
  dat2 <- copy(dat)
  mm <- model.matrix(forms[[3]],dat2[,A:=1])-model.matrix(forms[[3]],dat2[,A:=0]) 
  dat$CRTE <- exp(mm %*% pars$Y$beta)
  #summary(mm %*% pars$Y$beta)
  return(dat)
}

test1 <- copula_obs(n = 10000, rho = 0.25)
mean(test1$A)
x = test1[, !c("A", "Y", "CRTE"), with=FALSE]
forest.hat = causal_forest(x, test1$Y, test1$A, num.trees=2)
hist(forest.hat$W.hat)
summary(forest.hat$W.hat)

test2 <- copula_obs(10000, 0)
crtes = c(0, summary(test2$CRTE))
for(rho in seq(0.25, 1, 0.25)){
  test2 <- copula_obs(10000, rho)
  crtes = rbind(crtes, c(rho, summary(test2$CRTE)))
}
crtes

# Plot the distribution of true CRTE by rho level
df <- data.frame(
  value = c(copula_obs(10000, 0.25)$CRTE, 
            copula_obs(10000, 0.5)$CRTE, 
            copula_obs(10000, 0.75)$CRTE),
  group = factor(rep(c("0.25", "0.5", "0.75"), each = 10000))
)
vline_value <- mean(copula_obs(10000, 0)$CRTE)

ggplot(df, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) + xlim(0, 2) +   theme_minimal() +
  geom_vline(xintercept = vline_value, color = "black", 
             linetype = "dashed", size = 1) + 
  scale_fill_manual(values = c("blue", "red", "green")) + 
  labs(title = "Overlapping Density Curves",
       x = "Value", y = "Density", fill = "Group")
