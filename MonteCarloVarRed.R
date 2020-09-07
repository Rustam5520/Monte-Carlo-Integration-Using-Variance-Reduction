Monte_Carlo_Advanced <- function(g, lb, ub, n = 10^5, plot = FALSE, 
                                 var.reduc = FALSE,   # "var.reduc" which takes 
                                 #  a boolean (True or False).
                                 type = "Antithetic", # "type" which takes a string 
                                 # from user as "Antithetic" or "Control".
                                 cv = "NULL") {       # I've desided to add extra parameter 'cv'
  # which takes a list of 'Control Variates' 
  # ex: c(f1, f2...), whenever user                                                         # choses this approach
  #Generatiing X's.-------------------------------------------------------|
  x.classic <- runif(n, lb, ub)
  #Monte Carlo Integration, estimator and variance.-----------------------|
  exp.classic <- mean(sapply(x.classic, g)) * (ub - lb)
  var.classic <- mean(sapply(x.classic, g)^2) - mean(sapply(x.classic, g))^2
  #Finding R integrate function's estiomation.----------------------------|
  exp.R.integrate <- integrate(g, lb, ub)
  #Defining (Error) difference between Real and Estimated values.---------|
  error.classic_vs_R.integrate <- abs(exp.classic - exp.R.integrate$value)
  
  div_line <- strrep("-", 74) # Just a line that will be used to divide output parts. 
  
  if (!var.reduc){
    print(matrix(c(exp.classic, var.classic), nrow = 1,
                 dimnames =  list("Classical Approach", c("Estimated Mean", "Variance"))))
    cat(div_line, "\n|Note: In order to reduce variance you may want to use `Antithetic variables` or `Control Variates` methods|")
  }
  
  # In case of user defining wrong cv (Control Variates) parameter--------|
  if(var.reduc & type == "Control" & !is.list(cv)) {
    cat("Please specify your Control Variates in a list. ex: cv = c(f1, f2...)")
  }
  #-----------------------------------------------------------------------|
  if(var.reduc & type == "Antithetic"){
    x.ant <- runif(n / 2, lb, ub)
    
    exp.ant.pt1 <- mean(sapply(x.ant, g)) * (ub - lb)
    exp.ant.pt2 <- mean(sapply(ub - x.ant, g)) * (ub - lb)
    exp.ant <- (exp.ant.pt1 + exp.ant.pt2) / 2
    
    z1.1 <- g(x.ant)
    z1.2 <- g(ub - x.ant)
    
    var.ant <- (var(z1.1) + var(z1.2) + 2 * cov(z1.1, z1.2)) / 4
    
    result.ant <- matrix(c(exp.classic, exp.ant,
                           var.classic, var.ant),nrow=2,
                         dimnames = list(c("Classical Method","Antithetic variables"),
                                         c("Estimated Mean","Variance")))
    print(result.ant)
    var.red.ant = 100 * ((var.classic - var.ant) / var.classic)
    cat(div_line, "\n|The Antithetic variables approach achived",
        round(var.red.ant, 2), "%", "reduction in variance.|")
  }
  if(var.reduc & type == "Control" & is.list(cv)){
    x.cv <- runif(floor(n / length(cv)), lb, ub)
    if (length(cv) == 1){
      c <- -cov(g(x.cv), cv[[1]](x.cv)) / var(cv[[1]](x.cv))
      
      exp.f <- (1/(ub - lb)) * integrate(cv[[1]], lb, ub)$value
      
      exp.cv.1 <- mean(g(x.cv) + c * (cv[[1]](x.cv) - exp.f)) * (ub - lb)
      var.cv.1 <- var(g(x.cv) + c * (cv[[1]](x.cv) - exp.f))
      
      result.cv.1 <- matrix(c(exp.classic, exp.cv.1,
                              var.classic, var.cv.1),nrow=2,
                            dimnames = list(c("Classical Approach","Control Variates"),
                                            c("Estimated mean","Variance")))
      var.red.cv = 100 * ((var.classic - var.cv.1) / var.classic)
      print(result.cv.1)
      cat(div_line, "\n|The control variates approach achived",
          round(var.red.cv, 2), "%", "reduction in variance.|")  
    }
    else {
      # Matrix of g and f's of X--------------------------------------------|
      g_by_fs <- c(g(x.cv))
      for (i in 1:length(cv)) {
        g_by_fs <- cbind(g_by_fs, cv[[i]](x.cv))
      }
      # --------------------------------------------------------------------|
      # Vecror of Estimators of Control Variates----------------------------|
      mu <- c(1)
      for (i in 1:length(cv)) {
        mu <- c(mu, integrate(cv[[i]], lb, ub)$value)
      }
      # --------------------------------------------------------------------|    
      # Fitting Linear Model------------------------------------------------|
      L <- lm(g_by_fs[, 1] ~ g_by_fs[, -1])
      
      exp.cv <- sum(L$coefficients * mu)
      var.cv <- summary(L)$sigma^2
      
      result.cv <- matrix(c(exp.classic, exp.cv,
                            var.classic, var.cv),nrow=2,
                          dimnames = list(c("Classical Approach","Control variates"),
                                          c("Estimated mean","Variance")))
      var.red.cv = 100 * ((var.classic - var.cv) / var.classic)
      print(result.cv)
      cat(div_line, "\n|The Control Variates approach achived",
          round(var.red.cv, 2), "%", "reduction in variance.|")
    }
  }
  
  # Drawing plot.-----------------------------------------------------------|
  if(plot){
    y <- seq(lb, ub, 0.001)
    y.low <- rep(0, times = length(y))
    plot(y, g(y), type = "l", lwd = 2.5, main = body(g))
    lines(y, y.low, col = 'green')
    lines(y, g(y), col = 'red')
    polygon(c(y, rev(y)), c(g(y), rev(y.low)),col = "grey50", border = NA)
  }
  # ------------------------------------------------------------------------| 
}
