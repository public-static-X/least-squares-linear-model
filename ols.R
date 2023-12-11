
rm(list = ls())



ols <- function(formula, data, weights = NULL) {
  tryCatch({
    
    if (missing(data)) {
      X <- model.matrix(formula)
      Y <- model.frame(formula) |> model.response()
      data = model.frame(formula)
    } else {
      X <- model.matrix(formula, data = data)
      Y <- model.frame(formula, data = data) |> model.response()
    }
    
    if (!is.null(weights)) {
      W <- diag(weights)
      Coefficients <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% as.matrix(Y)
    } else {
      Coefficients <- solve(t(X) %*% X) %*% t(X) %*% as.matrix(Y)
    }
    
    Coefs <- as.vector(Coefficients); names(Coefs) <- rownames(Coefficients)
    
    # Create a list with a structure similar to lm output
    result <- list(
      designMatrix = X,
      coefficients = Coefs,
      weights = weights,
      response = as.vector(Y),
      residuals = as.vector(Y - X %*% Coefficients),
      fitted.values = as.vector(X %*% Coefficients),
      call = match.call(),
      terms = terms(formula, data = data),
      qr = as.list(qr(X)),
      rank = qr(X)$rank,
      df.residual = length(Y) - qr(X)$rank,
      contrasts = NULL,  # You can modify this part based on your needs
      xlevels = NULL     # You can modify this part based on your needs
    )
    
  }, error = function(e) {
    print(e)
  })
  
  structure(class = "ols", result)
}

# Weighted Least Squares (WLS) function
wls <- function(formula, data, weights) {
  ols(formula, data, weights)
}

# Example usage:
# wls_result <- wls(y ~ x1 + x2, data = your_data_frame, weights = your_weights_vector)


print.ols <- function(model) {
  cat("\nCall:\n", deparse(model$call), "\n")
  cat("\nCoefficients:\n")
  print(model$coefficients |> round(4), quote = FALSE)
}

plot.ols <- function(model,...) {
  
  plot(model$fitted.values, scale(model$residuals),
       col = (abs(model$residual)+2)^2 , pch=19,
       ylab = "Standardised residuals",
       xlab = "Predicted values from linear model")
  lin <- ols(model$residuals ~ model$fitted)
  df <- data.frame(X = model$fitted.values,
                   Y = lin$fitted.values)
  df <- df[order(df$X),]
  lines(df$X,df$Y,col = "red", lwd=4)
  abline(h = 0,lty=2)
  
  plot(density(model$residuals, bw=0.2),...)
  
  rm(df); rm(lin)
  
}

summary.ols <- function(model) {
  
  errorVariance <- sum(model$residuals ^2) / model$df.residual
  cov_scaled <- 
    if (is.null(model$weights)) {
      # OLS covariance matrix
      errorVariance * solve(t(model$designMatrix) %*% model$designMatrix)
    } else {
      # WLS covariance matrix
      weightMatrix <- diag(model$weights)
      weightedMSE <- sum(model$weights * model$residuals^2) / model$df.residual
      errorVariance <- weightedMSE
      errorVariance * solve(t(model$designMatrix) %*% weightMatrix %*% model$designMatrix)
    }
  
  
  VAR <- cov_scaled |> diag()
  SE <- sqrt(VAR)
  
  T_stat <- model$coefficients/ SE
  p_values <- 2*pt(-abs(T_stat),df = model$df.residual)
  
  Mat <- cbind(model$coefficients,SE,T_stat,p_values) 
  
  
  
  summar <- list()
  
  summar$coefficients <- Mat
  summar$call <- model$call
  
  summar$sigma <- sqrt(errorVariance)
  summar$cov_scaled <- cov_scaled
  
  colnames(summar$coefficients) <- c("Estimate","Std. Err","T value","Pr(> |t|)")
  
  class(summar) <- "summary.ols"
  
  return (summar)
   
}

formulate <- function(response, predictors) {
  return(paste(response,"~",paste(predictors,collapse = "+")) |> as.formula())
}

print.summary.ols <- function(modelSummary) {
  cat("\nCall:\n", deparse(modelSummary$call), "\n\n")
  cat("Coefficients:\n")

  p_values <- modelSummary$coef[,4]
  Significance <- ifelse(p_values <= 0.001, "***",
                   ifelse(p_values <= 0.01, "**", 
                          ifelse(p_values <= 0.05, "*", 
                                 ifelse(p_values < 0.1, ".", "") ) ) ) 
  
  TAB <- cbind(modelSummary$coef |> round(4),Significance)
  TAB[,4] <- ifelse(TAB[,4] == 0, "<.0001",TAB[,4])
  
  TAB[,1:3] <- sapply(TAB[,1:3], function(x) {
                    ifelse(as.numeric(x) > 0,paste0(" ",x),x)
                })
  print(TAB, quote = FALSE)
}

vcov.ols <- function(model) {
  errorVariance <- sum(model$residuals ^2) / model$df.residual
  VAR <- errorVariance * t(model$designMatrix) %*% model$designMatrix |> solve() 
  return (VAR)
}

vcov.summary.ols <- function(Summary) {
  print(Summary$cov_scaled)
}

confint.ols <- function(model, alpha = 0.05) {
  summar <- summary(model)
  mat <- cbind( coef(model) + qt(alpha/2,df = model$df.residual) *  summar$coef[,2],
                coef(model) - qt(alpha/2,df = model$df.residual) *  summar$coef[,2] )
  colnames(mat) <- c("2.5%","97.5%")
  
  mat
}

# Example usage of ols function

mod <- ols(Sepal.Length ~ ., data = iris)
mod
lin <- lm(Sepal.Length ~ ., data = iris)
vcov(mod)

plot(mod)
S <- summary(mod)
S

wls <- ols(Sepal.Length ~ ., data = iris, weights = c(rep(1,100),rep(3.5,50)))
lin_wls <- lm(Sepal.Length ~ ., data = iris, weights = c(rep(1,100),rep(3.5,50)))

