#' @title A Reference Class to perform Linear Regression
#'
#' @name linreg
#'
#' @description A Class to fit linear models. It can be used to carry out regression.
#'
#' @field formula An object of class formula
#' @field data A data frame from which variables are taken to formula
#'
#' @examples lr_object <- linreg$new(formula = Sepal.Length~Petal.Width, data = iris)
#'
#' lr_object$print()
#' lr_object$plot()
#' lr_object$resid()
#' lr_object$pred()
#' lr_object$coef()
#'
#' @references \href{https://liu.se/insidan/kommunikationsstod/grafiskprofil?l=en}{LiU Theme}
#'
#' @import ggplot2 gridExtra
#'
#' @export linreg
#'
#' @exportClass linreg

linreg = setRefClass(
  "linreg",
  fields = (
    list(
      formula = 'formula',
      data = 'data.frame',
      df_name = 'character',
      use.qr = 'logical',
      reg.coef = 'matrix',
      fitted.vals = 'matrix',
      residuals = 'matrix',
      df = 'numeric',
      res_var = 'matrix',
      var_beta = 'matrix',
      t_beta = 'matrix'
    )
  ),
  methods = list(
    initialize = function(formula,data, ..., use.qr = F){
      "Assign 'formula', 'data' and 'use.qr' upon object creation "
      formula <<- formula
      data <<- data
      use.qr <<- use.qr
      df_name <<- deparse(substitute(data))
      get_beta_coef()
      get_fitted_values()
      get_residuals()
      get_df()
      get_res_var()
      get_var_beta()
      get_t_beta()

    },
    yield_x  = function(){
      "Returns X(Independent Variables) along with the intercept"

      return(model.matrix(formula,data))
    },
    yield_y =function(){
      "Returns Y(dependent variable)"

      return(data[all.vars(formula)[1]])
    },
    get_ixtx = function(){
      'Returns Inverse(T(X)*X) i.e., inverse of dot product of t(x) and x.'
      x = yield_x()
      # Use LU Decomposition
      if (use.qr == F){

        return(solve(crossprod(x,x)))
      }
      # Use QR decomposition
      else{
        return(qr.solve(crossprod(x,x)))
      }
    },
    get_beta_coef = function(){
      'Returns the Beta Coefficient values'

      y = as.matrix(yield_y())
      tx = t(as.matrix(yield_x()))
      ixtx = get_ixtx()
      txy = tx %*% y
      # Calculate Beta
      beta = (ixtx %*% (tx%*%as.matrix(y)))

      reg.coef <<- beta
    },
    get_fitted_values = function(){
      "Return the Predicted value for 'Dependent variable' "
      x = yield_x()
      beta = get_beta_coef()
      yhat =  x %*% beta

      fitted.vals <<- yhat
    },
    get_residuals = function(){
      "Returns the residuals for the dependent variable(Y)"
      residuals <<- as.matrix(yield_y() - fitted.vals)
    },
    get_df = function(){
      'Returns the degree of freedom'
      x = yield_x()
      df <<- nrow(x) - ncol(x)
    },
    get_res_var = function(){
      'Returns Residual of Variance'
      res_var <<- (t(residuals)%*%residuals)/df
    },
    get_var_beta = function(){
      'Returns Variance for each Beta coefficient.'
      var_beta <<- res_var %*% diag(get_ixtx())
      colnames(var_beta) <<- colnames(get_ixtx())
    },
    get_t_beta = function(){
      'Returns t-statistic for each beta coefficient'
      t_beta <<- t(reg.coef)/sqrt(var_beta)
    },

    print = function(){
      'Prints out the Coefficients and Coefficient names'
      cat("linreg(formula = ", format(formula), ", data = ", df_name,")\n\nCoefficients:\n", sep = "")

      print.table(t(reg.coef))
    },
    resid = function(){
      'Returns the vector of Residuals e'
      return(residuals)
    },
    pred = function(){
      'Returns the predicted value y_hat'
      return(fitted.vals)
    },
    coef = function(){
      'Returns the coefficient as a named vector'
      coef_vector = as.vector(reg.coef)
      names(coef_vector) = rownames(reg.coef)
      return(coef_vector)
    },

    plot = function(){
      "Plots the 'Residual vs Fitted' and 'Scale-Location' plots respectively "
      x_label = as.character(paste0("Fitted Values\n",deparse(formula)))

      plot1 <- ggplot(data = data, aes(x=fitted.vals, y=residuals)) +
        geom_point() +
        xlab(x_label) + ylab("Residuals") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(colour = 'turquoise'),
              panel.background = element_rect(fill = 'turquoise'))+
        geom_smooth(method = 'loess', formula = y~x, color = 'red', se=T) +
        ggtitle("Residuals vs Fitted")


      sqrt_res_abs <- sqrt(abs(as.vector(residuals)/as.vector(res_var)))
      plot2 <- ggplot(data = data, aes(x=fitted.vals, y=sqrt_res_abs))+
        geom_point()+
        xlab(x_label) +
        ylab(expression(sqrt("|Standardized residuals|"))) +
        geom_abline(slope = 0, intercept = 0,linetype = "dotted") +
        geom_smooth(method = "loess", formula = y ~ x, color = 'red', se=T) +
        ggtitle("Scale-Location")+
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(colour = 'green'),
              panel.background = element_rect(fill = 'green'))

      grid.arrange(plot1,plot2, nrow=2)

    },
    summary = function(){
      'Function to summarise Inferential statistics summary similar to summary.lm()'
      my.sign = function(x) {
        if (x>=0 && x < 0.001){
          return('***')
        }
        else if (x>=0.001 && x <0.01){
          return('**')
        }
        else if (x>= 0.01 && x <0.05){
          return('*')
        }
        else if (x>=0.05 && x < 0.1){
          return('.')
        }
        else{
          return(' ')
        }
      }
      cat("Call:\n")
      cat("linreg(formula = ", format(formula), ", data = ", df_name, ")", sep = "")
      cat("\nCoefficients:\n")
      coef = unname(round(reg.coef,5))
      row_names = rownames(reg.coef)

      stdev = round(as.vector(sqrt(unname(var_beta))), 5)

      t_vals =  round(as.vector(unname(t_beta)),2)

      p_value =  2*pt(t_vals, df = df, lower.tail = FALSE)

      p_value_round = p_value

      p_value_round = formatC(p_value_round, format = 'e', digits = 2)

      significance = sapply(p_value, my.sign)

      summary_df = data.frame(cbind(coef, stdev,t_vals, p_value_round, significance), row.names = row_names)

      colnames(summary_df) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", 'Significance')

      print.data.frame(summary_df)
      cat(paste("\nResidual standard error: ", round(sqrt(res_var), 4), " on ", df," degrees of freedom", sep = ""))
    }


  )

)
