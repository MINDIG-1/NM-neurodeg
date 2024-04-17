library(splitTools)

model_fit_metrics <- function(observed, predicted) {

    # EXPV
    mean_observed <- mean(observed)
    expv <- 1 - sum((observed - predicted)^2) / sum((observed - mean_observed)^2)

    # MSLL
    epsilon <- 1e-10  # Adjust this value as needed
    log_loss <- -mean(log(pmax(predicted, epsilon)) * observed + log(pmax(1 - predicted, epsilon)) * (1 - observed))
    msll <- log_loss / log(0.5) # Standardize to [0, 1]

    # SMSE
    mse <- mean((observed - predicted)^2)
    baseline_mse <- mean((observed - mean(observed))^2)
    smse <- mse / baseline_mse  # Standardize to [0, 1]

    # Create a list to store the outputs
    outputs <- list(expv = expv, msll = msll, smse = smse)
    return(outputs)

}

data_sub_norm <- function (df, col) {
    df <- dplyr::select(df, global_id,age, sex, dataset,group, all_of(col))
    names(df) <- c("global_id", "age", "sex", "site", "group", "y")
    # df$y <- (df$y - mean(df$y)) / sd(df$y)

    return (df)
}


df_split <- function(df, strata_columns, p)
{
    ir <- df[strata_columns]
    y <- multi_strata(ir)
    inds <- partition(y, p = c(keep = p, notkeep = 1-p), split_into_list = FALSE)
    df_keep <- df[inds=='keep', ]
    df_notkeep <- df[inds=='notkeep', ]
    return(list(df1 = df_keep, df2 = df_notkeep))

}


gZscore <- function (obj, newdata)
{
  predData <- newdata[, -which(names(newdata) == "y")]
  # mu <- predict(obj,what = "mu", newdata = predData, type = "response")
  # sigma <- predict(obj,what = "sigma", newdata = predData, type="response")
  # nu <- predict(obj,what = "nu", newdata = predData, type="response")
  # tau <- predict(obj,what = "tau", newdata = predData, type="response")
if ("mu"%in%obj$parameters )
    {if ( is.null(obj$mu.fix))
      mu <- predict(obj,what = "mu",   newdata = predData, type = "response")
  else if (obj$mu.fix==TRUE) mu <- rep(fitted(obj, "mu")[1], length(xvar))
  }
if ("sigma"%in%obj$parameters)
   {if (is.null(obj$sigma.fix))
      sigma <- predict(obj,what = "sigma",   newdata = predData, type = "response")
  else if (obj$sigma.fix==TRUE) sigma <- rep(fitted(obj, "sigma")[1], length(xvar))
  }
if ("nu"%in%obj$parameters )
   { if  (is.null(obj$nu.fix))
      nu <- predict(obj,what = "nu",   newdata = predData, type = "response")
  else if (obj$nu.fix==TRUE) nu <- rep(fitted(obj, "nu")[1], length(xvar))
   }
if ("tau"%in%obj$parameters )
   { if (is.null(obj$tau.fix))
      tau <- predict(obj,what = "tau",   newdata = predData, type = "response")
  else if (obj$tau.fix==TRUE) tau <- rep(fitted(obj, "tau")[1], length(xvar))
  }


  yval <- newdata$y
  lpar <- length(obj$parameters)
  fname <- obj$family[1]
  pfun <- paste("p",fname,sep="")

  if(lpar==1)
          {newcall <-call(pfun,yval,mu=mu) }
  else if(lpar==2)
          {newcall <-call(pfun,yval,mu=mu,sigma=sigma) }
  else if(lpar==3)
          {newcall <-call(pfun,yval,mu=mu,sigma=sigma,nu=nu) }
  else
          {newcall <-call(pfun,yval,mu=mu,sigma=sigma,nu=nu,tau=tau) }

  cdf <- eval(newcall)
  rqres <- qnorm(cdf)
  return (rqres)
}
