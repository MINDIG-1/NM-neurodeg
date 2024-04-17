vars <- ls()
rm(list = vars)

library("gamlss")
library(ggplot2)
library("dplyr")
library(caret)
source("/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/utils.R")


#--------------------------------------------------Functions
Cov_test <- function(m_formula, si_formula, nu_formula, tau_formula, DistFam, cov, data_tr, data_te, result_path, fb){
  formulas <- c("mu_formula", "si_formula", "nu_formula", "tau_formula")
  i <- 1
  data_tr = na.omit(data_tr)
  m0 <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, tau.formula = tau_formula, family = DistFam, data = data_tr) #, n.cyc = 50)
  metrics_test <- model_fit_metrics(data_te$y, m0$mu.fv)

  pfname <- paste(result_path, paste("Model", i-1,cov,fb, "fit.png", sep = "_"), sep = "/")
  png(filename = pfname, width = 10, height = 8, units = 'in', res = 300)
  plot(m0)
  metrics_text <- sprintf(", AIC=%.4f, BIC=%.4f, EXPV=%.4f, MSLL=%.4f, SMSE=%.4f",
            m0$aic, m0$sbc, metrics_test$expv, metrics_test$msll, metrics_test$smse)
  mtext(paste("Model fit metrics,", fb, metrics_text), outer = TRUE, line = -1)
  dev.off()

  score_df <- data.frame(m_i = i-1, AIC_score = m0$aic, BIC_score = m0$sbc, mu_formula = deparse(mu_formula), si_formula=deparse(si_formula), nu_formula=deparse(nu_formula), tau_formula = deparse(tau_formula), DistFam=DistFam,
                         expv_test= metrics_test$expv, msll_test=metrics_test$msll, smse_test=metrics_test$smse )
  
  for (form in formulas){
  print(form)
  tryCatch({
  eval(parse(text = paste(form, "<- update(", form, ", ~ . + ",cov,")")))
  m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, tau.formula = tau_formula, family = DistFam, data = data_tr, n.cyc = 50)
  
  metrics_test <- model_fit_metrics(data_te$y, m$mu.fv)
  if(model$converged){
  temp_df <- data.frame(m_i= i, AIC_score = m$aic, BIC_score = m$sbc, mu_formula = deparse(mu_formula), si_formula=deparse(si_formula), nu_formula=deparse(nu_formula), tau_formula = deparse(tau_formula), DistFam=DistFam,
                        expv_test= metrics_test$expv, msll_test=metrics_test$msll, smse_test=metrics_test$smse )

  score_df <- rbind(score_df,temp_df)
  }
  

  pfname <- paste(result_path, paste("Model", i,cov,fb, "fit.png", sep = "_"), sep = "/")
  png(filename = pfname, width = 10, height = 8, units = 'in', res = 300)
  plot(m)
  metrics_text <- sprintf(", AIC=%.4f, BIC=%.4f, EXPV=%.4f, MSLL=%.4f, SMSE=%.4f",
            m$aic, m$sbc, metrics_test$expv, metrics_test$msll, metrics_test$smse)
  mtext(paste("Model fit metrics,", fb, metrics_text), outer = TRUE, line = -1)
  dev.off()
  }, error = function(e) {
    # Handling error for Code block 1
    warning(paste("Error in iteration", i," in ",cov, ":", conditionMessage(e)))
    temp_df <- data.frame(m_i= 1e10, AIC_score = 1e10, BIC_score = 1e10, mu_formula = deparse(mu_formula), si_formula=deparse(si_formula), nu_formula=deparse(nu_formula), tau_formula = deparse(tau_formula), DistFam=DistFam,
                        expv_test= 0, msll_test=0, smse_test=0 )
    if (length(dev.list()) > 0) dev.off()

  })
  i <- i + 1
  
}
return(score_df)

}

#--------------------------------------------------Paths
df_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Featuresdf/df_rel_pwr_avgch_allsites_EC.csv"
score_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Distnpoly"
result_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/CovEffect"

data_all <- read.csv(df_path)
data_all <- subset(data_all, age >= 40)

data_HC_all <- subset(data_all, group == "HC")

strata_columns <- c("age", "sex", "dataset")

data_HC<- df_split(data_HC_all,strata_columns, 0.8)
data_HC_tr <- data_HC$df1
data_HC_te <- data_HC$df2

f_bands <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
columns_to_plot <- c('pwr_rel_avg_0', 'pwr_rel_avg_1', 'pwr_rel_avg_2', 'pwr_rel_avg_3', 'pwr_rel_avg_4')

dist_list <- c("PE","GT","GG","GB2","BCT","exGAUS","JSU","SEP2","SEP3","SEP4","SHASH","SHASHo","ST3","ST4")

param_files <- list.files(score_path)
bic_files <- param_files[grep("^BIC", param_files)]
bic_files <- bic_files[-1]

phase2 <- TRUE

for (i in seq_along(columns_to_plot)){

  print(f_bands[i])
 
  cov_results_df <- data.frame()
  cov_sex_df <- data.frame()

  col <- columns_to_plot[i]
  
  params <- read.csv(paste(score_path, bic_files[grep(col, bic_files)], sep="/"))
  params <- params[1,]
  mu_formula <- as.formula(params$mu_formula)
  si_formula <- as.formula(params$si_formula)
  nu_formula <- as.formula(params$nu_formula)
  tau_formula <- as.formula("~1")
  DistFam <- params$DistFamily
  formulas <- c("mu_formula", "si_formula", "nu_formula", "tau_formula")


  data_tr <- data_sub_norm(data_HC_tr, col)
  data_te <- data_sub_norm(data_HC_te, col)

  data_tr2 <- data_tr
  data_tr2['group'] <- "HC_tr"
  data_te2 <- data_te
  data_te2['group'] <- "HC_te"

  data_temp <- rbind(data_tr2, data_te2)
  p_data_c <- ggplot(data = data_temp, aes(x=age, y=y, color = factor(group))) + 
                geom_point()

  cpname <- paste(result_path, paste0("Datascatter",f_bands[i], ".png"), sep = "/")
  ggsave(filename = cpname, plot = p_data_c, width = 10, height = 8)

  #------------------------------------------------Covariate :: Sex :: Fixed
  cov <- "factor(sex)"
  aic_sex <- Cov_test(mu_formula,si_formula, nu_formula, tau_formula,DistFam, cov, data_tr, data_te, result_path, f_bands[i])
  sorted_aic_sex <- aic_sex[order(aic_sex$AIC_score),]
  sorted_bic_sex <- aic_sex[order(aic_sex$BIC_score),]

  sorted_bic_sex <- sorted_bic_sex [1,]

  cov_sex_df <- rbind(cov_sex_df,sorted_bic_sex)
  cov_file_name <- paste0("cov_sex_",col,".csv")
  write.csv(cov_sex_df, paste(result_path, cov_file_name, sep="/"))


  #------------------------------------------------Covariate :: Site :: Fixed
  mu_formula <- as.formula(sorted_bic_sex$mu_formula)
  si_formula <- as.formula(sorted_bic_sex$si_formula)
  nu_formula <- as.formula(sorted_bic_sex$nu_formula)
  tau_formula <- as.formula(sorted_bic_sex$tau_formula)

  cov <- "factor(site)"
  aic_site <- Cov_test(mu_formula,si_formula, nu_formula, tau_formula,DistFam, cov, data_tr, data_te, result_path, f_bands[i])
  sorted_aic_site <- aic_site[order(aic_site$AIC_score),]
  sorted_bic_site <- aic_site[order(aic_site$BIC_score),]
  sorted_bic_site <- sorted_bic_site [1,]


  # ------------------------------------------------Covariate :: Site :: Random

  mu_formula <- as.formula(sorted_bic_site$mu_formula)
  si_formula <- as.formula(sorted_bic_site$si_formula)
  nu_formula <- as.formula(sorted_bic_site$nu_formula)
  tau_formula <- as.formula(sorted_bic_site$tau_formula)

  cov <- "random(factor(site))"
  aic_r_site <- Cov_test(mu_formula,si_formula, nu_formula, tau_formula,DistFam, cov, data_tr, data_te, result_path, f_bands[i])
  sorted_aic_r_site <- aic_r_site[order(aic_r_site$AIC_score),]
  sorted_bic_r_site <- aic_r_site[order(aic_r_site$BIC_score),]
  sorted_bic_r_site <- sorted_bic_r_site [1,]


  #-------------------------------------------------Final model
  cov_results_df <- rbind(cov_results_df,sorted_bic_r_site)

  if (phase2)
  {

    mu_formula <- as.formula(cov_results_df$mu_formula)
    si_formula <- as.formula(cov_results_df$si_formula)
    nu_formula <- as.formula(cov_results_df$nu_formula)
    tau_formula <- as.formula(cov_results_df$tau_formula)
tryCatch ({
    for (d in dist_list)
    {
      data_tr = na.omit(data_tr)
      model <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula,family = d, data = data_tr, n.cyc = 50)
      if (model$sbc < cov_results_df$BIC_score && model$converged)
      {
        cov_results_df$DistFam <- model$family[1]
      }
    }
  }, error = function(e) {

          })
    cov_file_name <- paste0("cov_result_",col,".csv")
  write.csv(cov_results_df, paste(result_path, cov_file_name, sep="/"))

  }

}
