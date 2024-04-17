vars <- ls()
rm(list = vars)


{
library("gamlss")
library("ggplot2")
library("dplyr")
library("caret")
library("caTools")
source("/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/utils.R")
}

#---------------------------------------------------Paths
{

df_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Featuresdf/df_rel_pwr_avgch_allsites_EC.csv"
result_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Zs"
score_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/CovEffect"
all_ch_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/AllchAllsites"

}


f_bands <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
columns_to_plot <- c('pwr_rel_avg_0', 'pwr_rel_avg_1', 'pwr_rel_avg_2', 'pwr_rel_avg_3', 'pwr_rel_avg_4')
n_bands <- 5

param_files <- list.files(score_path)
csv_files <- param_files[grep("^cov_result", param_files)]

avg_bool <- FALSE
all_channel <- TRUE
nch <- 19

if (avg_bool){

  #-------------------------------------------------------- data :: Train and Test :: Stratify on Sex
    data_all <- read.csv(df_path)
    data_all <- dplyr::select(data_all, global_id, age, sex, dataset,group, all_of(columns_to_plot ))
    data_all <- subset(data_all, age >= 40)
    
    data_HC_all <- subset(data_all, group == "HC")
    data_PD_all <- subset(data_all, group %in% c("PD", "PD Dementia", "PD without cognitive impairment", "PD-MCI"))
    data_PD_all <- mutate(data_PD_all, group = "PD")
    data_AD_all <- subset(data_all, group %in% c("AD", "At Risk", "AD Dementia", "a-MCI", "Alzheimer"))
    data_AD_all <- mutate(data_AD_all, group = "AD")

    strata_columns_HC <- c("age", "sex", "dataset")
    data_HC <- df_split(data_HC_all,strata_columns_HC, 0.8)
    data_HC_tr <- data_HC$df1
    data_HC_te <- data_HC$df2
    data_PD <- data_PD_all
    data_AD <- data_AD_all

    # Open a PDF device for saving the plot
  pdf(paste(result_path,"histogram_density_plot.pdf", sep="/"))

ll <- list()
  for (i in seq_along(columns_to_plot))
  {   #-------- Model Params
      col <-  columns_to_plot[i]
      params <- read.csv(paste(score_path, csv_files[grep(col, csv_files)], sep = "/"))
      mu_formula <- as.formula(params$mu_formula)
      si_formula <- as.formula(params$si_formula)
      nu_formula <- as.formula(params$nu_formula)
      tau_formula <- as.formula(params$nu_formula)
      DistFam <- params$DistFam

      data_tr <- data_sub_norm(data_HC_tr, col)
      data_tr$group <- "HC_tr"
      data_te <- data_sub_norm(data_HC_te, col)
      data_te$group <- "HC_te"
      data_te_PD <- data_sub_norm(data_PD, col)
      data_te_AD <- data_sub_norm(data_AD, col)

      data_tr = na.omit(data_tr)
      m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, tau.formula = tau_formula, family = DistFam, data = data_tr, n.cyc =50)

      zs_tr <- gZscore(m, data_tr)
      zs_te <- gZscore(m, data_te)
      zs_PD <- gZscore(m,data_te_PD)
      zs_AD <- gZscore(m,data_te_AD)
      
      data_tr$zs <- zs_tr
      data_te$zs <- zs_te
      data_te_PD$zs <- zs_PD
      data_te_AD$zs <- zs_AD

    combined_data <- bind_rows(data_tr, data_te, data_te_PD, data_te_AD)

    zs_avg <- combined_data[, c('zs', 'group')]
    zs_avg_ED <- zs_avg[zs_avg$zs > 2 | zs_avg$zs < -2, ]
    ll[[i]] <- table(zs_avg_ED$group)*100/table(zs_avg$group)
        
    desired_order <- c("HC_tr", "HC_te", "PD", "AD")
    combined_data$group <- factor(combined_data$group, levels = desired_order)
    p <- ggplot(combined_data, aes(x = zs, fill = group)) +
          geom_density(alpha = 0.5) +
          labs(title = paste0("Z-Scores Density Plot-Average Relative Power of ", f_bands[i]),
              x = "Z-Score", y = "Density") +
          scale_fill_manual(values = c("HC_tr" = "blue",
                                      "HC_te" = "red",
                                      "PD" = "green",
                                      "AD" = "purple")) +
          theme_light() +
          facet_wrap(~ group, ncol = 1)

    p_name <- col
    ggsave(file.path(result_path, paste0(p_name, ".png")), plot = p)

    zs_save <- paste(result_path,paste0("Avg-",f_bands[i]), sep='/')
    if (!dir.exists(zs_save))
    {
      dir.create(zs_save)
    }
    write.csv(combined_data, paste(zs_save,"zsAVG.csv", sep='/'))

    # Create the initial plot canvas
    ks_tr <- ks.test(zs_tr, "pnorm")
    ks_te <- ks.test(zs_te, "pnorm")
    shap_tr <- shapiro.test(zs_tr)
    shap_te <- shapiro.test(zs_te)

    # Add the histogram
    hist(zs_tr, breaks = 20, freq = FALSE)
    title <- paste("Zscore histogram of Training data of", f_bands[i])
    mtext(paste("Ks_tr: ",ks_tr$p.value," shap_tr: ",shap_tr$p.value), outer = TRUE, line = -1)

    hist(zs_te, breaks = 20, freq = FALSE)
    title <- paste("Zscore histogram of Training data of", f_bands[i])
    mtext(paste("Ks_te: ",ks_te$p.value," shap_te: ",shap_te$p.value), outer = TRUE, line = -1)

  }

  dev.off()
  combined_table <- do.call(rbind, ll)
  df_combined <- as.data.frame(combined_table)
  df_combined['fb'] <- f_bands
  print(df_combined)

}

dataset_list <- list.files(all_ch_path)
col_names_to_drop <- c("subject_id", "handedness", "scores")

n_s <- 1
n_e <- n_bands
if (all_channel)
{
  column_names <- paste0("X", 1:nch)
 
  for (i in n_s:n_e)
  {
    gDev_df_list <- list()

    df <- read.csv(paste(all_ch_path, dataset_list[i], sep = '/'))  #2310
    df <- df[, -which(names(df) %in% col_names_to_drop)]
    df <- subset(df, age > 40)
    data_HC <- subset(df, group == "HC")
    data_PD <- subset(df, group %in% c("PD", "PD Dementia", "PD without cognitive impairment", "PD-MCI"))
    data_PD <- mutate(data_PD, group = "PD")
    data_AD <- subset(df, group %in% c("AD", "At Risk", "AD Dementia", "a-MCI", "Alzheimer"))
    data_AD <- mutate(data_AD, group = "AD")


    strata_columns <- c("age", "sex", "dataset")
    data_HC <- df_split(data_HC,strata_columns, 0.8)
    data_HC_tr <- data_HC$df1
    data_HC_te <- data_HC$df2


    params <- read.csv(paste(score_path, csv_files[i], sep = "/"))
    mu_formula <- as.formula(params$mu_formula)
    si_formula <- as.formula(params$si_formula)
    nu_formula <- as.formula(params$nu_formula)
    tau_formula <- as.formula(params$nu_formula)
    DistFam <- params$DistFam

    for (j in seq_along(column_names)){
      col <- column_names[j]
      data_tr <- data_sub_norm(data_HC_tr, col)
      data_te <- data_sub_norm(data_HC_te, col)
      data_te_PD <- data_sub_norm(data_PD, col)
      data_te_AD <- data_sub_norm(data_AD ,col)

tryCatch({

  data_tr = na.omit(data_tr)
  m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = DistFam, data = data_tr, n.cyc =50)

  }, error = function(e) {

  data_tr = na.omit(data_tr)
  m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula, family = "ST4", data = data_tr, n.cyc =50)
})

temp_df <- data.frame(col = col, score = m$G.deviance, fb = f_bands[i], DistFam = m$family[1] )
gDev_df_list[[j]] <- temp_df


      zs_tr <- gZscore(m,data_tr)
      zs_te <- gZscore(m,data_te)
      zs_PD <- gZscore(m,data_te_PD)
      zs_AD <- gZscore(m,data_te_AD)
      
      col_name <- paste0("zs_", col)
      data_HC_tr[col_name] <- zs_tr
      data_HC_tr$group <- "HC_tr"
      data_HC_te[col_name] <- zs_te
      data_HC_te$group <- "HC_te"
      data_PD[col_name] <-zs_PD
      data_AD[col_name] <-zs_AD

    }
    gDev_df <- do.call(rbind, gDev_df_list)
    combined_data <- rbind(data_HC_tr, data_HC_te, data_PD, data_AD)
    

    zs_save <- paste(result_path,paste0("Allch-",f_bands[i]), sep = "/")
    if (!dir.exists(zs_save))
    {
      dir.create(zs_save)
    }
    write.csv(combined_data, paste(zs_save,"zsAllCh.csv", sep='/'))

    write.csv(gDev_df, paste(zs_save,"globalDev.csv", sep='/'))

  }


  }
  
