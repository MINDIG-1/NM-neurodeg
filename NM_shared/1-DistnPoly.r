vars <- ls()
rm(list = vars)

library("gamlss")
library("ggplot2")
library("dplyr")
source("/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/utils.R")


df_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Featuresdf/df_rel_pwr_avgch_allsites_EC.csv"
result_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Distnpoly/"

data_all <- read.csv(df_path)
data_all <- subset(data_all, age >= 40)
data_HC_all <- subset(data_all, group == "HC")

columns_to_plot <- c('pwr_rel_avg_0', 'pwr_rel_avg_1', 'pwr_rel_avg_2', 'pwr_rel_avg_3', 'pwr_rel_avg_4')

f_bands <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")

dist_train <- TRUE

dist_list <- c("PE","GT","GG","GB2","BCT","exGAUS","JSU","SEP2","SEP3","SEP4","SHASH","SHASHo","ST3","ST4")
n_c <- 50

FP.SET <- matrix(c(1,0,0,
                   1,1,0,
                   1,1,1,
                   2,1,0),
                   byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","nu")))

df_all <- data.frame(Feature = character(), iFP = character(), DistFamily = character(),
                     AIC_Score = numeric(), BIC_Score = numeric(),  mu_formula=character(), si_formula= character(), nu_formula=character(), stringsAsFactors = FALSE)

if (dist_train)
{

# Loop over each band, family, additive term, dof (I used M3 format)
for (column in columns_to_plot){

  df <- data.frame(Feature = character(), iFP = character(), DistFamily = character(),
                   AIC_Score = numeric(), BIC_Score = numeric(), mu_formula=character(), si_formula= character(), nu_formula=character(),stringsAsFactors = FALSE)

  for (d in dist_list){


      for(iFP in 1:NROW(FP.SET) ) {


        if(FP.SET[iFP,"mu"]>0){
          mu_formula <- as.formula(paste("y ~ 1 +fp(age, ", FP.SET[iFP,"mu"], ")"))
        }else{
          mu_formula <- as.formula(paste("y ~ 1 "))
        }
        if(FP.SET[iFP,"sigma"]>0){
          sigma_formula <- as.formula(paste("~ 1 +fp(age, ", FP.SET[iFP,"sigma"], ")"))
        }else{
          sigma_formula <- as.formula(paste("~ 1 "))
        }
        if(FP.SET[iFP,"nu"]>0){
          nu_formula <- as.formula(paste("~ 1 + fp(age, ", FP.SET[iFP,"nu"], ")"))
        }else{
          nu_formula <- as.formula(paste("~ 1"))
        }

          avgch <- data_sub_norm(data_HC_all, column)

          tryCatch({
            print(d)
            avgch = na.omit(avgch)
            model <- gamlss(mu_formula, sigma.formula = sigma_formula, nu.formula = nu_formula,family = d, data = avgch, n.cyc = n_c)
            if (model$converged){
            temp_df <- data.frame(Feature = column, iFP=paste(as.character(FP.SET[iFP,]), collapse = ""),
                                DistFamily = d, AIC_Score = model$aic, BIC_Score = model$sbc,mdlDev = model$G.deviance,
                                mu_formula=deparse(mu_formula), si_formula= deparse(sigma_formula), nu_formula=deparse(nu_formula),
                                stringsAsFactors = FALSE)
            } else{
              temp_df <- data.frame(Feature = column, iFP=paste(as.character(FP.SET[iFP,]), collapse = ""),
                                DistFamily = paste0(d," - Not Converged"), AIC_Score = numeric(), BIC_Score = numeric() ,mdlDev = model$G.deviance,
                                mu_formula=deparse(mu_formula), si_formula= deparse(sigma_formula), nu_formula=deparse(nu_formula),
                                stringsAsFactors = FALSE)
            }                    
            df <- rbind(df, temp_df)
            df_all <- rbind(df_all, temp_df)
          }, error = function(e) {
            error_list <- conditionMessage(e)
          })
        }
    }
  
  # csv of sorted BIC with params for each band separatly
  sorted_BIC <- df[order(df$BIC_Score), ]
  write.csv(sorted_BIC, paste0(result_path,"BIC_scores_", column, ".csv"))
  sorted_AIC <- df[order(df$AIC_Score), ]
  write.csv(sorted_AIC, paste0(result_path,"AIC_scores_", column, ".csv"))

}

# csv of sorted BIC with params for all bands
sorted_BIC_all <- df_all[order(df_all$BIC_Score), ]
write.csv(sorted_BIC_all, paste0(result_path,"BIC_scores_all.csv"))

df <- sorted_BIC_all
} else {
  df <- read.csv(paste0(result_path,"BIC_scores_all.csv"))
}
feature_mapping <- c("pwr_rel_avg_0" = "Average Relative Power-delta",
                    "pwr_rel_avg_1" = "Average Relative Power-theta",
                    "pwr_rel_avg_2" = "Average Relative Power-alpha",
                    "pwr_rel_avg_3" = "Average Relative Power-beta",
                    "pwr_rel_avg_4" = "Average Relative Power-gamma")


# Use mutate and gsub to replace feature values
df <- df %>%
  mutate(Feature = recode(Feature, !!!feature_mapping))


min_bic <- df %>%
  group_by(Feature) %>%
  summarize(min_BIC_Score = min(BIC_Score))


# Merge the minimum BIC_Score back to the original data
df <- left_join(df, min_bic, by = "Feature")

  df <- df %>%
  mutate(BIC_Score = BIC_Score - min_BIC_Score)


df <- df[df$BIC_Score < 1000, ]

# ggplot code
p <- ggplot(df, aes(DistFamily, BIC_Score)) +
  geom_bar(aes(fill = DistFamily), stat = "identity") +
  facet_wrap(~Feature, scales = "free_y", ncol = 1)+
  labs(x = "DistFamily", y = "BIC_Score") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

# Calculate the minimum BIC_Score for each facet
min_bic <- df %>%
  group_by(Feature) %>%
  slice(which.min(BIC_Score))

# Add geom_text to mark the minimum BIC_Score
p + geom_text(data = min_bic, aes(label = "Min"), vjust = -0.5)

distplot <- paste(result_path, "DistTest.png", sep="/")
ggsave(filename = distplot, plot = p, width = 20, height = 16)

print(p)
