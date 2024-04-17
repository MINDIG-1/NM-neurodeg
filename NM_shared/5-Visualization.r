vars <- ls()
rm(list = vars)


{
library("gamlss")
library(ggplot2)
library("dplyr")
library("caret")
library(caTools)
source("/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/utils.R")
require(GGally)
}

#---------------------------------------------------Functions
newDataGen <- function(samplesize, data) {
    original_data <- data[, c("sex", "site")]
    original_data$sex <- factor(original_data$sex)
    original_data$site <- factor(original_data$site)

    proportions <- table(original_data$sex, original_data$site) / nrow(original_data)

    # Create an empty data frame for the new dataset
    new_data <- data.frame(sex = character(0), site = character(0))

    for (sex_value in levels(original_data$sex)) {
        for (site_value in levels(original_data$site)) {
            n_to_sample <- round(samplesize * proportions[sex_value, site_value])

            subset_data <- subset(original_data, sex == sex_value & site == site_value)

            # Make sure n_to_sample is not greater than the available rows
            n_to_sample <- min(n_to_sample, nrow(subset_data))

            sampled_indices <- sample(nrow(subset_data), n_to_sample, replace = TRUE)
            new_data <- rbind(new_data, subset_data[sampled_indices, ])
        }
    }

    return(new_data)
}

#---------------------------------------------------Paths
{

df_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Featuresdf/df_rel_pwr_avgch_allsites_EC.csv"
result_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Visualization"
score_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/CovEffect"
zs_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Zs"

}

f_bands <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
columns_to_plot <- c('pwr_rel_avg_0', 'pwr_rel_avg_1', 'pwr_rel_avg_2', 'pwr_rel_avg_3', 'pwr_rel_avg_4')

param_files <- list.files(score_path)
csv_files <- param_files[grep("\\.csv$", param_files)]
vmax <- 1.1
vmin <- 0
txt_size <- 32


#-------------------------------------------------------- data :: Train and Test :: Stratify on Sex
data_all <- read.csv(df_path)
data_all <- subset(data_all, age >= 40)

# Specify the old and new group names
old_PD_names <- c("PD", "PD Dementia", "PD without cognitive impairment", "PD-MCI")
old_AD_names <- c("AD", "At Risk", "AD Dementia", "a-MCI", "Alzheimer")
new_PD_names <- "PD"
new_AD_names <- "AD"

data_all <- data_all %>%
  mutate(group = ifelse(group %in% old_PD_names, new_PD_names, group))
data_all <- data_all %>%
  mutate(group = ifelse(group %in% old_AD_names, new_AD_names, group))

print(data_all$group)

groups_tokeep <- c("HC", "AD", "PD")
data_all <- data_all[data_all$group %in% groups_tokeep, ]

data_HC_all <- subset(data_all, group == "HC")

strata_columns <- c("age", "sex", "dataset")
data_HC <- df_split(data_HC_all,strata_columns, 0.8)
data_HC_tr <- data_HC$df1
data_HC_te <- data_HC$df2

min_age <- min(data_HC_all$age)
max_age <- max(data_HC_all$age)

avg_zs_files <- list.files(zs_path)
avg_zs_files <- avg_zs_files[grep("Avg-", avg_zs_files)]


for (i in seq_along(columns_to_plot))
{   #----------------------------------------------- Model Params
    col <-  columns_to_plot[i]
    params <- read.csv(paste(score_path, csv_files[grep(paste0("result_",col), csv_files)], sep = "/"))
    mu_formula <- as.formula(params$mu_formula)
    si_formula <- as.formula(params$si_formula)
    nu_formula <- as.formula(params$nu_formula)
    tau_formula <- as.formula(params$nu_formula)
    DistFam <- params$DistFam

    avg_zs_file <- avg_zs_files[grep(f_bands[i], avg_zs_files)]
    print(avg_zs_file)
    zs_avg <- read.csv(paste(zs_path, avg_zs_file, "zsAVG.csv" ,sep = "/"))
    zs_avg_HC <- subset(zs_avg, group %in% c("HC_tr", "HC_te"))

    grp_name <- c("HC_tr", "HC_te", "PD", "AD")
    grp_val <- c(nrow(zs_avg[zs_avg$group == 'HC_tr', ]), nrow(zs_avg[zs_avg$group == 'HC_te', ]), nrow(zs_avg[zs_avg$group == 'PD', ]), nrow(zs_avg[zs_avg$group == 'AD', ]))
    ntot_grp <- data.frame(grp_name, grp_val)

    zs_avg <- zs_avg[zs_avg$zs > 2 | zs_avg$zs < -2, ]

    grp_val <- c(nrow(zs_avg[zs_avg$group == 'HC_tr', ]), nrow(zs_avg[zs_avg$group == 'HC_te', ]), nrow(zs_avg[zs_avg$group == 'PD', ]), nrow(zs_avg[zs_avg$group == 'AD', ]))
    ndev_grp <- data.frame(grp_name, grp_val)
    print('*********')
    grp_val <- ndev_grp['grp_val']/ntot_grp['grp_val']*100
    percdev_grp <- data.frame(grp_name, grp_val)
    percdev_file_name <- paste0("percdev_",col,".csv")
    write.csv(percdev_grp, paste(result_path, percdev_file_name, sep="/"))

    # Wilcox test: pval for signif diff of ED between HCte and clinical grp
    result_PD <- wilcox.test(as.numeric(unlist(na.omit(zs_avg[zs_avg$group == 'HC_te', ]['zs']))), as.numeric(unlist(na.omit(zs_avg[zs_avg$group == 'PD', ]['zs']))))
    result_AD <- wilcox.test(as.numeric(unlist(na.omit(zs_avg[zs_avg$group == 'HC_te', ]['zs']))), as.numeric(unlist(na.omit(zs_avg[zs_avg$group == 'AD', ]['zs']))))
    print(result_PD$p.value)
    print(result_AD$p.value)
    pval_ndev <- data.frame(c("PD", "AD"), c(result_PD$p.value, result_AD$p.value))
    write.csv(pval_ndev, paste(result_path,  paste0("pvaldev_",col,".csv"), sep="/"))

    zs_avg_HC <- zs_avg_HC[, !(names(zs_avg_HC) %in% c("zs", "X"))]
    data_all_p <- data_sub_norm(data_all, col)
    data_all_p <- data_all_p[data_all_p$group != "HC", ]
    data_all_p <- rbind(data_all_p,zs_avg_HC)

    data_tr <- data_sub_norm(data_HC_tr, col)
    data_te <- data_sub_norm(data_HC_te, col)

    data_tr = na.omit(data_tr)
    m <- gamlss(mu_formula, sigma.formula = si_formula, nu.formula = nu_formula,tau.formula = tau_formula, family = DistFam, data = data_tr, n.cyc =20)

tryCatch({

    Qua025 <- getQuantile(m, term="age", quantile = 0.025, n.points = 3000, plot = FALSE)
    c<-curve(Qua025, min_age, max_age,  n = 3000)
    c_age <- c$x
    c025 <- c$y
    Qua50 <- getQuantile(m, term="age", quantile = 0.5, n.points = 3000, plot = FALSE)
    c<-curve(Qua50, min_age, max_age,  n = 3000)
    c50 <- c$y
    Qua975 <- getQuantile(m, term="age", quantile = 0.975, n.points = 3000, plot = FALSE)
    c<-curve(Qua975, min_age, max_age,  n = 3000)
    c975 <- c$y

    pred_df <- data.frame(age = c_age, c975=c975, c50=c50, c025=c025)

    group_order <- c("HC_tr", "HC_te", "PD", "AD")
    group_colors <- c("HC_tr" = "dodgerblue", "HC_te" = "red", "PD" = "green", "AD" = "darkorchid")
    p_data_c <- ggplot() + 
                geom_point(data = data_all_p, aes(x=age, y=y, color = factor(group, levels = group_order))) +
                geom_point(data = zs_avg, aes(x = age, y= y , color = factor(group, levels = group_order)), shape = 21, size = 3)+
                geom_line(data = pred_df, aes(x=age, y=c975), color = "blue", linetype = "dashed")+
                geom_line(data = pred_df, aes(x=age, y=c50),  color = "blue" )+
                geom_line(data = pred_df, aes(x=age, y=c025), color = "blue", linetype = "dashed")+ 
                geom_line(data = data.frame(x = data_tr$age, y = m$mu.fv), aes(x=x, y=y))+
                ylim(vmin, vmax)+
                labs(x = "Age", y = paste0("Relative Power - ",f_bands[i])) +
                scale_color_manual(values = group_colors) 
    
    cpname <- paste(result_path, paste0("DatanMeanPlot_",f_bands[i], ".png"), sep = "/")
    ggsave(filename = cpname, plot = p_data_c, width = 10, height = 8)


    data_all_p_HCtrte <- data_all_p[data_all_p$group %in% c("HC_te", "HC_tr", "HC"), ]
    zs_avg_HCtrte <- subset(zs_avg, group %in% c("HC_tr", "HC_te", "HC"))
    group_order <- c("HC_tr", "HC_te")
    group_colors <- c("HC_tr" = "dodgerblue", "HC_te" = "red")
    p_data_c <- ggplot() + 
                geom_point(data = data_all_p_HCtrte, aes(x=age, y=y, color = factor(group, levels = group_order))) +
                geom_point(data = zs_avg_HCtrte, aes(x = age, y= y , color = factor(group, levels = group_order)), shape = 21, size = 3)+
                geom_line(data = pred_df, aes(x=age, y=c975), color = "blue", linetype = "dashed")+
                geom_line(data = pred_df, aes(x=age, y=c50),  color = "blue" )+
                geom_line(data = pred_df, aes(x=age, y=c025), color = "blue", linetype = "dashed")+ 
                geom_line(data = data.frame(x = data_tr$age, y = m$mu.fv), aes(x=x, y=y))+
                ylim(vmin, vmax)+
                labs(x = "Age", y = paste0("Relative Power - ",f_bands[i])) +
                scale_color_manual(values = group_colors) 
    
    cpname <- paste(result_path, paste0("DatanMeanPlot_",f_bands[i], "_HC.png"), sep = "/")
    ggsave(filename = cpname, plot = p_data_c, width = 10, height = 8)

    
    data_all_p_PDAD <- data_all_p[data_all_p$group %in% c("PD", "AD"), ]
    zs_avg_PDAD <- subset(zs_avg, group %in% c("PD", "AD"))
    group_order <- c("PD", "AD")
    group_colors <- c("PD" = "green", "AD" = "darkorchid")
    p_data_c <- ggplot() + 
                geom_point(data = data_all_p_PDAD, aes(x=age, y=y, color = factor(group, levels = group_order))) +
                geom_point(data = zs_avg_PDAD, aes(x = age, y= y , color = factor(group, levels = group_order)), shape = 21, size = 3)+
                geom_line(data = pred_df, aes(x=age, y=c975), color = "blue", linetype = "dashed")+
                geom_line(data = pred_df, aes(x=age, y=c50),  color = "blue" )+
                geom_line(data = pred_df, aes(x=age, y=c025), color = "blue", linetype = "dashed")+ 
                geom_line(data = data.frame(x = data_tr$age, y = m$mu.fv), aes(x=x, y=y))+
                ylim(vmin, vmax)+
                labs(x = "Age", y = paste0("Relative Power - ",f_bands[i])) +
                scale_color_manual(values = group_colors) 
    
    cpname <- paste(result_path, paste0("DatanMeanPlot_",f_bands[i], "_PDAD.png"), sep = "/")
    ggsave(filename = cpname, plot = p_data_c, width = 10, height = 8)

}, error = function(e) {
    print(conditionMessage(e))
    })
    
}
