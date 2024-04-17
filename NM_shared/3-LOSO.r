vars <- ls()
rm(list = vars)

#---------------------------------------------------libraries
{
  library(ggplot2)
  library(dplyr)
  library("gamlss")
  library(caret)

  source("/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/utils.R")
}

#---------------------------------------------------Paths
{

df_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Featuresdf/df_rel_pwr_avgch_allsites_EC.csv"
site_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Featuresdf"
plot_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/LOSO"
score_path <- "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/CovEffect"

}

# --------------------------------------------------Data and Results Paths

data_all <- read.csv(df_path)  # HC data
data_all <- subset(data_all, age >= 40)

data_HC_all <- subset(data_all, group == "HC")

min_age <- min(data_HC_all$age)
max_age <- max(data_HC_all$age)

n_s <- 200

#-------------------------------------------------Variables
{
  columns_to_plot <- c('pwr_rel_avg_0', 'pwr_rel_avg_1', 'pwr_rel_avg_2', 'pwr_rel_avg_3', 'pwr_rel_avg_4')
  f_bands <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")

  site_names <- c('andreas_miltiadous', 'avila', 'basel', 'cavanagh', 'CreaPark', 'gorsev', 'istanbul', 'leb2017', 'leb2019', 'marseille', 'railo', 'schalk', 'SRM', 'zenodoFuglsang')
  f <- 'df_rel_pwr_avgch_info_EC.csv'
  file_names <- paste(site_names, f, sep = "/")

  distinct_colors <- rainbow(length(site_names)+1)

param_files <- list.files(score_path)
aic_files <- param_files[grep("^cov_sex", param_files)]

}
corr_df <- data.frame(fb = character(), df_omitted = numeric(), correlation = numeric())

for (j in seq_along(columns_to_plot)){

    col <- columns_to_plot[j]
    params <- read.csv(paste(score_path, aic_files[grep(col, aic_files)], sep="/"))
    params <- params[1,]

    data_all_copy <-  data_HC_all
    data_all_copy <- data_sub_norm(data_all_copy, col)
    p <- ggplot(data = data_all_copy, aes(x = age, y = y, color = factor(site)), alpha = 0.4) + geom_point() + geom_point(size = 3) + scale_colour_manual(values = distinct_colors) +
    theme(plot.background = element_rect(fill = "white"))
    m_mu <-  data.frame()
    data_all_copy = na.omit(data_all_copy)
    m_all <- gamlss(as.formula(params$mu_formula), sigma.formula = as.formula(params$si_formula), nu.formula = as.formula(params$nu_formula), family = params$DistFamily, data = data_all_copy, n.cyc=50)
    
    Qua50 <- getQuantile(m_all, term="age", quantile = 0.5, n.points = n_s, plot = FALSE)
    c<-curve(Qua50, min_age, max_age,  n = n_s)
    mu_all <- c$y
    p <- p + geom_line(data=data.frame(x = c$x, y = mu_all, site = "All"), aes(x=x, y=y)) + theme(plot.background = element_rect(fill = "white"))

    # In each loop, one site is discarded
    for (i in 1:length(file_names)) {

        combined_data <- data.frame()
        # Read all files except the current one
        files_to_read <- file.path(site_path, file_names[-i])

        for (file in files_to_read){
          current_data <- read.csv(file)
          current_data <- subset(current_data, age >= 40)
          current_data <- subset(current_data, group == 'HC')
          combined_data <- rbind(combined_data, data_sub_norm(current_data,col))
        }

        # Read and combine the remaining files
        current_data_tr <- combined_data
        #Training the gamlss model
        current_data_tr = na.omit(current_data_tr)
        m <- gamlss(as.formula(params$mu_formula), sigma.formula = as.formula(params$si_formula), nu.formula = as.formula(params$nu_formula), family = params$DistFamily, data = current_data_tr, n.cyc=50)

        Qua50 <- getQuantile(m, term="age", quantile = 0.5, n.points = n_s, plot = FALSE)
         c<-curve(Qua50, min_age, max_age,  n = n_s)
         mu_site <- c$y

        correlation <- cor(mu_all, mu_site)
        corr_df_temp <- data.frame(fb = col, df_omitted=site_names[i], correlation= correlation)
        corr_df <- rbind(corr_df, corr_df_temp)

        # plotting the fitted mu on top of the scatter plot of all HC data
        m_mu <- data.frame(x = c$x, y = mu_site, site =site_names[i] )
        p <- p + geom_line(data = m_mu, aes(x=x,y=y, color = factor(site)), linewidth = 1, alpha = 0.7 ) + theme(plot.background = element_rect(fill = "white"))

    }

datamu<- paste(plot_path, paste0("MeanTrack",f_bands[j],".png"), sep="/")
print(p )
ggsave(filename = datamu, plot = p, width = 10, height = 8)
}

write.csv(corr_df, paste(plot_path, "corr_sites.csv", sep="/"))

