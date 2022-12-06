## TO DO when it works
# * Try formula variable
# * DONE: prefix in variable
# * starting positions plink and prsice
# * target_prs.lasso <- read.table('IBD_GSA_fin_prs_lasso1.txt', header = T)[,1:3] ==> variable file name

project_dir = "/staging/leuven/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/"
output_data_dir = paste(project_dir, "output_data/", sep="")

# GET THEM FROM SNAKEMAKE
target_prefix <- "IBD_GSA_fin"
external_prefix <- "1000G_EUR_fin"
pval_thr <- "5e-8,1e-5,0.01,0.05,0.1,0.5" #p-values used in PLINK and PRSice

# =========================================================================================================================================================================================

library(dplyr)
library(ggplot2)
library(data.table)
library(ROCR)
library(pROC)
library(pscl)
library(stringr)

# =========================================================================================================================================================================================
setwd(output_data_dir)

plink_thr <- read.table("a_1_julia_plink/target_data/range_list", header=F)[,1]

target_plink_files <- c()
external_plink_files <- c()
for (i in 1:length(plink_thr)) {
    target_plink_files[i] <- paste(target_prefix,'.', plink_thr[i], '.profile', sep="")
    external_plink_files[i] <- paste(external_prefix,'.', plink_thr[i], '.profile', sep="")
}

# =========================================================================================================================================================================================

no_plink_thr = length(plink_thr)
no_prsice_thr = (str_count(pval_thr, ',')) + 1
# no_lasso_thr = 
# no_ldpred_thr =  

# =========================================================================================================================================================================================

setwd(project_dir)

# Phenotypes file
phenotype <- read.table("data/final_phenotypes.txt", header=T)
colnames(phenotype) <- c("FID", "IID","pheno") 

# PC file
pcs_file = paste('data/', target_prefix, '.eigenvec', sep="")
pcs <- read.table(pcs_file, header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 

covariate_file = paste('data/', target_prefix, '.cov', sep="")
covariate <- read.table(covariate_file, header=T)

# =========================================================================================================================================================================================
### PLINK ###
# =========================================================================================================================================================================================

# Get TARGET scores
setwd(paste(output_data_dir,'a_1_julia_plink/target_data',sep=""))

target_prs.plink <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(target_prs.plink) <- c('FID', 'IID', 'pheno')
target_prs.plink <- read.table(target_plink_files[1], header = T)[,1:3]

for (i in 1:length(plink_thr)) {
    target_prs.plink_file <- read.table(target_plink_files[i], header = T)
    target_prs.plink_file <- target_prs.plink_file %>% select(1,2,6)
    score_col_name <- paste('SCORE', plink_thr[i], sep="")
    colnames(target_prs.plink_file) <- c('FID', 'IID', score_col_name)
    target_prs.plink <- merge(target_prs.plink, target_prs.plink_file, by = c("FID", "IID"))    
}    

# Get EXTERNAL scores
setwd(paste(output_data_dir,'a_1_julia_plink/external_data',sep=""))

external_prs.plink <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(external_prs.plink) <- c('FID', 'IID')
external_prs.plink <- read.table(external_plink_files[1], header = T)[,1:2]

for (i in 1:length(plink_thr)) {
    external.prs.plink_file <- read.table(external_plink_files[i], header = T)
    external.prs.plink_file <- external.prs.plink_file %>% select(1,2,6)
    score_col_name_ext <- paste('SCORE', plink_thr[i], sep="")
    colnames(external.prs.plink_file) <- c('FID', 'IID', score_col_name_ext)
    external_prs.plink <- merge(external_prs.plink, external.prs.plink_file, by = c("FID", "IID"))
}    

cat("\nPLINK target and external scores uploaded into single dataframes\n")

no_plink_thr = length(plink_thr)
start_pos_ext = 3
start_pos_target = 4

ext_scores.plink <- external_prs.plink[,start_pos_ext:(start_pos_ext+no_plink_thr-1)]

ext_mean.plink <- as.data.frame(colMeans(ext_scores.plink))
ext_sd.plink <- as.data.frame(apply(ext_scores.plink, 2, sd))

colnames(ext_mean.plink) <- c("ext_mean.plink")
colnames(ext_sd.plink) <- c("ext_sd.plink")

cat("\nPLINK external scores' mean and standard deviation calculated")
cat("\nPLINK starting standardization...")

# Before standardization

std_names <- colnames(target_prs.plink)

for (sc in c(start_pos_target:(start_pos_target+no_plink_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.plink)[sc], "_std", sep="")
}

std_prs.plink <- target_prs.plink
colnames(std_prs.plink) <- std_names

for (sc in c(start_pos_target:(start_pos_target+no_plink_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.plink[(sc-start_pos_target+1),1]
    stdv <- ext_sd.plink[(sc-start_pos_target+1),1]
    std_prs.plink[,sc] <- (target_prs.plink[,sc] - mean) / stdv
}

cat("\nPLINK standardization complete")

std_prs_ph.plink <- std_prs.plink
std_prs_ph.plink$PHENO <- std_prs_ph.plink$PHENO - 1

# vector with the explained variance for each threshold 
explained_var.plink <- c()
for (i in 1:length(plink_thr)) {
    
    score_var <- paste('SCORE', plink_thr[i], '_std', sep="")
    
    log_model.plink <- glm(std_prs_ph.plink$PHENO ~ pull(std_prs_ph.plink, score_var), family = binomial(link = "logit"))
    explained_var.plink[i] <- with(summary(log_model.plink), 1 - deviance/null.deviance)
}

#obtain the highest explained variance 
best_prs.plink <- which.max(explained_var.plink)

#obtain the best prs by indexing from which model had the best explained variance
best_target.plink <- std_prs_ph.plink %>% select(1, 2, 3, (start_pos_target+best_prs.plink-1))
best_target.plink <- best_target.plink[complete.cases(best_target.plink), ]
colnames(best_target.plink) <- c('FID', 'IID', 'pheno', 'std_prs')

#recode phenotype into a factor
best_target.plink$pheno <- as.factor(best_target.plink$pheno)

cat("\nPLINK scores for best threshold saved")

setwd(output_data_dir)

#need to write table for best prs
write.table(best_target.plink, "005_comparison/best_prs_plink.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nPLINK finished")

# =========================================================================================================================================================================================
### PRSICE ###
# =========================================================================================================================================================================================

setwd(output_data_dir)

# Get TARGET scores

target_prs.prsice <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(target_prs.prsice) <- c('FID', 'IID')
target_prs.prsice_file <- read.table("a_2_julia_prsice/target_data/IBD_GSA_fin.all_score", header = T)
target_prs.prsice <- merge(target_prs.prsice_file[,1:2], phenotype, by=c("FID", "IID"))
target_prs.prsice <- merge(target_prs.prsice, target_prs.prsice_file, by = c("FID", "IID"))

# Get EXTERNAL scores
external_prs.prsice <- read.table("a_2_julia_prsice/external_data/1000G_EUR_fin.all_score", header = T)

cat("\nPRSice target and external scores uploaded into single dataframes")

no_prsice_thr = 7
start_pos_ext = 3
start_pos_target = 4

ext_scores.prsice <- external_prs.prsice[,start_pos_ext:(start_pos_ext+no_prsice_thr-1)]

ext_mean.prsice <- as.data.frame(colMeans(ext_scores.prsice))
ext_sd.prsice <- as.data.frame(apply(ext_scores.prsice, 2, sd))

colnames(ext_mean.prsice) <- c("ext_mean.prsice")
colnames(ext_sd.prsice) <- c("ext_sd.prsice")

cat("\nPRSice external scores' mean and standard deviation calculated")
cat("\nPRSice starting standardization...")

# Before standardization

std_names <- colnames(target_prs.prsice)

for (sc in c(start_pos_target:(start_pos_target+no_prsice_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.prsice)[sc], "_std", sep="")
}

std_prs.prsice <- target_prs.prsice
colnames(std_prs.prsice) <- std_names

for (sc in c(start_pos_target:(start_pos_target+no_prsice_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.prsice[(sc-start_pos_target+1),1]
    stdv <- ext_sd.prsice[(sc-start_pos_target+1),1]
    std_prs.prsice[,sc] <- (target_prs.prsice[,sc] - mean) / stdv
}

cat("\nPRSice standardization complete")

std_prs_ph.prsice <- std_prs.prsice
std_prs_ph.prsice$pheno <- std_prs_ph.prsice$pheno - 1

explained_var.prsice <- c()
for (i in 1:no_prsice_thr) {    
    score_var <- colnames(std_prs_ph.prsice)[start_pos_target+i-1]
    
    log_model.prsice <- glm(std_prs_ph.prsice$pheno ~ pull(std_prs_ph.prsice, score_var), family = binomial(link = "logit"))
    explained_var.prsice[i] <- with(summary(log_model.prsice), 1 - deviance/null.deviance)
}

#obtain the highest explained variance 
best_prs.prsice <- which.max(explained_var.prsice)

#obtain the best prs by indexing from which model had the best explained variance
best_target.prsice <- std_prs_ph.prsice %>% select(1, 2, 3, (start_pos_target+best_prs.prsice-1))
best_target.prsice <- best_target.prsice[complete.cases(best_target.prsice), ]
colnames(best_target.prsice) <- c('FID', 'IID', 'pheno', 'std_prs')

#recode phenotype into a factor
best_target.prsice$pheno <- as.factor(best_target.prsice$pheno)

cat("\nPRSice scores for best threshold saved")

setwd(output_data_dir)

#need to write table for best prs
write.table(best_target.prsice, "005_comparison/best_prs_prsice.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nPRSice finished")


# =========================================================================================================================================================================================
### LASSOSUM ###
# =========================================================================================================================================================================================

setwd(output_data_dir)

lasso_thr = 4

# Get TARGET scores
setwd(paste(output_data_dir,'a_3_julia_lasso/target_data',sep=""))

target_prs.lasso <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(target_prs.lasso) <- c('FID', 'IID', 'pheno')
target_prs.lasso <- read.table('IBD_GSA_fin_prs_lasso1.txt', header = T)[,1:3]

for (i in 1:lasso_thr) {
    lasso_file_name = paste(target_prefix, '_prs_lasso', i, '.txt', sep="")
    target_prs.lasso_file <- read.table(lasso_file_name, header = T)
    target_prs.lasso_file <- target_prs.lasso_file %>% select(1,2,5)
    score_col_name <- paste('best.pgs', i, sep="")
    colnames(target_prs.lasso_file) <- c('FID', 'IID', score_col_name)
    target_prs.lasso <- merge(target_prs.lasso, target_prs.lasso_file, by = c("FID", "IID"))    
}

#read in scores
# target.prs.lasso.1 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso1.txt", header = T)
# target.prs.lasso.2 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso2.txt", header = T)
# target.prs.lasso.3 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso3.txt", header = T)
# target.prs.lasso.4 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso4.txt", header = T)

# Get EXTERNAL scores
setwd(paste(output_data_dir,'a_3_julia_lasso/external_data',sep=""))

external_prs.lasso <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(external_prs.lasso) <- c('FID', 'IID')
external_prs.lasso <- read.table('1000G_EUR_fin_prs_lasso1_external.txt', header = T)[,1:2]

for (i in 1:lasso_thr) {
    lasso_file_name = paste(external_prefix, '_prs_lasso', i, '_external.txt', sep="")
    external_prs.lasso_file <- read.table(lasso_file_name, header = T)
    external_prs.lasso_file <- external_prs.lasso_file %>% select(1,2,4)
    score_col_name_ext <- paste('best.pgs', i, sep="")
    colnames(external_prs.lasso_file) <- c('FID', 'IID', score_col_name_ext)
    external_prs.lasso <- merge(external_prs.lasso, external_prs.lasso_file, by = c("FID", "IID"))
} 

# external_prs.lasso.1 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso1_external.txt", header = TRUE)
# external_prs.lasso.1 <- external_prs.lasso.1[complete.cases(external_prs.lasso.1), ]
# external_prs.lasso.2 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso2_external.txt", header = TRUE)
# external_prs.lasso.2 <- external_prs.lasso.2[complete.cases(external_prs.lasso.2), ]
# external_prs.lasso.3 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso3_external.txt", header = TRUE)
# external_prs.lasso.3 <- external_prs.lasso.3[complete.cases(external_prs.lasso.3), ]
# external_prs.lasso.4 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso4_external.txt", header = TRUE)
# external_prs.lasso.4 <- external_prs.lasso.4[complete.cases(external_prs.lasso.4), ]

cat("\nlassoSum target and external scores uploaded into single dataframes")

no_lasso_thr = lasso_thr
start_pos_ext = 3
start_pos_target = 4

ext_scores.lasso <- external_prs.lasso[,start_pos_ext:(start_pos_ext+no_lasso_thr-1)]

ext_mean.lasso <- as.data.frame(colMeans(ext_scores.lasso))
ext_sd.lasso <- as.data.frame(apply(ext_scores.lasso, 2, sd))

# ext_mean.lasso.1 <- mean(external_prs.lasso.1$best.pgs)
# ext_sd.lasso.1 <- sd(external_prs.lasso.1$best.pgs)
# ext_mean.lasso.2 <- mean(external_prs.lasso.2$best.pgs)
# ext_sd.lasso.2 <- sd(external_prs.lasso.2$best.pgs)
# ext_mean.lasso.3 <- mean(external_prs.lasso.3$best.pgs)
# ext_sd.lasso.3 <- sd(external_prs.lasso.3$best.pgs)
# ext_mean.lasso.4 <- mean(external_prs.lasso.4$best.pgs)
# ext_sd.lasso.4 <- sd(external_prs.lasso.4$best.pgs)

colnames(ext_mean.lasso) <- c("ext_mean.prsice")
colnames(ext_sd.lasso) <- c("ext_sd.prsice")

cat("\nlassoSum external scores' mean and standard deviation calculated")
cat("\nlassoSum starting standardization...")

# Before standardization

std_names <- colnames(target_prs.lasso)

for (sc in c(start_pos_target:(start_pos_target+no_lasso_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.lasso)[sc], "_std", sep="")
}

std_prs.lasso <- target_prs.lasso
colnames(std_prs.lasso) <- std_names

for (sc in c(start_pos_target:(start_pos_target+no_lasso_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.lasso[(sc-start_pos_target+1),1]
    stdv <- ext_sd.lasso[(sc-start_pos_target+1),1]
    std_prs.lasso[,sc] <- (target_prs.lasso[,sc] - mean) / stdv
}

# #standardize target
# target.prs.lasso.1$std_prs <- (target.prs.lasso.1$best.pgs - ext_mean.lasso.1) / ext_sd.lasso.1
# target.prs.lasso.2$std_prs <- (target.prs.lasso.2$best.pgs - ext_mean.lasso.2) / ext_sd.lasso.2
# target.prs.lasso.3$std_prs <- (target.prs.lasso.3$best.pgs - ext_mean.lasso.3) / ext_sd.lasso.3
# target.prs.lasso.4$std_prs <- (target.prs.lasso.4$best.pgs - ext_mean.lasso.4) / ext_sd.lasso.4

cat("\nlassoSum standardization complete")

std_prs_ph.lasso <- std_prs.lasso
std_prs_ph.lasso$pheno <- std_prs_ph.lasso$pheno - 1

# #update the phenotype 
# target.prs.lasso.1$pheno <- target.prs.lasso.1$pheno - 1
# target.prs.lasso.2$pheno <- target.prs.lasso.2$pheno - 1
# target.prs.lasso.3$pheno <- target.prs.lasso.3$pheno - 1
# target.prs.lasso.4$pheno <- target.prs.lasso.4$pheno - 1

explained_var.lasso <- c()
for (i in c(start_pos_target:(start_pos_target+no_lasso_thr-1))) {
    
    score_var <- colnames(std_prs_ph.lasso)[i]
    
    log_model.lasso <- glm(std_prs_ph.lasso$pheno ~ pull(std_prs_ph.lasso, score_var), family = binomial(link = "logit"))
    explained_var.lasso[i] <- with(summary(log_model.lasso), 1 - deviance/null.deviance)
}

#obtain the highest explained variance 
best_prs.lasso <- which.max(explained_var.lasso)

#obtain the best prs by indexing from which model had the best explained variance
best_target.lasso <- std_prs_ph.lasso %>% select(1, 2, 3, (start_pos_target+best_prs.lasso-1))
best_target.lasso <- best_target.lasso[complete.cases(best_target.lasso), ]
colnames(best_target.lasso) <- c('FID', 'IID', 'pheno', 'std_prs')

#recode phenotype into a factor
best_target.lasso$pheno <- as.factor(best_target.lasso$pheno)

cat("\nlassoSum scores for best threshold saved")

# #construct logistic regression to see which is the best of the four
# log_model.lasso.1 <- glm(pheno ~ std_prs, data = target.prs.lasso.1, family = binomial(link = "logit"))
# log_model.lasso.2 <- glm(pheno ~ std_prs, data = target.prs.lasso.2, family = binomial(link = "logit"))
# log_model.lasso.3 <- glm(pheno ~ std_prs, data = target.prs.lasso.3, family = binomial(link = "logit"))
# log_model.lasso.4 <- glm(pheno ~ std_prs, data = target.prs.lasso.4, family = binomial(link = "logit"))

# #make a list of the explained variance
# explained_var.lasso <- c(with(summary(log_model.lasso.1), 1 - deviance/null.deviance),
# 						with(summary(log_model.lasso.2), 1 - deviance/null.deviance), 
# 						with(summary(log_model.lasso.3), 1 - deviance/null.deviance),
# 						with(summary(log_model.lasso.4), 1 - deviance/null.deviance))

# #obtain the highest explained variance 
# best_prs.lasso <- which.max(explained_var.lasso)

# #make a list of the prs datasets
# targets.lasso <- list(target.prs.lasso.1, target.prs.lasso.2, target.prs.lasso.3, target.prs.lasso.4)

# #obtain the best prs by indexing from which model had the best explained variance
# best_target.lasso <- as.data.frame(targets.lasso[best_prs.lasso])
# best_target.lasso <- best_target.lasso[complete.cases(best_target.lasso), ]

# #recode phenotype into a factor
# best_target.lasso$pheno <- as.factor(best_target.lasso$pheno)

setwd(output_data_dir)


#need to write table for best prs
write.table(best_target.lasso, "005_comparison/best_prs_lassoSum.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nlassoSum finished")



# =========================================================================================================================================================================================
### LDPRED ###
# =========================================================================================================================================================================================


# =========================================================================================================================================================================================
### BOXPLOTS ###
# =========================================================================================================================================================================================

# PLINK 

# #make boxplot
# ggplot(best_target.plink, aes(y = std_prs, group = pheno, fill = pheno, alpha = 0.5)) +
# labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
# geom_boxplot() + guides (alpha = "none") + theme_bw() + 
# theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
#       axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases"))

# # PRSice

# #make boxplot
# ggplot(best_target.prsice, aes(y = std_prs, group = pheno, fill = pheno, alpha = 0.5)) +
# labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
# geom_boxplot() + guides (alpha = "none") + theme_bw() + 
# theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
#       axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases"))

# # lassoSum

# #make boxplot
# ggplot(best_target.lasso, aes(y = std_prs, group = pheno, fill = pheno, alpha = 0.5)) +
# labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
# geom_boxplot() + guides (alpha = "none") + theme_bw() + 
# theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
#       axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases"))






# =========================================================================================================================================================================================
### Run the ROC ###
# =========================================================================================================================================================================================

# setwd('/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/output_data/005_comparison')
# # best_target.prsice <- read.table("best_prs_prsice.txt", header=T)
# # best_target.lasso <- read.table("best_prs_lassoSum.txt", header=T)

# head(best_target.prsice)
# head(best_target.lasso)


# prsice.roc <- roc(best_target.prsice$pheno, best_target.prsice$std_prs, smooth = F)
# lasso.roc <- roc(best_target.lasso$pheno, best_target.lasso$std_prs, smooth = F)

# #do the bootstrapping 
# prsice.roc.ci <- ci.auc(prsice.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
# lasso.roc.ci <- ci.auc(lasso.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)

# #plot the roc curve
# #to add more curves just add the roc output to the list
# #example: `list(lasso=lasso.roc, PLINK = plink.roc)

# ggroc(list(prsice=prsice.roc, lassoSum=lasso.roc)) + 
# labs(title = "Reciever Operator Curve Across\nConstructed Logisitc Models from PRS Tools ",y = "Sensitivity (True Positive Rate)", x = "Specificty (1 - False Positive Rate)", col = "Tool") +
# theme_bw() + theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12)) +
# geom_abline(intercept = 1, slope = 1)


# #write a function to allow for rounding
# scaleFUN <- function(x) sprintf("%.3f", x)

# #create a dataframe for the auc metric 
# performance_metrics_auc <- data.frame(matrix(ncol = 4, nrow = 0))


# #enter the estimates
# plink_perform_auc_FAKE <- c('PLINK_fake', 0.6734, 0.73241, 0.7876)
# prsice_perform_auc <- c("PRSice",prsice.roc.ci[1],prsice.roc.ci[2],prsice.roc.ci[3])
# lasso_perform_auc <- c("lassoSum",lasso.roc.ci[1],lasso.roc.ci[2],lasso.roc.ci[3])
# ldpred_perform_auc_FAKE <- c('LDpred_fake', 0.6734, 0.73241, 0.7876)

# #combine with the empty dataframe
# performance_metrics_auc <- rbind(performance_metrics_auc, 
#                                  plink_perform_auc_FAKE, prsice_perform_auc, 
#                                  lasso_perform_auc, ldpred_perform_auc_FAKE)
# #rename the columns of the dataframe
# colnames(performance_metrics_auc) <- c("tool", "lower_95", "estimate", "upper_95")

# #make the dataframe columns numeric
# performance_metrics_auc$lower_95 <- as.numeric(performance_metrics_auc$lower_95)
# performance_metrics_auc$estimate <- as.numeric(performance_metrics_auc$estimate)
# performance_metrics_auc$upper_95 <- as.numeric(performance_metrics_auc$upper_95)
# head(performance_metrics_auc)

# #make the r-sqr plot
# ggplot(performance_metrics_auc, aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
# labs(title = bquote("Bootstrapped Estimation of AUC Value Across Tools"), y = bquote(R^2~"(Percentage of Explained Variation)"), x = '', fill = 'Tool') +
# theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.02) +
# theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())+
# geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 



