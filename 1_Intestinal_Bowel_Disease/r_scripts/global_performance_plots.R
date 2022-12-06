## BOXPLOTS, ROC, AUC, R2

print("Entering the R script")

library(dplyr)
library(ggplot2)
library(data.table)
library(ROCR)
library(pROC)
library(pscl)
library(stringr)

# =========================================================================================================================================================================================

# project_dir = snakemake@params[[1]]
project_dir = "/staging/leuven/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/"
output_data_dir = paste(project_dir, "output_data/", sep="")

setwd(paste(output_data_dir,'005_good_comparison_2',sep=""))

# best_target.plink_file <- snakemake@input[[1]]
# best_target.prsice_file <- snakemake@input[[2]]
# best_target.lasso_file <- snakemake@input[[3]]
# # best_target.ldpred_file <- snakemake@input[[4]]

# print(best_target.plink_file)
# print("See if we get the whole path or only the name")

# best_target.plink <- read.table(best_target.plink_file, header=T)
# best_target.prsice <- read.table(best_target.prsice_file, header=T)
# best_target.lasso <- read.table(best_target.lasso_file, header=T)
# # best_target.ldpred <- read.table(best_target.ldpred_file, header=T)

best_target.plink <- read.table("best_prs_plink.txt", header=T)
best_target.prsice <- read.table("best_prs_prsice.txt", header=T)
best_target.lasso <- read.table("best_prs_lassoSum.txt", header=T)
best_target.ldpred <- best_target.lasso

# =========================================================================================================================================================================================
### BOXPLOTS ###
# =========================================================================================================================================================================================

tmp.plink <- best_target.plink
tmp.plink$tool <- "PLINK"
tmp.prsice <- best_target.prsice
tmp.prsice$tool <- "PRSice"
tmp.lasso <- best_target.lasso
tmp.lasso$tool <- "lassoSum"
tmp.ldpred <- best_target.ldpred
tmp.ldpred$tool <- "LDpred"

tmp.all <- rbind(tmp.plink, tmp.prsice, tmp.lasso, tmp.ldpred)

tmp.all$pheno <- as.factor(tmp.all$pheno)

#make boxplot
tmp.all %>%
    arrange(std_prs) %>%
    mutate(tool = factor(tool, levels=c('PLINK', 'PRSice', 'lassoSum', 'LDpred'))) %>%
    ggplot( aes(y = std_prs, group = pheno, fill = pheno, alpha = 0.5)) +
    labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
    geom_boxplot() + guides (alpha = "none") + theme_bw() + 
    theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases")) +
    facet_wrap(~ tool, ncol = 4)



# =========================================================================================================================================================================================
### Run the ROC ###
# =========================================================================================================================================================================================

plink.roc <- roc(best_target.plink$pheno, best_target.plink$std_prs, smooth = F)
prsice.roc <- roc(best_target.prsice$pheno, best_target.prsice$std_prs, smooth = F)
lasso.roc <- roc(best_target.lasso$pheno, best_target.lasso$std_prs, smooth = F)
# ldpred.roc <- roc(best_target.ldpred$pheno, best_target.ldpred$std_prs, smooth = F)

#do the bootstrapping 
plink.roc.ci <- ci.auc(plink.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
prsice.roc.ci <- ci.auc(prsice.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
lasso.roc.ci <- ci.auc(lasso.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
# ldpred.roc.ci <- ci.auc(ldpred.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)


# =====================================================================
## when LDpred results => delete
ldpred.roc <- lasso.roc
ldpred.roc.ci <- lasso.roc.ci
# =====================================================================

#plot the roc curve
#to add more curves just add the roc output to the list
#example: `list(lasso=lasso.roc, PLINK = plink.roc)

ggroc(list(PLINK=plink.roc, PRSice=prsice.roc, lassoSum=lasso.roc, LDpred=ldpred.roc)) + 
  labs(title = "Reciever Operator Curve Across\nConstructed Logisitc Models from PRS Tools ",y = "Sensitivity (True Positive Rate)", x = "Specificty (1 - False Positive Rate)", col = "Tool") +
  theme_bw() + theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12)) +
  geom_abline(intercept = 1, slope = 1)


# =========================================================================================================================================================================================
### Performance metrics: AUC and R2 ###
# =========================================================================================================================================================================================

### AUC ###

#write a function to allow for rounding
scaleFUN <- function(x) sprintf("%.3f", x)

#create a dataframe for the auc metric 
performance_metrics_auc <- data.frame(matrix(ncol = 4, nrow = 0))

#enter the estimates
plink_perform_auc <- c('PLINK', plink.roc.ci[1],plink.roc.ci[2],plink.roc.ci[3])
prsice_perform_auc <- c("PRSice",prsice.roc.ci[1],prsice.roc.ci[2],prsice.roc.ci[3])
lasso_perform_auc <- c("lassoSum",lasso.roc.ci[1],lasso.roc.ci[2],lasso.roc.ci[3])
ldpred_perform_auc <- c('LDpred',ldpred.roc.ci[1],ldpred.roc.ci[2],ldpred.roc.ci[3])

#combine with the empty dataframe
performance_metrics_auc <- rbind(performance_metrics_auc, plink_perform_auc, prsice_perform_auc, lasso_perform_auc, ldpred_perform_auc)
#rename the columns of the dataframe
colnames(performance_metrics_auc) <- c("tool", "lower_95", "estimate", "upper_95")

#make the dataframe columns numeric
performance_metrics_auc$lower_95 <- as.numeric(performance_metrics_auc$lower_95)
performance_metrics_auc$estimate <- as.numeric(performance_metrics_auc$estimate)
performance_metrics_auc$upper_95 <- as.numeric(performance_metrics_auc$upper_95)

#make the auc plot
performance_metrics_auc %>%
    arrange(estimate) %>%
    mutate(tool = factor(tool, levels=c('PLINK', 'PRSice', 'lassoSum', 'LDpred'))) %>%
    ggplot( aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
    labs(title = bquote("Bootstrapped Estimation of AUC Value Across Tools"), y = "AUC", x = '', fill = 'Tool') +
    theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.02) +
    theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 


### R2 ###

#bootstrap
n <- nrow(best_target.plink)
resamples <- 1000

bootstrap_pseudor2_plink <- sapply(1:resamples, function(j) {
  bootstraps <- sample(c(1:n), n, TRUE)
  pR2(glm(best_target.plink$pheno[bootstraps] ~ best_target.plink$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
})

bootstrap_pseudor2_prsice <- sapply(1:resamples, function(j) {
  bootstraps <- sample(c(1:n), n, TRUE)
  pR2(glm(best_target.prsice$pheno[bootstraps] ~ best_target.prsice$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
})

bootstrap_pseudor2_lasso <- sapply(1:resamples, function(j) {
  bootstraps <- sample(c(1:n), n, TRUE)
  pR2(glm(best_target.lasso$pheno[bootstraps] ~ best_target.lasso$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
})

# bootstrap_pseudor2_ldpred <- sapply(1:resamples, function(j) {
#   bootstraps <- sample(c(1:n), n, TRUE)
#   pR2(glm(best_target.ldpred$pheno[bootstraps] ~ best_target.ldpred$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
# })

# =====================================================================
## when LDpred results => delete
bootstrap_pseudor2_ldpred <- bootstrap_pseudor2_lasso
# =====================================================================

#save the upper and lower bounds, as the estimate
lower_95_plink <- mean(bootstrap_pseudor2_plink) - 1.96 * sd(bootstrap_pseudor2_plink)
r_sqr_est_plink <-mean(bootstrap_pseudor2_plink)
upper_95_plink <- mean(bootstrap_pseudor2_plink) + 1.96 * sd(bootstrap_pseudor2_plink)

lower_95_prsice <- mean(bootstrap_pseudor2_prsice) - 1.96 * sd(bootstrap_pseudor2_prsice)
r_sqr_est_prsice <-mean(bootstrap_pseudor2_prsice)
upper_95_prsice <- mean(bootstrap_pseudor2_prsice) + 1.96 * sd(bootstrap_pseudor2_prsice)

lower_95_lasso <- mean(bootstrap_pseudor2_lasso) - 1.96 * sd(bootstrap_pseudor2_lasso)
r_sqr_est_lasso <-mean(bootstrap_pseudor2_lasso)
upper_95_lasso <- mean(bootstrap_pseudor2_lasso) + 1.96 * sd(bootstrap_pseudor2_lasso)

lower_95_ldpred <- mean(bootstrap_pseudor2_ldpred) - 1.96 * sd(bootstrap_pseudor2_ldpred)
r_sqr_est_ldpred <-mean(bootstrap_pseudor2_ldpred)
upper_95_ldpred <- mean(bootstrap_pseudor2_ldpred) + 1.96 * sd(bootstrap_pseudor2_ldpred)


#create an empty dataframe
performance_metrics_rsqr <- data.frame(matrix(ncol = 4, nrow = 0))

#enter the estimates
plink_perform_rsqr <- c('PLINK', lower_95_plink, r_sqr_est_plink ,upper_95_plink)
prsice_perform_rsqr <- c("PRSice",lower_95_prsice, r_sqr_est_prsice ,upper_95_prsice)
lasso_perform_rsqr <- c("lassoSum",lower_95_lasso, r_sqr_est_lasso ,upper_95_lasso)
ldpred_perform_rsqr <- c("LDpred",lower_95_ldpred, r_sqr_est_ldpred ,upper_95_ldpred)


#combine with the empty dataframe
performance_metrics_rsqr <- rbind(performance_metrics_rsqr, plink_perform_rsqr, prsice_perform_rsqr, lasso_perform_rsqr, ldpred_perform_rsqr)
#rename the columns of the dataframe
colnames(performance_metrics_rsqr) <- c("tool", "lower_95", "estimate", "upper_95")

#make the dataframe columns numeric
performance_metrics_rsqr$lower_95 <- as.numeric(performance_metrics_rsqr$lower_95)
performance_metrics_rsqr$estimate <- as.numeric(performance_metrics_rsqr$estimate)
performance_metrics_rsqr$upper_95 <- as.numeric(performance_metrics_rsqr$upper_95)


#make the r-sqr plot
# ggplot(performance_metrics_rsqr, aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
# labs(title = bquote("Bootstrapped Estimation of"~R^2~"Value Across Tools"), y = bquote(R^2~"(Percentage of Explained Variation)"), x = '', fill = 'Tool') +
# theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.005) +
# theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())+
# geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 

#make the r-sqr plot
performance_metrics_rsqr %>%
    arrange(estimate) %>%
    mutate(tool = factor(tool, levels=c('PLINK', 'PRSice', 'lassoSum', 'LDpred'))) %>%
    ggplot(aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
    labs(title = bquote("Bootstrapped Estimation of"~R^2~"Value Across Tools"), y = bquote(R^2~"(Percentage of Explained Variation)"), x = '', fill = 'Tool') +
    theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.005) +
    theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 


