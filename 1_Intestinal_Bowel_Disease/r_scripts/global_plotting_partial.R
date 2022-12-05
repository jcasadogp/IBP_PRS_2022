## Try formula variable

library(dplyr)
library(ggplot2)
library(data.table)
library(ROCR)
library(pROC)
library(pscl)

#no_plink_thr = 
no_prsice_thr = 7
# no_lasso_thr = 
# no_ldpred_thr = 


# Project directory
setwd('/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/')

# Phenotypes file
phenos <- read.table("data/final_phenotypes.txt", header=T)
colnames(phenos) <- c("IID", "IID.1", "phenos")

# =========================================================================================================================================================================================
### PLINK ###
# =========================================================================================================================================================================================


p.threshold <- c(5e-8,1e-5,0.01,0.05,0.1,0.5)

phenotype <- read.table("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/data/final_phenotypes.txt", header=T)
colnames(phenotype) <- c("FID", "IID","pheno") 

pcs <- read.table("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/data/IBD_GSA_fin.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
covariate <- read.table("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/data/IBD_GSA_fin.cov", header=T)

pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))

pheno$pheno <- pheno$pheno - 1

null.model <- glm(pheno~., data=pheno[,!colnames(pheno)%in%c("FID","IID")],  family = binomial(link = "logit"))

summary(null.model)$coeff["PC1",]

#summary(null.model)
null.r2 <- with(summary(null.model), 1 - deviance/null.deviance)

prs.result <- NULL
for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/snake_isaac/IBD_GSA_fin.",i,".profile"), header=T)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns

    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    #model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    model <- glm(pheno~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")],  family = binomial(link = "logit"))
    # model R2 is obtained as 
    #model.r2 <- summary(model)$r.squared
    model.r2 <- with(summary(model), 1 - deviance/null.deviance)
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    #print(summary(model))
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}

 prs <- read.table(paste0("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/snake_isaac/IBD_GSA_fin.",5e-8,".profile"), header = T)

pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))

library(ggplot2)

#recode phenotype into a factor
pheno.prs$pheno <- as.factor(pheno.prs$pheno)

ggplot(pheno.prs, aes(y = (SCORE - mean(SCORE)) / sd(SCORE), group = pheno, fill = pheno, alpha = 0.5)) +
labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
geom_boxplot() + guides (alpha = "none") + theme_bw() + 
theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases"))





p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Read in the phenotype file 
phenotype <- read.table("EUR.height", header=T)
# Read in the PCs
pcs <- read.table("EUR.eigenvec", header=F)
# The default output from plink does not include a header
# To make things simple, we will add the appropriate headers
# (1:6 because there are 6 PCs)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
# Read in the covariates (here, it is sex)
covariate <- read.table("EUR.cov", header=T)
# Now merge the files
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
# We can then calculate the null model (model with PRS) using a linear regression 
# (as height is quantitative)
null.model <- lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
# And the R2 of the null model is 
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL
for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("EUR.",i,".profile"), header=T)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    # model R2 is obtained as 
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}
# Best result is:
prs.result[which.max(prs.result$R2),]
q() # exit R

phenotype <- read.table("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/data/final_phenotypes.txt", header=T)
colnames(phenotype) <- c("FID", "IID","pheno") 
pcs <- read.table("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/data/IBD_GSA_fin.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 

scores <- read.table("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/output_data/002_prsice_2/target_data/IBD_GSA_fin.all_score", header=T)



#pheno <- merge(phenotype, covariate, by=c("FID", "IID"))
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
pheno$pheno <- pheno$pheno - 1


new_scores <- merge(phenotype, scores, by=c("FID", "IID"))
new_scores$pheno <- new_scores$pheno-1


prs.result <- NULL

new_scores=new_scores[,!colnames(new_scores)%in%c("FID","IID")]

for(i in 3:ncol(scores) - 2){
    #model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    model <- glm(new_scores$pheno~ new_scores[,i+1],  family = binomial(link = "logit"))
    # model R2 is obtained as 
    #model.r2 <- summary(model)$r.squared
    model.r2 <- with(summary(model), 1 - deviance/null.deviance)
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    #print(summary(model))
    prs.coef <- summary(model)$coeff[]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}



new_scores[,!colnames(new_scores)%in%c("FID","IID")]










# =========================================================================================================================================================================================
### PRSICE ###
# =========================================================================================================================================================================================

target_prs.prsice <- read.table("output_data/a_2_julia_prsice/target_data/IBD_GSA_fin.all_score", header = T)
target_prs.prsice <- merge(target_prs.prsice, phenos, by = "IID")

# >>> FILL WITH THE POSITION IN WHICH THE SCORES START
starting_scores_pos = 3

external_prs.prsice <- read.table("output_data/a_2_julia_prsice/external_data/1000G_EUR_fin.all_score", header = T)


ext_scores.prsice <- external_prs.prsice[,starting_scores_pos:(starting_scores_pos+no_prsice_thr-1)]

ext_mean.prsice <- as.data.frame(colMeans(ext_scores.prsice))
ext_sd.prsice <- as.data.frame(apply(ext_scores.prsice, 2, sd))

colnames(ext_mean.prsice) <- c("ext_mean.prsice")
colnames(ext_sd.prsice) <- c("ext_sd.prsice")

#ext_mean.prsice
#ext_sd.prsice

# ..... standardization .....

std_names <- colnames(target_prs.prsice)

for (sc in c(starting_scores_pos:(starting_scores_pos+no_prsice_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.prsice)[sc], "_std", sep="")
}

std_prs.prsice <- target_prs.prsice
colnames(std_prs.prsice) <- std_names

for (sc in c(starting_scores_pos:(starting_scores_pos+no_prsice_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.prsice[(sc-starting_scores_pos+1),1]
    stdv <- ext_sd.prsice[(sc-starting_scores_pos+1),1]
    std_prs.prsice[,sc] <- (target_prs.prsice[,sc] - mean) / stdv
}

# ..... update the phenotype values .....
std_prs_ph.prsice <- std_prs.prsice
std_prs_ph.prsice$phenos <- std_prs_ph.prsice$phenos - 1


plot(std_prs_ph.prsice[,3], std_prs_ph.prsice$phenos)

# construct logistic regression to see which is the best of the four
log_model.prsice.1 <- glm(phenos ~ Pt_5e.08_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))
log_model.prsice.2 <- glm(phenos ~ Pt_1e.05_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))
log_model.prsice.3 <- glm(phenos ~ Pt_0.01_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))
log_model.prsice.4 <- glm(phenos ~ Pt_0.05_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))
log_model.prsice.5 <- glm(phenos ~ Pt_0.1_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))
log_model.prsice.6 <- glm(phenos ~ Pt_0.5_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))
log_model.prsice.7 <- glm(phenos ~ Pt_1_std, data = std_prs_ph.prsice, family = binomial(link = "logit"))

#make a list of the explained variance
explained_var.prsice <- c(with(summary(log_model.prsice.1), 1 - deviance/null.deviance), 
                          with(summary(log_model.prsice.2), 1 - deviance/null.deviance), 
                          with(summary(log_model.prsice.3), 1 - deviance/null.deviance), 
                          with(summary(log_model.prsice.4), 1 - deviance/null.deviance), 
                          with(summary(log_model.prsice.5), 1 - deviance/null.deviance), 
                          with(summary(log_model.prsice.6), 1 - deviance/null.deviance), 
                          with(summary(log_model.prsice.7), 1 - deviance/null.deviance))

#obtain the highest explained variance 
best_prs.prsice <- which.max(explained_var.prsice)

#obtain the best prs by indexing from which model had the best explained variance
best_target.prsice <- std_prs_ph.prsice %>% select(1, 2, (starting_scores_pos+best_prs.prsice-1), (starting_scores_pos+no_prsice_thr+1)) #Pt_5e.08
best_target.prsice <- best_target.prsice[complete.cases(best_target.prsice), ]
colnames(best_target.prsice) <- c('IID', 'FID', 'std_prs', 'phenos')
head(best_target.prsice)

#recode phenotype into a factor
best_target.prsice$phenos <- as.factor(best_target.prsice$phenos)

#make boxplot
ggplot(best_target.prsice, aes(y = std_prs, group = phenos, fill = phenos, alpha = 0.5)) +
labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
geom_boxplot() + guides (alpha = "none") + theme_bw() + 
theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases"))

setwd("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/output_data/")

#need to write table for best prs
write.table(best_target.prsice, "005_comparison/best_prs_prsice.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)



# =========================================================================================================================================================================================
### LASSOSUM ###
# =========================================================================================================================================================================================

setwd("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/output_data/")

#read in scores
target.prs.lasso.1 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso1.txt", header = T)
target.prs.lasso.2 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso2.txt", header = T)
target.prs.lasso.3 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso3.txt", header = T)
target.prs.lasso.4 <- read.table("a_3_julia_lasso/target_data/IBD_GSA_fin_prs_lasso4.txt", header = T)

external_prs.lasso.1 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso1_external.txt", header = TRUE)
external_prs.lasso.1 <- external_prs.lasso.1[complete.cases(external_prs.lasso.1), ]
external_prs.lasso.2 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso2_external.txt", header = TRUE)
external_prs.lasso.2 <- external_prs.lasso.2[complete.cases(external_prs.lasso.2), ]
external_prs.lasso.3 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso3_external.txt", header = TRUE)
external_prs.lasso.3 <- external_prs.lasso.3[complete.cases(external_prs.lasso.3), ]
external_prs.lasso.4 <- read.table("a_3_julia_lasso/external_data/1000G_EUR_fin_prs_lasso4_external.txt", header = TRUE)
external_prs.lasso.4 <- external_prs.lasso.4[complete.cases(external_prs.lasso.4), ]

ext_mean.lasso.1 <- mean(external_prs.lasso.1$best.pgs)
ext_sd.lasso.1 <- sd(external_prs.lasso.1$best.pgs)
ext_mean.lasso.2 <- mean(external_prs.lasso.2$best.pgs)
ext_sd.lasso.2 <- sd(external_prs.lasso.2$best.pgs)
ext_mean.lasso.3 <- mean(external_prs.lasso.3$best.pgs)
ext_sd.lasso.3 <- sd(external_prs.lasso.3$best.pgs)
ext_mean.lasso.4 <- mean(external_prs.lasso.4$best.pgs)
ext_sd.lasso.4 <- sd(external_prs.lasso.4$best.pgs)

#standardize target
target.prs.lasso.1$std_prs <- (target.prs.lasso.1$best.pgs - ext_mean.lasso.1) / ext_sd.lasso.1
target.prs.lasso.2$std_prs <- (target.prs.lasso.2$best.pgs - ext_mean.lasso.2) / ext_sd.lasso.2
target.prs.lasso.3$std_prs <- (target.prs.lasso.3$best.pgs - ext_mean.lasso.3) / ext_sd.lasso.3
target.prs.lasso.4$std_prs <- (target.prs.lasso.4$best.pgs - ext_mean.lasso.4) / ext_sd.lasso.4

#update the phenotype 
target.prs.lasso.1$pheno <- target.prs.lasso.1$pheno - 1
target.prs.lasso.2$pheno <- target.prs.lasso.2$pheno - 1
target.prs.lasso.3$pheno <- target.prs.lasso.3$pheno - 1
target.prs.lasso.4$pheno <- target.prs.lasso.4$pheno - 1

#construct logistic regression to see which is the best of the four
log_model.lasso.1 <- glm(pheno ~ std_prs, data = target.prs.lasso.1, family = binomial(link = "logit"))
log_model.lasso.2 <- glm(pheno ~ std_prs, data = target.prs.lasso.2, family = binomial(link = "logit"))
log_model.lasso.3 <- glm(pheno ~ std_prs, data = target.prs.lasso.3, family = binomial(link = "logit"))
log_model.lasso.4 <- glm(pheno ~ std_prs, data = target.prs.lasso.4, family = binomial(link = "logit"))

#make a list of the explained variance
explained_var.lasso <- c(with(summary(log_model.lasso.1), 1 - deviance/null.deviance),
						with(summary(log_model.lasso.2), 1 - deviance/null.deviance), 
						with(summary(log_model.lasso.3), 1 - deviance/null.deviance),
						with(summary(log_model.lasso.4), 1 - deviance/null.deviance))

#obtain the highest explained variance 
best_prs.lasso <- which.max(explained_var.lasso)

#make a list of the prs datasets
targets.lasso <- list(target.prs.lasso.1, target.prs.lasso.2, target.prs.lasso.3, target.prs.lasso.4)

#obtain the best prs by indexing from which model had the best explained variance
best_target.lasso <- as.data.frame(targets.lasso[best_prs.lasso])
best_target.lasso <- best_target.lasso[complete.cases(best_target.lasso), ]

#recode phenotype into a factor
best_target.lasso$pheno <- as.factor(best_target.lasso$pheno)

#make boxplot
ggplot(best_target.lasso, aes(y = std_prs, group = pheno, fill = pheno, alpha = 0.5)) +
labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
geom_boxplot() + guides (alpha = "none") + theme_bw() + 
theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases"))

setwd("/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/output_data/")

#need to write table for best prs
write.table(best_target.lasso, "005_comparison/best_prs_lassoSum.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)



# =========================================================================================================================================================================================
### LDPRED ###
# =========================================================================================================================================================================================












# =========================================================================================================================================================================================
### Run the ROC ###
# =========================================================================================================================================================================================

setwd('/lustre1/project/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/output_data/005_comparison')
# best_target.prsice <- read.table("best_prs_prsice.txt", header=T)
# best_target.lasso <- read.table("best_prs_lassoSum.txt", header=T)

head(best_target.prsice)
head(best_target.lasso)


prsice.roc <- roc(best_target.prsice$phenos, best_target.prsice$std_prs, smooth = F)
lasso.roc <- roc(best_target.lasso$pheno, best_target.lasso$std_prs, smooth = F)

#do the bootstrapping 
prsice.roc.ci <- ci.auc(prsice.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
lasso.roc.ci <- ci.auc(lasso.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)

#plot the roc curve
#to add more curves just add the roc output to the list
#example: `list(lasso=lasso.roc, PLINK = plink.roc)

ggroc(list(prsice=prsice.roc, lassoSum=lasso.roc)) + 
labs(title = "Reciever Operator Curve Across\nConstructed Logisitc Models from PRS Tools ",y = "Sensitivity (True Positive Rate)", x = "Specificty (1 - False Positive Rate)", col = "Tool") +
theme_bw() + theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12)) +
geom_abline(intercept = 1, slope = 1)


#write a function to allow for rounding
scaleFUN <- function(x) sprintf("%.3f", x)

#create a dataframe for the auc metric 
performance_metrics_auc <- data.frame(matrix(ncol = 4, nrow = 0))


#enter the estimates
plink_perform_auc_FAKE <- c('PLINK_fake', 0.6734, 0.73241, 0.7876)
prsice_perform_auc <- c("PRSice",prsice.roc.ci[1],prsice.roc.ci[2],prsice.roc.ci[3])
lasso_perform_auc <- c("lassoSum",lasso.roc.ci[1],lasso.roc.ci[2],lasso.roc.ci[3])
ldpred_perform_auc_FAKE <- c('LDpred_fake', 0.6734, 0.73241, 0.7876)

#combine with the empty dataframe
performance_metrics_auc <- rbind(performance_metrics_auc, 
                                 plink_perform_auc_FAKE, prsice_perform_auc, 
                                 lasso_perform_auc, ldpred_perform_auc_FAKE)
#rename the columns of the dataframe
colnames(performance_metrics_auc) <- c("tool", "lower_95", "estimate", "upper_95")

#make the dataframe columns numeric
performance_metrics_auc$lower_95 <- as.numeric(performance_metrics_auc$lower_95)
performance_metrics_auc$estimate <- as.numeric(performance_metrics_auc$estimate)
performance_metrics_auc$upper_95 <- as.numeric(performance_metrics_auc$upper_95)
head(performance_metrics_auc)

#make the r-sqr plot
ggplot(performance_metrics_auc, aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
labs(title = bquote("Bootstrapped Estimation of AUC Value Across Tools"), y = bquote(R^2~"(Percentage of Explained Variation)"), x = '', fill = 'Tool') +
theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.02) +
theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 



