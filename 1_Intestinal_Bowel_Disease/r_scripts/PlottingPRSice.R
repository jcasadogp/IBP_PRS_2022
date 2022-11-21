setwd('/staging/leuven/stg_00092/IBP_PRSproject/')

library(ROCR)
library(ggplot2)
library(dplyr)

pheno1 <- read.table('/staging/leuven/stg_00092/IBP_PRSproject/data/final_phenotypes.txt', head=TRUE)
prsice_best <- read.table('/staging/leuven/stg_00092/IBP_PRSproject/data/002_prsice_results/final_IBP.best', head=TRUE)

pheno2 <- pheno1 %>% select("IID", "IBDvsCON")
colnames(pheno2) <- c("IID", "phenotype")
# head(pheno2)

prsice_best <- na.omit(prsice_best)
# colnames(prsice_best)
# head(prsice_best)

prsice_prs <- merge(prsice_best, pheno2, by = "IID", all.x = TRUE)
prsice_prs[prsice_prs$IID %in% c('9263', '9535', '5924'),]

nrow(prsice_prs[(is.na(prsice_prs$pheno)),])

controls <- prsice_prs %>% filter(phenotype == 1)
cases <- prsice_prs %>% filter(phenotype == 2)

mean_cases <- mean(na.omit(cases$PRS)) 
iqr_cases <- IQR(na.omit(cases$PRS))

mean_controls <- mean(na.omit(controls$PRS))
iqr_controls <- IQR(na.omit(controls$PRS))

ybot <- ifelse(mean_controls - (iqr_controls * 1.5) < mean_cases - (iqr_cases * 1.5), iqr_controls*1.5, iqr_cases*1.5)
ytop <- ifelse(mean_controls + (iqr_controls * 1.5) > mean_cases + (iqr_cases * 1.5), iqr_controls*1.5, iqr_cases*1.5)

prsice_prs$phenotype <- as.factor(prsice_prs$phenotype)
# head(prsice_prs)

