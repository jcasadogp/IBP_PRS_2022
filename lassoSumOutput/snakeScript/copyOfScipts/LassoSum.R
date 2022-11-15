#this is the script for running a LassoSum analysis 


#the install command needs to be run every new interactive session
library(devtools)
#install_github("tshmak/lassosum")

library(dplyr)
library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)

#invoke 2 threads
cl <- makeCluster(2)

#read in the files
sum.stat <- "/lustre1/project/stg_00092/IBP_PRSproject/IBD_basefile_hg38_PRS.txt"
bfile <- "/lustre1/project/stg_00092/IBP_PRSproject/final_IBP"

# Read in and process the covariates 
#we do not have any covariates
covariate <- fread("/lustre1/project/stg_00092/IBP_PRSproject/final_IBP.cov") 
pcs <- fread("/lustre1/project/stg_00092/IBP_PRSproject/final_IBP.eigenvec") 

#index pcs so only have 8 columns 
pcs <- pcs[,1:8]

#then set the column names 
setnames(pcs, colnames(pcs), c("FID","IID", paste0("PC",1:6)))
covariate$FID <- as.character(covariate$FID)
covariate$IID <- as.character(covariate$IID)

# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs, by = c("FID", "IID"))

#save reference to human 38 genome
ld.file <- "EUR.hg38"

# output prefix
prefix <- "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/final_IBD"

# Read in the target phenotype file
target.pheno <- fread("/lustre1/project/stg_00092/IBP_PRSproject/final_phenotypes.txt")

#our names are messed up, rename
colnames(target.pheno)[1] = "FID"

# Read in the summary statistics
ss <- fread(sum.stat)

########modifying our data format
##need to add a column of N the number of variants
#used in calculating the effect size estimates
#for now I guess just use the number of variants 
ss$N = nrow(ss)

#need to add a column for OR
#need to take the 'effect' column and raise it
#to the exponential 
ss$OR = exp(ss$Effect)

#########
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P.value == 0]

# Transform the P-values into correlation
cor <- p2cor(p = ss$P.value,
        n = ss$N, 
        sign = log(ss$OR)
        )

#obtain fam file
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline( 
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$Allele2,
    A2 = ss$Allele1,
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file, 
    cluster=cl
)
# we save the validation output
# we use validation rather than pseudovalidation, 
# because we want to allow for covariates 
target.res <- validate(out, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov))

#save the best best s, lambda and pgs 
best_s <- target.res$best.s
best_lambda <- target.res$best.lambda
results_pgs <- target.res$results.table

#save the validation table from the results
validation_table <- target.res$validation.result

#save the validation plot
pdf("/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/validation_plot.pdf") 
#Model validation refers to the process of confirming that the model actually achieves its intended purpose
plot(target.res)#, ylab = "Validation (correlation)")
# Close the pdf file
dev.off() 
#the best s and lambda will correspond to the highest validation
#value 

#save the result from the validation
write.table(best_s, quote = FALSE, "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/best_s_val.txt")
write.table(best_lambda, quote = FALSE, "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/best_lambda_val.txt")
write.table(results_pgs, quote = FALSE, "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/pgs_results.txt")
write.table(validation_table, quote = FALSE, "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/validation_results.txt")


# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
write.table(r2, quote = FALSE, "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/final_IBD_max_rsqr.txt")



#############################################################################################
#############################1000 Genomes####################################################
#############################################################################################
#NOW run the same but analysis but for the 1000 genome project

sum.stat.1kg <- "/lustre1/project/stg_00092/IBP_PRSproject/IBD_basefile_hg38_PRS.txt"
bfile.1kg <- "/lustre1/project/stg_00092/IBP_PRSproject/PCA2/1kG_MDS9"
race.file <- "/lustre1/project/stg_00092/IBP_PRSproject/PCA2/race_1kG14.txt"

# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well

#save reference to human 38 genome
ld.file.1kg <- "EUR.hg38"

# output prefix
prefix.1kg <- "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/1kG_results"

# Read in the summary statistics
ss.1kg <- fread(sum.stat.1kg)

#read in race file for 1kG
race <- fread(race.file, header = F)
colnames(race) <- c("FID","IID",'ancestry')
eur <- race %>% filter(ancestry == 'EUR')
eur <- eur[,1:2]

########modifying our data format
##need to add a column of N the number of variants
#used in calculating the effect size estimates
#for now I guess just use the number of variants 
ss.1kg$N = nrow(ss.1kg)

#need to add a column for OR
#need to take the 'effect' column and raise it
#to the exponential 
ss.1kg$OR = exp(ss.1kg$Effect)

#########


# Remove P-value = 0, which causes problem in the transformation
ss.1kg <- ss.1kg[!P.value == 0]

# Transform the P-values into correlation
cor.1kg <- p2cor(p = ss.1kg$P.value,
        n = ss.1kg$N, 
        sign = log(ss.1kg$OR)
        )

#obtain fam file
fam.1kg <- fread(paste0(bfile.1kg, ".fam"))
fam.1kg[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]


# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing

out.1kg <- lassosum.pipeline( 
    cor = cor.1kg,
    chr = ss.1kg$CHR,
    pos = ss.1kg$BP,
    A1 = ss.1kg$Allele2,
    A2 = ss.1kg$Allele1,
    keep.test = eur,
    ref.bfile = bfile.1kg,
    test.bfile = bfile.1kg,
    LDblocks = ld.file.1kg, 
    cluster=cl
)
# we have to extract the PGS associated with the best s and best lambda

x = 1
y = 1
for(val in out$s)
{
  if(val == target.res$best.s)
  {
    print(val)     
    for (v in out$lambda)
    {
      if(v == target.res$best.lambda)
      {
        print(v)
        prs.target = c(x,y)
      }
      y = y + 1
    }
  }
  x = x + 1
}

#save the 1kg pgs 
results_pgs.1kg <- out.1kg$pgs[[prs.target[1]]][,prs.target[2]]

#save the result from the validation
write.table(results_pgs.1kg, quote = FALSE, "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/pgs_results.1kg.txt")


