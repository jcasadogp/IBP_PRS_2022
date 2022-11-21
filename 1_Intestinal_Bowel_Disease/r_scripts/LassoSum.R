#this is the script for running a LassoSum analysis 


#the install command needs to be run every new interactive session
library(devtools)
install_github("tshmak/lassosum")

library(dplyr)
library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)

#invoke 2 threads
cl <- makeCluster(2)

#read in the files
sum.stat <- snakemake@input[[1]]
bfile <- snakemake@params[[1]]

# Read in and process the covariates 
#we do not have any covariates
covariate <- fread(snakemake@input[[2]]) 
pcs <- fread(snakemake@input[[3]]) 

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
#prefix <- "/lustre1/project/stg_00092/IBP_PRSproject/lassoSumOutput/final_IBD"

# Read in the target phenotype file
target.pheno <- fread(snakemake@input[[4]])

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
pdf(snakemake@output[[1]]) 
#Model validation refers to the process of confirming that the model actually achieves its intended purpose
plot(target.res)#, ylab = "Validation (correlation)")
# Close the pdf file
dev.off() 
#the best s and lambda will correspond to the highest validation
#value 

#save the result from the validation
write.table(best_s, quote = FALSE, snakemake@output[[2]])
write.table(best_lambda, quote = FALSE, snakemake@output[[3]])
write.table(results_pgs, quote = FALSE, snakemake@output[[4]])
write.table(validation_table, quote = FALSE, snakemake@output[[5]])


# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
write.table(r2, quote = FALSE, snakemake@output[[6]])
